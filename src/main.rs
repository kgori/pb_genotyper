#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(dead_code)]

use anyhow::Error;
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::fs::File;
use std::path::PathBuf;

mod cigar;
mod regions;

use cigar::{check_cigar_overlap, get_read_position, ReadPosition, get_insertion_at_position};

fn get_chrom_names(bamfile: &std::path::Path) -> Result<Vec<String>, Error> {
    let bam = Reader::from_path(bamfile)?;
    let header = bam.header();
    let chroms = header.target_names();
    let chroms = chroms
        .iter()
        .map(|x| String::from_utf8(x.to_vec()))
        .collect::<Result<Vec<_>, _>>();
    match chroms {
        Ok(chroms) => Ok(chroms),
        Err(e) => Err(Error::msg(format!("Error reading chromosome names: {}", e))),
    }
}

#[derive(Debug, Clone)]
struct Variant {
    chrom: String,
    pos: i64,
    ref_base: String,
    alt_base: String,
    is_germline: bool,
    is_indel: bool,
}

impl Variant {
    fn is_insertion(&self) -> bool {
        self.ref_base.len() == 1 && self.alt_base.len() > 1 &&
            self.is_indel &&
            self.ref_base.as_bytes()[0] == self.alt_base.as_bytes()[0]
    }

    fn is_deletion(&self) -> bool {
        self.ref_base.len() > 1 && self.alt_base.len() == 1 &&
            self.is_indel &&
            self.ref_base.as_bytes()[0] == self.alt_base.as_bytes()[0]
    }
}

fn open_tsv(file_path: &str) -> Result<csv::Reader<File>, Error> {
    let file = File::open(file_path)?;
    let reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(file);
    Ok(reader)
}

// Loads variants from a TSV file. Variant positions are 1-based.
fn load_variants(file_path: &str) -> Result<Vec<Variant>, Error> {
    let mut reader = open_tsv(file_path)?;
    reader.records().map(|record| {
        let record = record?;
        let chrom = record.get(0).ok_or_else(|| Error::msg("missing chrom"))?;
        let pos = record.get(1).ok_or_else(|| Error::msg("missing pos"))?;
        let ref_base = record.get(2).ok_or_else(|| Error::msg("missing ref_base"))?;
        let alt_base = record.get(3).ok_or_else(|| Error::msg("missing alt_base"))?;
        let is_germline = record.get(4).ok_or_else(|| Error::msg("missing is_germline"))?;
        let is_indel = record.get(5).ok_or_else(|| Error::msg("missing is_indel"))?;
        Ok(Variant {
            chrom: chrom.to_string(),
            pos: pos.parse()?,
            ref_base: ref_base.to_string(),
            alt_base: alt_base.to_string(),
            is_germline: is_germline.to_lowercase().parse()?,
            is_indel: is_indel.to_lowercase().parse()?,
        })
    }).collect()
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
}

// Loads regions from a TSV file. Region positions are 1-based, inclusive ranges,
// meaning that both the start and end positions are included in the region.
fn load_regions(file_path: &str) -> Result<Vec<Region>, Error> {
    let mut reader = open_tsv(file_path)?;
    reader.records().map(|record| {
        let record = record?;
        let chrom = record.get(0).ok_or_else(|| Error::msg("missing chrom"))?;
        let start = record.get(1).ok_or_else(|| Error::msg("missing start"))?;
        let end = record.get(2).ok_or_else(|| Error::msg("missing end"))?;
        let name = record.get(3).ok_or_else(|| Error::msg("missing name"))?;
        Ok(Region {
            chrom: chrom.to_string(),
            start: start.parse()?,
            end: end.parse()?,
            name: name.to_string(),
        })
    }).collect()
}


fn get_region_with_name(regions: &Vec<Region>, name: &str) -> Option<Region> {
    let found = regions.iter().find(|r| r.name == name);
    match found {
        Some(region) => Some(region.clone()),
        None => None,
    }
}


fn fetch_bam_reads_from_region(bam: &mut IndexedReader, region: &Region) -> Result<(), Error> {
    bam.fetch((region.chrom.as_str(), region.start, region.end))?;
    Ok(())
}


fn get_variants_from_region(variants: &Vec<Variant>, region: &Region) -> Vec<Variant> {
    variants.iter()
        .filter(|v| v.chrom == region.chrom && v.pos >= region.start && v.pos <= region.end).cloned().collect()
}

fn read_overlaps_region(record: &Record, region: &Region) -> bool {
    let overlap = check_cigar_overlap(record, (region.start - 1).into(), region.end.into());
    overlap
}


fn main() -> Result<(), Error> {
    let variants_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/code/pb_explore/data/variants.tsv";
    let variants = load_variants(variants_path)?;
    println!("{:?}", variants[0]);
    
    let regions_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/code/pb_explore/data/ht_regions.tsv";
    let mut regions = load_regions(regions_path)?;
    regions.sort_by(|a, b| a.name.cmp(&b.name));
    println!("{:?}", regions[0]);

    let my_region = get_region_with_name(&regions, "F").unwrap();

    let my_variants = get_variants_from_region(&variants, &my_region);
    println!("{:?}", my_variants[0]);

    let bam_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/alignments/2169Tb.1.hifi_reads.aligned.bam";
    let mut bam = IndexedReader::from_path(bam_path).unwrap();
    let mut record = Record::new();

    fetch_bam_reads_from_region(&mut bam, &my_region)?;

    while let Some(result) = bam.read(&mut record) {
        if record.qname() == b"m84093_240411_112803_s2/223871912/ccs" {
            break;
        }
    }
    
    println!("record = {:?}", record);
    println!("record.qname = {:?}", String::from_utf8(record.qname().to_vec()));
    println!("pos = {:?}", record.pos());
    println!("is supplementary? = {:?}", record.is_supplementary());

    assert!(read_overlaps_region(&record, &my_region));
    for var in &my_variants {
        if var.pos >= record.pos() && var.pos <= (record.pos() + record.seq_len() as i64) {
            let result = genotype_variant(&record, var);
            println!("variant = {:?}", var);
            println!("result = {:?}", result);
        }
    }
    Ok(())
}

#[derive(Debug)]
enum GenotypeResult {
    Reference(String),
    Alternative(String),
    ThirdAllele(String),
    OutOfRange,
    Deleted,
}

fn genotype_insertion(record: &Record, variant: &Variant) -> Result<GenotypeResult, Error> {
    let insertion = get_insertion_at_position(record, (variant.pos - 1).into());
    match insertion {
        Some(insertion) => {
            let seq = String::from_utf8(insertion.seq.clone())?;
            if seq == variant.alt_base {
                Ok(GenotypeResult::Alternative(seq))
            }
            else {
                Ok(GenotypeResult::ThirdAllele(seq))
            }
        },
        None => Ok(GenotypeResult::Reference(variant.ref_base.clone())),
    }
}

fn genotype_snp(record: &Record, variant: &Variant) -> Result<GenotypeResult, Error> {
    let read_position = get_read_position(record, (variant.pos - 1).into());
    let ref_base = variant.ref_base.chars().next().unwrap();
    let alt_base = variant.alt_base.chars().next().unwrap();
    match read_position {
        ReadPosition::Match(i) => {
            let base = record.seq()[i] as char;
            if base == ref_base {
                Ok(GenotypeResult::Reference(base.to_string()))
            }
            else if base == alt_base {
                Ok(GenotypeResult::Alternative(base.to_string()))
            }
            else {
                Ok(GenotypeResult::ThirdAllele(base.to_string()))
            }
        },
        ReadPosition::Deletion => Ok(GenotypeResult::Deleted),
        ReadPosition::NotOverlapped => Ok(GenotypeResult::OutOfRange),
    }
}

fn genotype_deletion(record: &Record, variant: &Variant) -> Result<GenotypeResult, Error> {
    let pos = (variant.pos - 1) as i64;
    // First check that the first base of the deletion is present in the read
    let read_position = get_read_position(record, pos);
    match read_position {
        ReadPosition::Match(i) => {
            let base = record.seq()[i] as char;
            if base != variant.ref_base.chars().next().unwrap() {
                return Err(Error::msg("Read is ambiguous at first position of deletion"));
            }
        },
        _ => return Err(Error::msg("First base of indel not found in read")),
    };

    // Then check the next n bases are deleted
    let deletion_size = variant.ref_base.len() - variant.alt_base.len();
    let mut deletions_found = 0;
    for i in 1..=deletion_size {
        let i: i64 = i as i64;
        let read_position = get_read_position(record, pos + i);
        match read_position {
            ReadPosition::Deletion => deletions_found += 1,
            ReadPosition::Match(_) => {},
            _ => return Err(Error::msg("Read is ambiguous at deletion position")),
        }
    }
    if deletions_found == deletion_size {
        Ok(GenotypeResult::Alternative(variant.alt_base.clone()))
    }
    else {
        Ok(GenotypeResult::Reference(variant.ref_base.clone()))
    }
}


fn genotype_variant(record: &Record, variant: &Variant) -> Result<GenotypeResult, Error> {
    if variant.pos <= record.pos() || variant.pos > (record.pos() + record.seq_len() as i64) {
        return Ok(GenotypeResult::OutOfRange);
    }
    
    if variant.is_indel {
        if variant.is_insertion() {
            return genotype_insertion(record, variant);
        } else if variant.is_deletion() {
            return genotype_deletion(record, variant);
        } else {
            return Err(Error::msg("Indels with multiple bases in both ref and alt are not supported"));
        }
    }
    else {
        return genotype_snp(record, variant);
    }
}
