#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(dead_code)]
#![allow(unused_labels)]

use anyhow::Error;
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::fs::File;
use std::path::PathBuf;
use std::fmt;
use std::io::Write;
use std::collections::{HashMap, HashSet};

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

#[derive(Debug, Clone, Hash, Eq, PartialEq, Ord, PartialOrd)]
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

fn generate_variant_scores(variants: &Vec<Variant>) -> HashMap<Variant, i64> {
    let mut scores = HashMap::new();
    for (i, var) in variants.iter().enumerate() {
        scores.insert(var.clone(), i as i64);
    }
    scores
}

fn generate_sort_keys(read_names: &Vec<String>, results: &HashMap<String, HashMap<Variant, GenotypeResult>>, variant_scores: &HashMap<Variant, i64>) -> HashMap<String, i64> {
    let mut sorter = HashMap::new();
    for read_name in read_names {
        let entry = results.get(read_name);
        match entry {
            Some(entry) => {
                if entry.is_empty() {
                    sorter.insert(read_name.clone(), i64::MAX);
                }
                else {
                    let score = entry.keys().map(|k| variant_scores.get(k).unwrap_or(&i64::MAX)).min().unwrap();
                    sorter.insert(read_name.clone(), *score);
                }
            },
            None => {
                sorter.insert(read_name.clone(), i64::MAX);
            },
        }
    }
    sorter
}

fn main() -> Result<(), Error> {
    let variants_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/code/pb_explore/data/variants.tsv";
    let variants = load_variants(variants_path)?;
    
    let regions_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/code/pb_explore/data/ht_regions.tsv";
    let mut regions = load_regions(regions_path)?;
    regions.sort_by(|a, b| a.name.cmp(&b.name));

    let bam_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/alignments/2169Tb.1.hifi_reads.aligned.bam";
    let mut bam = IndexedReader::from_path(bam_path).unwrap();
    let mut record = Record::new();

    let mut results = HashMap::<String, HashMap<Variant, GenotypeResult>>::new();

    let mut all_read_names: Vec<String> = Vec::new();
    let mut seen_set: HashSet<String> = HashSet::new();
    
    for region in &regions {
        fetch_bam_reads_from_region(&mut bam, region)?;
        let my_variants = get_variants_from_region(&variants, region);
        while let Some(result) = bam.read(&mut record) {
            let key = String::from_utf8(record.qname().to_vec()).unwrap();
            
            if !seen_set.contains(&key) {
                all_read_names.push(key.clone());
                seen_set.insert(key.clone());
            }
            if !results.contains_key(&key) {
                results.insert(key.clone(), HashMap::new());
            }
            if !read_overlaps_region(&record, region) {
                continue;
            }
            for var in &my_variants {
                if var.pos >= record.pos() && var.pos <= (record.pos() + record.seq_len() as i64) {
                    let result = genotype_variant(&record, var)?;
                    match result {
                        GenotypeResult::OutOfRange => {},
                        _ => {
                            results.get_mut(&key).unwrap().insert(var.clone(), result);
                        },
                    }
                }
            }
        }
        eprintln!("region {} done", region.name);
    }


    // Sort the read_names
    eprintln!("Sorting read names");
    let variant_scores = generate_variant_scores(&variants);
    let sorter = generate_sort_keys(&all_read_names, &results, &variant_scores);
    all_read_names.sort_by(|a, b| sorter.get(a).cmp(&sorter.get(b)).then_with(|| a.cmp(b)));
    

    // Create an output directory
    std::fs::create_dir_all("output")?;

    
    // Write results
    let outfile_name = "output/variants.tsv";
    eprintln!("Writing to {}", outfile_name);
    let outfile = std::fs::File::create(outfile_name)?;
    let mut writer = std::io::BufWriter::new(outfile);
    let header = format!("CHROM\tPOS\tREF\tALT\tIS_GERMLINE\tIS_INDEL\t{}", &all_read_names.join("\t"));
    writeln!(writer, "{}", header)?;
    for var in &variants {
        let mut s = format!("{}\t{}\t{}\t{}\t{}\t{}\t", var.chrom, var.pos, var.ref_base, var.alt_base, var.is_germline, var.is_indel);
        for read_name in &all_read_names {
            let read_name_lookup = results.get(read_name);
            match read_name_lookup {
                Some(entry) => {
                    if !entry.contains_key(var) {
                        s.push_str(".\t");
                    } else {
                        s.push_str(&format!("{}", entry.get(var).unwrap()));
                        s.push_str("\t");
                    }
                },
                None => {
                    eprintln!("No entry for read {} in results", read_name);
                },
            }
        }

        writeln!(writer, "{}", s.trim_end())?;
    }
    println!("Wrote to {}", outfile_name);

    for k in results.keys() {
        let f = k.replace("/", "!");
        let outfile_name = format!("output/{}.tsv", f);
        println!("Writing to {}", outfile_name);
        let outfile = std::fs::File::create(outfile_name)?;
        let mut writer = std::io::BufWriter::new(outfile);
        
        let header = format!("CHROM\tPOS\tREF\tALT\tIS_GERMLINE\tIS_INDEL\t{}", k);
        writeln!(writer, "{}", header)?;
        
        let mut values = Vec::new();
        for val in results.get(k).unwrap() {
            values.push(val);
        }
        values.sort_by(|a, b| {
            a.0.chrom
                .cmp(&b.0.chrom)
                .then_with(|| a.0.pos.cmp(&b.0.pos))
        });
        for (var, result) in values {
            match result {
                //GenotypeResult::OutOfRange => {},
                _ => {
                    let data = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        var.chrom, var.pos, var.ref_base, var.alt_base, var.is_germline, var.is_indel, result);
                    writeln!(writer, "{}", data)?;
                },
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone)]
enum GenotypeResult {
    Reference(String),
    Alternative(String),
    ThirdAllele(String),
    OutOfRange,
    Ambiguous(String),
}

impl fmt::Display for GenotypeResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            GenotypeResult::Reference(base) => write!(f, "0"),
            GenotypeResult::Alternative(base) => write!(f, "1"),
            GenotypeResult::ThirdAllele(base) => write!(f, "x"),
            GenotypeResult::OutOfRange => write!(f, "."),
            GenotypeResult::Ambiguous(msg) => write!(f, "?"),
        }
    }
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
        ReadPosition::Deletion => Ok(GenotypeResult::Ambiguous("SNP is inside a deletion".to_owned())),
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
                return Ok(GenotypeResult::Ambiguous("Read is ambiguous at first position of deletion".to_owned()));
            }
        },
        _ => return Ok(GenotypeResult::Ambiguous("First base of indel not found in read".to_owned())),
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
            _ => return Ok(GenotypeResult::Ambiguous("Read contains a deletion of different size than the variant".to_owned())),
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
