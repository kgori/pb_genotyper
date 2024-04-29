#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(dead_code)]

use anyhow::Error;
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::fs::File;
use std::path::PathBuf;

mod cigar;
mod regions;

use cigar::{check_cigar_overlap, position_in_read, ReadPosition, insertion_at_position};

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

#[derive(Debug)]
struct Variant {
    chrom: String,
    pos: u32,
    ref_base: char,
    alt_base: char,
    is_germline: bool,
    is_indel: bool,
}

fn open_tsv(file_path: &str) -> Result<csv::Reader<File>, Error> {
    let file = File::open(file_path)?;
    let reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(file);
    Ok(reader)
}

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
            ref_base: ref_base.chars().next().ok_or_else(|| Error::msg("ref_base is empty"))?,
            alt_base: alt_base.chars().next().ok_or_else(|| Error::msg("alt_base is empty"))?,
            is_germline: is_germline.to_lowercase().parse()?,
            is_indel: is_indel.to_lowercase().parse()?,
        })
    }).collect()
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start: u32,
    end: u32,
    name: String,
}

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


fn main() -> Result<(), Error> {
    let variants_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/code/pb_explore/data/variants.tsv";
    let variants = load_variants(variants_path)?;
    println!("{:?}", variants[0]);
    
    let regions_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/code/pb_explore/data/ht_regions.tsv";
    let mut regions = load_regions(regions_path)?;
    regions.sort_by(|a, b| a.name.cmp(&b.name));
    println!("{:?}", regions[0]);

    let my_region = get_region_with_name(&regions, "C").unwrap();

    let bam_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/alignments/2169Tb.1.hifi_reads.aligned.bam";
    let mut bam = IndexedReader::from_path(bam_path).unwrap();
    let mut record = Record::new();

    let chrnames = get_chrom_names(std::path::Path::new(bam_path)).unwrap();
    
    // let f = bam.fetch(("21", 28669568, 28771175));
    // let f = bam.fetch(("1", 30237476, 33634687));
    let f = fetch_bam_reads_from_region(&mut bam, &my_region)?;
    println!("{:?}", f);

    while let Some(result) = bam.read(&mut record) {
        // let record = result.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        if qname == "m84093_240411_112803_s2/111285790/ccs" {
            println!("record = {:?}", record);
            println!("record.qname = {:?}", qname);
            println!("{}", String::from_utf8(record.seq().as_bytes().to_vec()).unwrap());
            break;
        }
    }
    println!("record = {:?}", record);
    println!("record.qname = {:?}", String::from_utf8(record.qname().to_vec()));
    println!("is supplementary? = {:?}", record.is_supplementary());
    println!("cigar = {:?}", record.cigar());
    println!("cigar string = {:?}", record.cigar().to_string());
    println!("pos = {:?}", record.pos());

    let overlap = check_cigar_overlap(&record, 30237476, 33634687);
    println!("overlap segment C = {:?}", overlap);
    let overlap = check_cigar_overlap(&record, 28669568, 28771175);
    println!("overlap segment D = {:?}", overlap);

    for refpos in 30243800..30243880 {
        let ins_found = insertion_at_position(&record, refpos);
        match ins_found {
            Some(ins) => {
                let seq = String::from_utf8(ins.seq.to_vec()).unwrap();
                println!("ins_found = {:?}, {:?}", ins, seq);
            },
            None => {}
        }
    }
    Ok(())
}
