#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(dead_code)]

use anyhow::Error;
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::fs::File;
use std::sync::Arc;
use datafusion::prelude::*;
use datafusion::datasource::listing;
use datafusion::datasource::file_format::parquet;
use datafusion::datasource::TableProvider;

mod cigar;
mod regions;

use cigar::{check_cigar_overlap, position_in_read};

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

use tokio;
#[tokio::main]
async fn main() -> Result<(), Error> {
    let ctx = SessionContext::new();
    let session_state = ctx.state();
    let csv_path = "/nfs/users/nfs_k/kg8/lustre/projects/ctvt_horizontal_transfer/2024-01-19_CTVT_Nepal_samples/pacbio/code/pb_explore/data/ht_regions.tsv";
    let csv_options = CsvReadOptions::new().delimiter(b'\t').has_header(true);
    let csv_df = ctx.read_csv(csv_path, csv_options).await?;

    let arrow_path = "/nfs/users/nfs_k/kg8/lustre/projects/ctvt_horizontal_transfer/2024-01-19_CTVT_Nepal_samples/nf_postprocess_cnvs_and_snvs/runs/filtered_20240201/out/results/dataset";
    let arrow_path =  listing::ListingTableUrl::parse(arrow_path)?;
    let file_format = parquet::ParquetFormat::default();
    let listing_options = listing::ListingOptions::new(Arc::new(file_format))
        .with_file_extension(".parquet");

    let resolved_schema = listing_options.infer_schema(&session_state, &arrow_path).await?;
    let config = listing::ListingTableConfig::new(arrow_path)
        .with_listing_options(listing_options)
        .with_schema(resolved_schema);
    let provider = Arc::new(listing::ListingTable::try_new(config)?);

    let join_expr = col("CHROM").eq(col("CHROM"))
        .and(col("POS").gt_eq(col("START")))
        .and(col("POS").lt_eq(col("END")));

// Apply join between provider and CSV DataFrame
    let joined_df = ctx
        .execute_physical_plan(provider)?
        .join(&ctx.collect(csv_df)?, join_expr, JoinType::Inner)?;

    let bam_path = "/lustre/scratch126/casm/team267ms/kg8/projects/pacbio/alignments/2169Tb.1.hifi_reads.aligned.bam";
    let mut bam = IndexedReader::from_path(bam_path).unwrap();
    let mut record = Record::new();

    let chrnames = get_chrom_names(std::path::Path::new(bam_path)).unwrap();
    
    //let f = bam.fetch(("21", 28669568, 28771175));
    let f = bam.fetch(("1", 30237476, 33634687));
    println!("{:?}", f);

    while let Some(result) = bam.read(&mut record) {
        // let record = result.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        if qname == "m84093_240411_112803_s2/102699030/ccs" {
            println!("record = {:?}", record);
            println!("record.qname = {:?}", qname);
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
    let readpos = position_in_read(&record, 30238190).unwrap();
    println!("readpos = {:?}", readpos);
    println!("base = {:?}", record.seq()[readpos] as char);
    println!("qual = {:?}", record.qual()[readpos]);
    Ok(())
}
