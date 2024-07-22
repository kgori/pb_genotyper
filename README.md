pb_genotyper
============

## This is part of the nuclear horizontal transfer in CTVT project.

This repository contains code used to genotype PacBio reads from CTVT samples
for a set of specified variant positions.

The code is written in Rust and will need a rust toolchain to compile.

To compile the code, run `cargo build --release` in the root directory of the
repository.

## Example usage
```rust
Usage: pb_genotyper --variants <VARIANTS> --regions <REGIONS> --bamfile <BAMFILE> --output <OUTPUT>

    Options:
      -v, --variants <VARIANTS>
      -r, --regions <REGIONS>
      -b, --bamfile <BAMFILE>
      -o, --output <OUTPUT>
      -h, --help                 Print help
      -V, --version              Print version
```
where `VARIANTS` is a TSV file containing a list of variant positions in the format
`<CHROM> <POS> <REF> <ALT> <is_germline> <is_indel>`, `REGIONS` is a TSV file containing
regions to load from the BAM file in the format `<CHROM> <START> <END> <ID>`, `BAMFILE`
is a BAM file containing the aligned reads to genotype, and `OUTPUT` is the output folder
to write the genotyped reads to.


### Example data
There are some example input files in `data/variants.tsv` and `data/ht_regions.tsv`.
The read data is too large to be included in the repository, but is available
from ENA \[submission in prep.\].
