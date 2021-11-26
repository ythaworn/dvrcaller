# `DVR calling`

Estimate the pattern of 68 direct variable repeats (DVRs) in _Mycobacterium tuberculosis_ (mtb) genome from short reads data.


## Input

Paired-end reads, with file names of the form `_{1,2}.fastq.gz`


## Usage

There are two steps:
1. Map reads to a reconstructed sequence containing 68 DVRs using `bwa mem`, implemented as a [nextflow](https://www.nextflow.io/) pipeline

`nextflow run read-mapping.nf --reads "<path>/*_{1,2}.fastq.gz" --outdir <dir>`

where `<path>` is the path to the directory containing the reads data, and `<dir>` is the output directory.

2. Call DVRs from depth of mapped reads using an R script

`Rscript dvr-calling.R <in_dir> <out_dir> --plot`

where `<in_dir>` is the directory containing read depth data from step 1, `<out_dir>` is the output directory, and `--plot` plots mapped reads across the DVR regions (useful for checking the inferred DVRs).


## Outputs

`dvr.csv` contains inferred DVR pattern for each sample and its summary, including the number of DVRs deleted, ranges of DVR deletions, quantiles (25%, 50% and 75%) of read depths at DVR deletions, and spoligotype.

Other files, `dvr-full.csv` and `dvr-filt.csv`, contain more information about read depth statistics at all 68 DVR regions.


## Requirements

The nextflow script `read-mapping.nf` requires nextflow, bwa and samtools.

The R script `dvr-calling.R` requires tidyverse, foreach and gridExtra.
