# DVR caller for Mtb

Estimate the pattern of 68 [direct variable repeats (DVRs)](https://doi.org/10.1128/jcm.35.4.907-914.1997) in _Mycobacterium tuberculosis_ (mtb) genome from short-read genome sequencing data.


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

`read-mapping.nf` requires `nextflow`, `bwa` and `samtools`.

`dvr-calling.R` requires `tidyverse`, `foreach` and `gridExtra` packages in `R`.


## Citation

Netikul T, Thawornwattana Y, Mahasirimongkol S, Yanai H, Maung HMW, Chongsuvivatwong V, Palittapongarnpim P. (2022) __Whole-genome single nucleotide variant phylogenetic analysis of _Mycobacterium tuberculosis_ lineage 1 in endemic regions of Asia and Africa.__ _Sci Rep._, 12(1):1565. [[doi:10.1038/s41598-022-05524-0]](https://doi.org/10.1038/s41598-022-05524-0)


