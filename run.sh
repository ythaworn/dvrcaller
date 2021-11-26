#!/bin/bash

# path to working directory
wd=./

# path to directory containing short reads
din=$wd/data/*_{1,2}.fastq.gz
dout_depth=$wd/out/depth
dout=$wd/out/dvr


timer=$SECONDS

nextflow run read-mapping.nf -resume --reads "$din" --outdir $dout_depth

Rscript dvr-calling.R $dout_depth $dout --plot

timer=$(($SECONDS-timer))
printf "Time used: %02d:%02d:%02d\n" "$((timer/3600))" "$((timer/60%60))" "$((timer%60))"
