#!/usr/bin/env nextflow

/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *
 * read-mapping.nf
 * 
 * Nextflow workflow for read mapping using BWA MEM
 * 
 * Requirements:
 * - nextflow
 * - bwa
 * - samtools
 * 
 * Created by Yuttapong Thawornwattana
 * 24/11/2021
 */


/*
 * Defines paths and parameters
 */
basedir = "."
run_id = "test_run"

// input fastq files
params.reads = "$basedir/data/*_{1,2}.fastq.gz"

// reference
params.ref = "$basedir/ref/14723_8_51.6.fasta"

// output dir
params.outdir = "$basedir/out"


/*
 * The reference genome file
 */
ref_genome = file(params.ref)



/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
  .fromFilePairs( params.reads )
  .ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
  .set{ read_pairs }



/*
 * Builds the genome index required by the mapping process
 */
process build_ref_index {
  input:
  file ref from ref_genome
      
  output:
  file "${ref}.*" into ref_index
        
  """
  bwa index -a is $ref
  """
}



/*
 * Maps each read-pair to the reference genome using bwa mem
 * Convert sam to bam
 */
process read_mapping {
  tag "$id"
  publishDir "$params.outdir/depth", mode: 'copy'

  input:
  file ref from ref_genome
  file ref_index
  set id, file(reads) from read_pairs

  output:
  file "${id}_depth.zip"

  """
  bwa mem -c 100 -R "@RG\\tID:${id}\\tSM:${id}\\tPL:Illumina" -M -T 50 $ref $reads > sam_file
  samtools fixmate -O bam sam_file bam_fixmate
  samtools sort -O bam -o ${id}.bam bam_fixmate
  samtools index ${id}.bam
  samtools depth ${id}.bam > ${id}_depth.txt
  zip -q ${id}_depth.zip ${id}_depth.txt
  rm ${id}_depth.txt ${id}.bam sam_file bam_fixmate
  """
}


/*
 * That's all!
 */
workflow.onComplete { 
  println ( workflow.success ? "Done :)" : "Something went wrong :(" )
}
