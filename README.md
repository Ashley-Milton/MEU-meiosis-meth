# MEU-meiosis-meth

## Overview
Scripts used to process and visualise EM-seq (whole genome DNA methylation) data. [Bismark](https://github.com/FelixKrueger/Bismark) and [methylKit](https://github.com/al2na/methylKit) were the primary tools used.

![MEU_methylation drawio](https://github.com/Ashley-Milton/MEU-meiosis-meth/assets/64052884/c0abda13-81af-489f-9d1e-e5f78be36121)

## Preparing a genome for use with Bismark tools (Bismark Genome Preparation)
The Bismark Genome Preparation script needs to be run to prepare the genome for alignments.
Requirements:
- Directory containing the genome (.fa or .fasta)

See [1_Bismark_Genome_Prep_generic.pbs](/1_Bismark_Genome_Prep_generic.pbs).

## Trimming raw reads (Trimmomatic)
[Trimmomatic](https://github.com/usadellab/Trimmomatic) is used for trimming sequencing reads.
Requirements:
- Directory containing the sequencing read files (.fastq.gz, e.g.)

See [2_Trimmomatic_all_loop_generic.pbs](/2_Trimmomatic_all_loop_generic.pbs).

## Aligning trimmed reads to a genome (Bismark Genome Alignment)
Bismark Genome Alignment aligns the reads to the genome and makes methylation calls.
Requirements:
- Directory containing the genome (.fa or .fasta) and two bisulfite genome subdirectories (generated after running Bismark Genome Preparation)
- Directory containing trimmed read files (.fq.gz, output from Trimmomatic)

See [3_Bismark_Alignment_1_generic.pbs](/3_Bismark_Alignment_1_generic.pbs).

## Merging aligned reads (samtools)
Aligned reads can be merged using samtools merge.
Requirements:
- Directory containing aligned read files (.bam, output from Bismark Genome Alignment)

See [4_samtools_merge_bams_1234_generic.pbs](/4_samtools_merge_bams_1234_generic.pbs).

## Removing duplicate reads (Bismark Deduplicate)
Duplicate reads are removed using Bismark Deduplicate.
Requirements:
- Directory containing merged aligned read files (.bam, output from samtools merge)

See [5_Bismark_Deduplicate_Loop_generic.pbs](/5_Bismark_Deduplicate_Loop_generic.pbs).

## Generating cytosine reports (Bismark Methylation Extractor)
Bismark Methylation Extractor can be used to generate genome-wide cytosine report output. _Every cytosine on both strands will be considered irrespective of whether they were actually covered by any reads in the experiment or not._
Requirements:
- Directory containing the genome (.fa or .fasta) that was used in Bismark Genome Alignment
- Directory containing deduplicated files (.bam, output from Bismark Deduplicate)

See [6_Bismark_Methylation_Extractor_1-4_generic.pbs](/6_Bismark_Methylation_Extractor_1-4_generic.pbs).

## Processing cytosine reports (methylKit)
Once cytosine reports have been generated, they can be processed and visualised in various ways using the R package [methylKit](https://github.com/al2na/methylKit).See [7_methylKit_generic.R](/7_methylKit_generic.R).

## Generating plots
Further custom plots can be generated using R packages such as [ggplot2](https://github.com/tidyverse/ggplot2). See [8_Plots_for_paper_generic.R](/8_Plots_for_paper_generic.R).
