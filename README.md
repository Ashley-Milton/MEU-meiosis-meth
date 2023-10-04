# MEU-meiosis-meth

## Overview
Scripts used to process and visualise EM-seq (whole genome DNA methylation) data. [Bismark](https://github.com/FelixKrueger/Bismark) and [methylKit](https://github.com/al2na/methylKit) were the primary tools used.

![MEU_methylation drawio](https://github.com/Ashley-Milton/MEU-meiosis-meth/assets/64052884/c0abda13-81af-489f-9d1e-e5f78be36121)

## Preparing a genome for use with Bismark tools (Bismark Genome Preparation)
The Bismark Genome Preparation script needs to be run to prepare the genome for alignments. Requirements:
- Directory containing the genome (.fa or .fasta)

## Trimming raw reads (Trimmomatic)

## Aligning trimmed reads to a genome (Bismark Genome Alignment)
Bismark Genome Alignment aligns the reads to the genome and makes methylation calls. Requirements:
- Directory containing the genome (.fa or .fasta) and two bisulfite genome subdirectories (generated after running Bismark Genome Preparation)
- Trimmed read files (generated after running Trimmomatic)

## Merging aligned reads (samtools)

## Removing duplicate reads (Bismark Deduplicate)

## Generating CpG reports (Bismark Methylation Extractor)
Bismark Methylation Extractor can be used to generate genome-wide cytosine report output.

## Processing CpG reports (methylKit)

## Generating plots
