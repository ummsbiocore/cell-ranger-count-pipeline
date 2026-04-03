Cell Ranger (v9.0.0) is a suite of pipelines for analyzing data from 10x Genomics Chromium single-cell RNA-seq experiments. It aligns reads, builds feature–barcode matrices, and supports downstream tasks like clustering and gene-expression analysis.

## Table of Contents
- [Features](#features)
- [Suggested pipelines based on library type](#suggested-pipelines-based-on-library-type)
- [Inputs](#inputs)
- [Adding Sequences To Reference Genome](#adding-sequences-to-reference-genome)
- [Example Datasets](#example-datasets)

## Features

Cellranger count is the core pipeline that:
- Takes FASTQ files as input.
- Aligns reads and filters low-quality ones.
- Identifies cellular barcodes and counts UMIs.
- Produces feature–barcode matrices used for clustering and gene-expression analysis.
- Optionally processes Feature Barcoding data alongside standard Gene Expression reads.

Please check following web site for detailed information:  https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count


## Suggested pipelines based on library type

| library type | Cell Ranger Pipeline |
| -------- | -------- |
| 3' Gene Expression     | count     |
| 3' Gene Expression + Antibody/CRISPR Guide Capture     | count     |
| Antibody Capture only     | count     |
| 3' Gene Expression + Cell Multiplexing (+ Antibody/CRISPR Guide Capture)     | multi     |
| Flex Gene Expression (+ Antibody/CRISPR Guide Capture)    | multi     |
| 5' V(D)J only     | vdj     |
| 5' Gene Expression only    | count     |
| Antibody Capture (+ Gene Expression)     | count     |
| CRISPR (FB) + Gene Expression     | count     |
| 5' Gene Expression + V(D)J (+ FB)    | multi     |
| 5' Gene Expression + V(D)J + Antigen Capture (BEAM) (+ Antibody Capture)     | multi     |

## Inputs
  - Reads
  - Mate (single or pair)
  - Alternatively, BCL folders and samplesheets (enable run_mkfastq option)

## Adding Sequences To Reference Genome
1. Choose genome_build you want to use as initial reference sequence.
2. Turn on the “run_Download_Genomic_Sources” option.
3. Activate the “add_sequences_to_reference” option and configure its settings. You can input your custom sequence in FASTA format. If you prefer to provide GTF file, it is also supported. If you don’t provide one, it will be automatically generated.
4. Enable the run_mkref. This step ensures the inclusion of custom FASTA/GTF files into the genome.


## Example Datasets
1. reads: https://www.viafoundry.com/test_data/fastq_10x_pbmc_1k_v3/ 
- merge following files from two lanes: 
	- pbmc_1k_v3_S1_L001_R1_001.fastq.gz,pbmc_1k_v3_S1_L001_R2_001.fastq.gz
	- pbmc_1k_v3_S1_L002_R1_001.fastq.gz,pbmc_1k_v3_S1_L002_R2_001.fastq.gz
2. Mate: pair
3. Genome build: human_hg39_gencode_v32_cellranger_v6
