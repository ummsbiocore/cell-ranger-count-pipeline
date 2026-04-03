Cell Ranger (v8.0.1) is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis. 

Steps:
  1. Cellranger count takes FASTQ files performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis.  cellranger count also processes Feature Barcoding data alongside Gene Expression reads.
  2. Cellranger aggr aggregates outputs from multiple runs of cellranger count, normalizing those runs to the same sequencing depth and then recomputing the feature-barcode matrices and analysis on the combined data. The aggr pipeline can be used to combine data from multiple samples into an experiment-wide feature-barcode matrix and analysis.

Please check following web site for detailed information:  https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

Inputs: 
  - Reads
  - Mate (single or pair)
  - Alternatively, BCL folders and samplesheets (enable run_mkfastq option)

#### Adding Sequences To Reference Genome
1. Choose genome_build you want to use as initial reference sequence.
2. Turn on the “run_Download_Genomic_Sources” option.
3. Activate the “add_sequences_to_reference” option and configure its settings. You can input your custom sequence in FASTA format. If you prefer to provide GTF file, it is also supported. If you don’t provide one, it will be automatically generated.
4. Enable the run_mkref. This step ensures the inclusion of custom FASTA/GTF files into the genome.

Program Versions:
  - Cellranger-v8.0.1
  - fastqc=0.11.8
  - multiqc=1.7
  - bcl2fastq2=2.20.0.422
  - samtools=1.3
  - igvtools=2.5.3
  - bedtools=2.30.0


