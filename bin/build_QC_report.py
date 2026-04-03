#!/usr/bin/env python3

import sys, os, argparse, subprocess, csv
from collections import defaultdict
from textwrap import dedent
from multiprocessing import Pool

def main(args):

	build_report(args)

def build_report(args):
	output = []
	output.append(print_header(args))
	output.append(setup())
	output.append(read_input())
	output.append(qc())
	output.append(doublet_removal())
	output.append(filtering())
	output.append(normalization())
	output.append(session_info())
	write_output(output, args.output_prefix)
	run_markdown(args.output_prefix)

def print_header(args):

	if args.remove_mitochondrial_genes:
		remove_mitochondrial_genes = 'TRUE'
	else:
		remove_mitochondrial_genes = 'FALSE'
	
	if args.remove_ribosomal_genes:
		remove_ribosomal_genes = 'TRUE'
	else:
		remove_ribosomal_genes = 'FALSE'

	return(dedent(
	'''
	---
	title: "Quality Control and Filtering Report"
	author: "Via Scientific"
	output: 
	  html_document:
	    toc: true
	    toc_float:
	      toc_collapsed: true
	      toc_depth: 3
	    code_folding: hide
	params:
	  output_prefix: "{}"
	  input_file: "{}"
	  metadata_file: "{}"
	  min_genes: {}
	  max_genes: {}
	  min_UMIs: {}
	  max_UMIs: {}
	  percent_mitochondrial_cutoff: {}
	  percent_ribosomal_cutoff: {}
	  remove_mitochondrial_genes: {}
	  remove_ribosomal_genes: {}
	  doublet_removal: {}
	  doublet_percentage: {}
	  normalization_method: "{}"
	  variable_features: {}
	---
	''').format(args.output_prefix, args.input_file, args.metadata_file, args.min_genes, args.max_genes, args.min_UMIs, args.max_UMIs, args.percent_mitochondrial_cutoff, args.percent_ribosomal_cutoff, args.remove_mitochondrial_genes, args.remove_ribosomal_genes, args.doublet_removal, args.doublet_percentage, args.normalization_method, args.variable_features).strip())

def setup():
	
	return(dedent(
	'''
	```{r setup, include=FALSE}
	suppressPackageStartupMessages({
	  library(dplyr)
	  library(Seurat)
	  library(ggplot2)
	  library(DropletUtils)
	  library(DoubletFinder)
	  library(clustree)
	  library(gridExtra)
	})
	
	get_present_pseudogenes <- function(seurat_object, pseudogenes) {
	  # Ensure that the Seurat object has rownames (gene names)
	  if (is.null(rownames(seurat_object))) {
	    stop("The Seurat object does not have gene names as rownames.")
	  }
	  gene_names <- rownames(seurat_object)
	  pseudogenes_in_dataset <- intersect(pseudogenes, gene_names)
	  return(pseudogenes_in_dataset)
	}
	
	violin_layers = list(
	  theme_classic(base_size = 10),
	  theme(plot.title = element_text(hjust = 0.5),
	        axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        axis.text.x = element_blank(),
	        axis.line.x = element_blank(),
	        )
	)
	```'''))

def read_input():
	
	return(dedent(
	'''
	# Read in samples
	
	If the data is a raw Count Matrix, meaning that the empty droplets are not filtered out by cellranger pipeline or Drop-Seq pipeline, the emptyDrops algorithm from DropUtils will be run in order to remove the empty droplets.
	
	```{r read in samples, message=FALSE}
	data = Read10X_h5(params$input_file)
	
	MultiModal=FALSE
	if (is.list(data)) {
	  MultiModal = TRUE
	  Raw = data
	  data = data[["Gene Expression"]]  
	}
	
	#If any cells with 0 in the dataset, raw input is true
	if (any(colSums(data)==0)) {
		RawInput = TRUE
	} else {
		RawInput = FALSE
	}
	
	if (params$remove_ribosomal_genes) {
	  if (any(grepl("^RP[SL]", rownames(data)))) {
	    rb.genes <- rownames(data)[grep("^RP[SL]", rownames(data))]
	    data = data[!rownames(data) %in% rb.genes,]
	  } else if (any(grepl("^Rp[sl]", rownames(data)))) {
	     rb.genes = rownames(data)[grep("^Rp[sl]", rownames(data))]
	     GTgenes = c("Gm42418","AY036118")
	     rb.genes = c(rb.genes, GTgenes)
	     data = data[!rownames(data) %in% rb.genes,]
	  }
	}
	
	mt.genes = c()
	if (params$remove_mitochondrial_genes) {
	  if (any(grepl("^MT-", rownames(data)))) {
	    mt.genes <- rownames(data)[grep("^MT-",rownames(data))]
	    data = data[!rownames(data) %in% mt.genes,]
	  } else if (any(grepl("^mt-", rownames(data)))) {
	    mt.genes = rownames(data)[grep("^mt-", rownames(data))]
	    data = data[!rownames(data) %in %mt.genes,]
	  }
	}
	
	if (RawInput) {
	  empty = emptyDrops(data[!rownames(data) %in% mt.genes,])
	  empty = data.frame(empty)
	  empty = empty[!is.na(empty$FDR),]
	  empty$DropIdentity = ifelse(empty$FDR<0.001, yes="Non Empty", no="Empty")
	
	  ggplot(empty, aes(x = DropIdentity, y = Total)) +
	    geom_bar(stat = 'identity') +
	    xlab("Droplet classification") +
	    ylab("Number of UMIs per cell") +
	    ggtitle("Empty Droplet classification")
	
	  data=data[, rownames(empty)[empty$FDR<0.05]]
	}
	
	data = CreateSeuratObject(data)
	data$sample = params$output_prefix
	data$orig.ident = params$output_prefix
	
	if (any(grepl("^MT-", rownames(data))) | any(grepl("^mt-", rownames(data)))) {
	    if (any(grepl("^MT-",rownames(data)))) {
	      data[["percent.mt"]] = PercentageFeatureSet(data, pattern="^MT-")
	      rb.genes = rownames(data)[grep("^RP[SL]", rownames(data))]
	      data[["percent.ribo"]] = PercentageFeatureSet(data, features = rb.genes)	
	    } else if (any(grepl("^mt-",rownames(data)))) {
	      data[["percent.mt"]] = PercentageFeatureSet(data,pattern="^mt-")
	      rb.genes = rownames(data)[grep("^Rp[sl]", rownames(data))]
	      present_pseudogenes = get_present_pseudogenes(data, c("Gm42418","AY036118") )
	      rb.genes = c( rb.genes, present_pseudogenes )
	      data[["percent.ribo"]] = PercentageFeatureSet(data, features = rb.genes)	
	    }
	} else {
	  data[["percent.mt"]] = 0
	  data[["percent.ribo"]] = 0
	}
	
	min_gene_threshold = quantile(data$nFeature_RNA, params$min_genes)
	max_gene_cutoff = quantile(data$nFeature_RNA, params$max_genes)
	min_umi_threshold = quantile(data$nCount_RNA, params$min_UMIs)
	max_umi_cutoff = quantile(data$nCount_RNA, params$max_UMIs)
	```'''))

def qc():
	
	return(dedent(
	'''
	# Pre-filtering Quality Control {.tabset .tabset-pills}
	
	Cells with very high mitochondrial and ribosomal contents will bias the downstream clustering and differential expression analysis.
	
	## Distributions
	
	Grey lines represent the selected cutoffs for the respective features.
	
	```{r prefiltering_QC, figure.width=8}
	prefilter_data = data.frame(gene=data$nFeature_RNA, umi=data$nCount_RNA, mitochondria=data$percent.mt, ribosomal=data$percent.ribo)
	
	genes = ggplot(prefilter_data, aes(x=1, y=gene)) +
	  violin_layers +
	  ggtitle("Genes") +
	  labs(y = "Number of Genes") +
	  geom_violin(fill='#fe7f65') +
	  geom_boxplot(width=.1, fill=NA) +
	  geom_hline(yintercept = min_gene_threshold, linetype=2, color='grey') +
	  geom_hline(yintercept = max_gene_cutoff, linetype=2, color='grey')
	
	umi = ggplot(prefilter_data, aes(x=1, y=umi)) +
	  violin_layers +
	  ggtitle("UMIs") +
	  labs(y = "Number of UMIs") +
	  geom_violin(fill='#6697c2') +
	  geom_boxplot(width=.1, fill=NA) +
	  geom_hline(yintercept = min_umi_threshold, linetype=2, color='grey') +
	  geom_hline(yintercept = max_umi_cutoff, linetype=2, color='grey')
	
	mitochondrial = ggplot(prefilter_data, aes(x=1, y=mitochondria)) +
	  violin_layers +
	  ggtitle("% Mitochondrial") +
	  labs(y = "Percent Mitochondiral Reads") +
	  geom_violin(fill='#409537') +
	  geom_boxplot(width=.1, fill=NA) +
	  geom_hline(yintercept = params$percent_mitochondrial_cutoff, linetype=2, color='grey')
	  
	ribosomal = ggplot(prefilter_data, aes(x=1, y=ribosomal)) +
	  violin_layers +
	  ggtitle("% Ribosomal") +
	  labs(y = "Percent Ribosomal Reads") +
	  geom_violin(fill='#fdb23c') +
	  geom_boxplot(width=.1, fill=NA) +
	  geom_hline(yintercept = params$percent_ribosomal_cutoff, linetype=2, color='grey')
	
	grid.arrange(genes, umi, mitochondrial, ribosomal, nrow=1, top='Distributions per Cell')
	```
	
	## Comparisons {.tabset}
	
	### Number of UMIs vs Number of Genes
	
	Grey lines represent the selected cutoffs for the respective features. Cells within the white box will be kept.
	
	``` {r prefiltering_umi_vs_gene}
	ggplot(prefilter_data, aes(x=gene, y=umi)) +
	  theme_classic() +
	  theme(panel.background = element_rect(fill = "grey95")) +
	  geom_rect(aes(xmin = min_gene_threshold, xmax = max_gene_cutoff, ymin = min_umi_threshold, ymax = max_umi_cutoff), fill='white') +
	  labs(x = "Number of Genes", y="Number of UMIs") +
	  geom_point() +
	  geom_vline(xintercept = min_gene_threshold, linetype=2, color='grey') +
	  geom_vline(xintercept = max_gene_cutoff, linetype=2, color='grey') +
	  geom_hline(yintercept = min_umi_threshold, linetype=2, color='grey') +
	  geom_hline(yintercept = max_umi_cutoff, linetype=2, color='grey')
	```
	
	``` {r prefiltering_mitochondria_vs, eval = isFALSE(all(prefilter_data$mitochondria == 0)), results='asis'}  
	cat('\\n### % Mitochondrial vs Other Counts\\n')
	
	cat("Grey lines represent the selected cutoffs for the respective features. Cells within the white box will be kept.")
	
	ggplot(prefilter_data, aes(x = gene, y = mitochondria)) +
	  theme_classic() +
	  theme(panel.background = element_rect(fill = "grey95")) +
	  geom_rect(aes(xmin = min_gene_threshold, xmax = max_gene_cutoff, ymin = 0, ymax = params$percent_mitochondrial_cutoff), fill='white') +
	  labs(x = "Number of Genes", y = "% Mitochondiral Reads") +
	  geom_point() +
	  geom_vline(xintercept = min_gene_threshold, linetype=2, color='grey') +
	  geom_vline(xintercept = max_gene_cutoff, linetype=2, color='grey') +
	  geom_hline(yintercept = params$percent_mitochondrial_cutoff, linetype=2, color='grey')
	
	ggplot(prefilter_data, aes(x = umi, y = mitochondria)) +
	  theme_classic() +
	  theme(panel.background = element_rect(fill = "grey95")) +
	  geom_rect(aes(xmin = min_umi_threshold, xmax = max_umi_cutoff, ymin = 0, ymax = params$percent_mitochondrial_cutoff), fill='white') +
	  labs(x = "Number of UMIs", y= "% Mitochondrial Reads") +
	  geom_point() +
	  geom_vline(xintercept = min_umi_threshold, linetype=2, color='grey') +
	  geom_vline(xintercept = max_umi_cutoff, linetype=2, color='grey') +
	  geom_hline(yintercept = params$percent_mitochondrial_cutoff, linetype=2, color='grey')
	```
	
	``` {r prefiltering_ribosomal_vs, eval = isFALSE(all(prefilter_data$ribosomal == 0)), results='asis'}
	cat('\\n### % Ribosomal vs Other Counts\\n')
	
	cat("Grey lines represent the selected cutoffs for the respective features. Cells within the white box will be kept.")
	
	ggplot(prefilter_data, aes(x = gene, y = ribosomal)) +
	  theme_classic() +
	  theme(panel.background = element_rect(fill = "grey95")) +
	  geom_rect(aes(xmin = min_gene_threshold, xmax = max_gene_cutoff, ymin = 0, ymax = params$percent_ribosomal_cutoff), fill='white') +
	  labs(x = "Number of Genes", y = "% Ribosomal Reads") +
	  geom_point() +
	  geom_vline(xintercept = min_gene_threshold, linetype=2, color='grey') +
	  geom_vline(xintercept = max_gene_cutoff, linetype=2, color='grey') +
	  geom_hline(yintercept = params$percent_ribosomal_cutoff, linetype=2, color='grey')
	
	ggplot(prefilter_data, aes(x = umi, y = ribosomal)) +
	  theme_classic() +
	  theme(panel.background = element_rect(fill = "grey95")) +
	  geom_rect(aes(xmin = min_umi_threshold, xmax = max_umi_cutoff, ymin = 0, ymax = params$percent_ribosomal_cutoff), fill='white') +
	  labs(x = "Number of UMIs", y= "% Ribosomal Reads") +
	  geom_point() +
	  geom_vline(xintercept = min_umi_threshold, linetype=2, color='grey') +
	  geom_vline(xintercept = max_umi_cutoff, linetype=2, color='grey') +
	  geom_hline(yintercept = params$percent_ribosomal_cutoff, linetype=2, color='grey')
	```'''))

def doublet_removal():
	
	return(dedent(
	'''
	# Doublet Classification {.tabset}
	
	Doublets, or sometimes called multiplets, are the droplets which include two or more cells. Including these droplets in the downstream analysis will bias the results because these droplets include gene expression profiles of more than 1 cell. DoubletFinder is used to classify the doublet. Doublet classification can be turned off in the run settings.
	
	```{r Doublet Removal, message=FALSE, results=FALSE, warning=FALSE, fig.show='hide'}
	
	DoubletRemovalHandle=as.logical(params$doublet_removal)
	
	data$Doublet.Classification="SingleLet"
	
	if (is.na(DoubletRemovalHandle)) {
	  if (MultiModal) {
	    if (any(grepl("Multiplexing", names(Raw)))) {
	      DoubletRemovalHandle=FALSE
	    } else {
	      DoubletRemovalHandle=TRUE
	    }
	  } else {
	    DoubletRemovalHandle=TRUE
	  }
	} 
	
	if (DoubletRemovalHandle) {
	  DoubletRemoval = NormalizeData(data)
	  DoubletRemoval = FindVariableFeatures(DoubletRemoval)
	  DoubletRemoval = ScaleData(DoubletRemoval)
	  DoubletRemoval = RunPCA(DoubletRemoval, npcs=100)
	  pc.changes = diff(diff(DoubletRemoval@reductions$pca@stdev))
	  pc.changes = abs(pc.changes)
	  pc.changes = which(pc.changes >= mean(pc.changes))
	  DoubletRemoval = FindNeighbors(DoubletRemoval, dims = 1:(max(pc.changes)+2), reduction = "pca")
	  DoubletRemoval = FindClusters(DoubletRemoval, resolution = seq(2.0, 0.1, -0.1))
	
	  names = paste0(DefaultAssay(DoubletRemoval), "_snn_res.")
	  SC3_Stability = clustree(DoubletRemoval, prefix = names)
	  SC3_Stability.results = SC3_Stability$data
	  SC3_Stability.results = SC3_Stability.results[, c(names, "sc3_stability")]
	  colnames(SC3_Stability.results)[1] = "resolution"
	  SC3_Stability.results.mean = aggregate(sc3_stability ~ resolution, SC3_Stability.results, mean)
	  colnames(SC3_Stability.results.mean)[2] = "sc3_stability_mean"
	  Idents(DoubletRemoval) = paste0(DefaultAssay(DoubletRemoval), "_snn_res.", max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean == max(SC3_Stability.results.mean$sc3_stability_mean)]))
	
	  DoubletRemoval$seurat_clusters = Idents(DoubletRemoval)
	
	  params_selection = paramSweep_v3(DoubletRemoval, PCs = 1:(max(pc.changes)+2), sct = FALSE)
	  params_selection = summarizeSweep(params_selection, GT = FALSE)
	  params_selection = find.pK(params_selection)
	
	  annotations = DoubletRemoval@meta.data$DoubletRemoval$seurat_clusters
	
	  homotypic.prop = modelHomotypic(annotations) 
	  nExp_poi = round(params$doublet_percentage * nrow(DoubletRemoval@meta.data))  
	  nExp_poi.adj = round(nExp_poi * (1 - homotypic.prop))
	
	  DoubletRemoval = doubletFinder_v3(DoubletRemoval, PCs = 1:(max(pc.changes)+2), pN = 0.25, pK = as.numeric(as.character(params_selection$pK))[params_selection$BCmetric == max(params_selection$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	  doublet.classification.name = colnames(DoubletRemoval@meta.data)[ncol(DoubletRemoval@meta.data)]
	  DoubletRemoval$Doublet.Classification = DoubletRemoval@meta.data[, doublet.classification.name]
	
	  data$Doublet.Classification = DoubletRemoval$Doublet.Classification
	}
	```
	
	```{r doublet_removal}
	
	if (DoubletRemovalHandle) {
	  colors = c("Singlet"="#6697c2", "Doublet"="#fe7f65")
	  
	  doublet_data = data.frame(gene=DoubletRemoval$nFeature_RNA, umi=DoubletRemoval$nCount_RNA, doublet_status=data$Doublet.Classification) %>%
	                 mutate(doublet_status = factor(doublet_status, levels=c("Singlet", "Doublet")))
	  
	  num_singlet = nrow(doublet_data %>% filter(doublet_status == 'Singlet'))
	  num_doublet = nrow(doublet_data %>% filter(doublet_status == 'Doublet'))
	  
	  axis_labels <- c("Singlet"=paste0('Singlet (n = ', num_singlet,')'), "Doublet"=paste0('Doublet (n = ', num_doublet,')'))
	  
	  gene_doublet = ggplot(doublet_data, aes(x=doublet_status, y=gene, fill=doublet_status)) +
	    theme_classic() +
	    theme(plot.title = element_text(hjust = 0.5),
	          legend.position = 'none',
	          axis.title.x = element_blank(),
	          axis.ticks.x = element_blank(),
	          axis.line.x = element_blank()) +
	    ggtitle("Gene Count") +
	    labs(y = "Number of Genes") +
	    scale_x_discrete(labels=axis_labels) +
	    scale_fill_manual(values=colors) +
	    geom_violin() +
	    geom_boxplot(width=.1)
	  
	  umi_doublet = ggplot(doublet_data, aes(x=doublet_status, y=umi, fill=doublet_status)) +
	    theme_classic() +
	    theme(plot.title = element_text(hjust = 0.5),
	          legend.position = 'none',
	          axis.title.x = element_blank(),
	          axis.ticks.x = element_blank(),
	          axis.line.x = element_blank()) +
	    ggtitle("UMI Count") +
	    labs(y = "Number of UMIs") +
	    scale_x_discrete(labels=axis_labels) +
	    scale_fill_manual(values=colors) +
	    geom_violin() +
	    geom_boxplot(width=.1)
	  
	  grid.arrange(gene_doublet, umi_doublet, nrow=1, top='Distributions per Cell')
	} else {
	  cat("Doublet removal was not performed.")
	}
	```'''))

def filtering():
	
	return(dedent(
	'''
	# Filtering Summary {.tabset .tabset-pills}
	
	```{r Filtering}
	
	if (DoubletRemovalHandle) {
	  data$Filtering = ifelse(
	      data$Doublet.Classification!='Doublet' &
	      data$nFeature_RNA > min_gene_threshold &
	      data$nFeature_RNA < max_gene_cutoff &
	      data$nCount_RNA > min_umi_threshold &
	      data$nCount_RNA < max_umi_cutoff &
	      data$percent.mt < params$percent_mitochondrial_cutoff &
	      data$percent.ribo < params$percent_ribosomal_cutoff,
	    yes="Keep",
	    no="Drop"
	  )
	} else {
	  data$Filtering = ifelse(
	      data$nFeature_RNA > min_gene_threshold &
	      data$nFeature_RNA < max_gene_cutoff &
	      data$nCount_RNA > min_umi_threshold &
	      data$nCount_RNA < max_umi_cutoff &
	      data$percent.mt < params$percent_mitochondrial_cutoff &
	      data$percent.ribo < params$percent_ribosomal_cutoff,
	    yes="Keep",
	    no="Drop")
	}
	
	filter_status = data.frame(pass_min_gene = ifelse(data$nFeature_RNA > min_gene_threshold, 1, 0)) %>%
	    mutate(pass_max_gene = ifelse(data$nFeature_RNA < max_gene_cutoff, 1, 0)) %>%
	    mutate(pass_min_umi = ifelse(data$nCount_RNA > min_umi_threshold, 1, 0)) %>%
	    mutate(pass_max_umi = ifelse(data$nCount_RNA < max_umi_cutoff, 1, 0)) %>%
	    mutate(pass_mitochondrial = ifelse(data$percent.mt < params$percent_mitochondrial_cutoff, 1, 0)) %>%
	    mutate(pass_ribosomal = ifelse(data$percent.ribo < params$percent_ribosomal_cutoff, 1, 0)) %>%
	    mutate(pass_doublet_detection = ifelse(DoubletRemovalHandle & data$Doublet.Classification == 'Doublet', 0, 1))
	
	filter_summary = filter_status %>%
	    mutate(signature = paste0(pass_min_gene, pass_max_gene, pass_min_umi, pass_max_umi, pass_mitochondrial, pass_ribosomal, pass_doublet_detection)) %>%
	    group_by(signature) %>%
	    summarise(count = n()) %>%
	    mutate(sample = params$output_prefix)
	
	write.table(filter_summary, file=paste0(params$output_prefix, '_filter_summary.tsv'), quote=FALSE, sep='	', row.names = FALSE)
	
	data = subset(data, subset = Filtering == "Keep")
	```
	
	## Overall Cell Fate
	
	```{r Filtering barplot}
	filter_df = filter_summary %>% 
	  mutate(status = ifelse(signature == '1111111', "Pass", "Filtered Out")) %>%
	  group_by(status) %>%
	  summarise(count = sum(count)) %>%
	  mutate(percent = round(count / sum(count)*100, 2)) %>%
	  mutate(label = paste0(count, ' (', percent, '%)')) %>%
	  mutate(status = factor(status, levels=c('Pass', 'Filtered Out')))
	
	ggplot(filter_df, aes(x=status, y=count, fill=status, label=label)) +
	  theme_classic() +
	  theme(axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        legend.position = 'none') +
	  labs(y = "Number of Cells") +
	  scale_fill_manual(values = c("#6697c2", "#fe7f65")) +
	  geom_bar(stat = 'identity') +
	  geom_text(vjust=-.5)
	```
	
	## Distributions after Filtering {.tabset}
	
	```{r postfiltering_QC, figure.width=8}
	postfilter_data = data.frame(gene=data$nFeature_RNA, umi=data$nCount_RNA, mitochondria=data$percent.mt, ribosomal=data$percent.ribo)
	
	genes = ggplot(postfilter_data, aes(x=1, y=gene)) +
	  violin_layers +
	  ggtitle("Genes") +
	  labs(y = "Number of Genes") +
	  geom_violin(fill='#fe7f65') +
	  geom_boxplot(width=.1, fill=NA)
	
	umi = ggplot(postfilter_data, aes(x=1, y=umi)) +
	  violin_layers +
	  ggtitle("UMIs") +
	  labs(y = "Number of UMIs") +
	  geom_violin(fill='#6697c2') +
	  geom_boxplot(width=.1, fill=NA)
	
	mitochondrial = ggplot(postfilter_data, aes(x=1, y=mitochondria)) +
	  violin_layers +
	  ggtitle("% Mitochondrial") +
	  labs(y = "Percent Mitochondiral Reads") +
	  geom_violin(fill='#409537') +
	  geom_boxplot(width=.1, fill=NA)
	  
	ribosomal = ggplot(postfilter_data, aes(x=1, y=ribosomal)) +
	  violin_layers +
	  ggtitle("% Ribosomal") +
	  labs(y = "Percent Ribosomal Reads") +
	  geom_violin(fill='#fdb23c') +
	  geom_boxplot(width=.1, fill=NA)
	
	grid.arrange(genes, umi, mitochondrial, ribosomal, nrow=1, top='Distributions per Cell')
	```
	
	## Individual Criteria Results {.tabset}
	
	### Gene Count Filter
	
	``` {r gene_count_filter_results}
	colors = c("Pass"="#6697c2", "Too Low"="#fe7f65", "Too High"="#fe7f65", "Error"='black')
	
	gene_df = filter_status %>% 
	  mutate(gene_status = case_when(pass_min_gene == 1 & pass_max_gene == 1 ~ "Pass",
	                                 pass_min_gene == 1 & pass_max_gene == 0 ~ "Too High",
	                                 pass_min_gene == 0 & pass_max_gene == 1 ~ "Too Low",
	                                 pass_min_gene == 0 & pass_max_gene == 0 ~ "Error")) %>%
	  group_by(gene_status) %>%
	  summarise(count = n()) %>%
	  mutate(gene_status = factor(gene_status, levels = c('Pass', 'Too High', 'Too Low', 'Error'))) %>%
	  mutate(percent = round(count / sum(count)*100, 2)) %>%
	  mutate(label = paste0(count, ' (', percent, '%)'))
	
	
	ggplot(gene_df, aes(x=gene_status, y=count, fill=gene_status, label=label)) +
	  theme_classic() +
	  theme(axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        legend.position = 'none') +
	  labs(y = "Number of Cells") +
	  scale_fill_manual(values = colors) +
	  geom_bar(stat = 'identity') +
	  geom_text(vjust=-.5)
	```
	
	### UMI Count Filter
	
	``` {r umi_count_filter_results}
	colors = c("Pass"="#6697c2", "Too Low"="#fe7f65", "Too High"="#fe7f65", "Error"='black')
	
	umi_df = filter_status %>% 
	  mutate(umi_status = case_when(pass_min_umi == 1 & pass_max_umi == 1 ~ "Pass",
	                                pass_min_umi == 1 & pass_max_umi == 0 ~ "Too High",
	                                pass_min_umi == 0 & pass_max_umi == 1 ~ "Too Low",
	                                pass_min_umi == 0 & pass_max_umi == 0 ~ "Error")) %>%
	  group_by(umi_status) %>%
	  summarise(count = n()) %>%
	  mutate(umi_status = factor(umi_status, levels = c('Pass', 'Too High', 'Too Low', 'Error'))) %>%
	  mutate(percent = round(count / sum(count)*100, 2)) %>%
	  mutate(label = paste0(count, ' (', percent, '%)'))
	
	
	ggplot(umi_df, aes(x=umi_status, y=count, fill=umi_status, label=label)) +
	  theme_classic() +
	  theme(axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        legend.position = 'none') +
	  labs(y = "Number of Cells") +
	  scale_fill_manual(values = colors) +
	  geom_bar(stat = 'identity') +
	  geom_text(vjust=-.5)
	```
	
	### % Mitochondrial Filter
	
	``` {r mitochondrial_filter_results}
	colors = c("Pass"="#6697c2", "Too Low"="#fe7f65", "Too High"="#fe7f65", "Error"='black')
	
	mitochondrial_df = filter_status %>% 
	  mutate(mitochondrial_status = ifelse(pass_mitochondrial == 1, "Pass", "Too High")) %>%
	  group_by(mitochondrial_status) %>%
	  summarise(count = n()) %>%
	  mutate(mitochondrial_status = factor(mitochondrial_status, levels = c('Pass', 'Too High', 'Too Low', 'Error'))) %>%
	  mutate(percent = round(count / sum(count)*100, 2)) %>%
	  mutate(label = paste0(count, ' (', percent, '%)'))
	
	ggplot(mitochondrial_df, aes(x=mitochondrial_status, y=count, fill=mitochondrial_status, label=label)) +
	  theme_classic() +
	  theme(axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        legend.position = 'none') +
	  labs(y = "Number of Cells") +
	  scale_fill_manual(values = colors) +
	  geom_bar(stat = 'identity') +
	  geom_text(vjust=-.5)
	```
	
	### % Ribosomal Filter
	
	``` {r ribosomal_filter_results}
	colors = c("Pass"="#6697c2", "Too Low"="#fe7f65", "Too High"="#fe7f65", "Error"='black')
	
	ribosomal_df = filter_status %>% 
	  mutate(ribosomal_status = ifelse(pass_ribosomal == 1, "Pass", "Too High")) %>%
	  group_by(ribosomal_status) %>%
	  summarise(count = n()) %>%
	  mutate(ribosomal_status = factor(ribosomal_status, levels = c('Pass', 'Too High', 'Too Low', 'Error'))) %>%
	  mutate(percent = round(count / sum(count)*100, 2)) %>%
	  mutate(label = paste0(count, ' (', percent, '%)'))
	
	ggplot(ribosomal_df, aes(x=ribosomal_status, y=count, fill=ribosomal_status, label=label)) +
	  theme_classic() +
	  theme(axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        legend.position = 'none') +
	  labs(y = "Number of Cells") +
	  scale_fill_manual(values = colors) +
	  geom_bar(stat = 'identity') +
	  geom_text(vjust=-.5)
	```'''))

def normalization():

	return(dedent(
	'''
	# Normalization

	As the violin plots shown in the QC section, the sequencing depth and coverage of each cell in a single cell RNA-Seq dataset vary significantly.
	
	The normalization step normalize the gene expression profile of each cell, which makes them comparable to each other in the downstream analysis.
	
	The SCTransform is recommended as it enhances the biological signature in the data, however it is quite time-consuming and memory-consuming. 
	
	The LogNormalize is very standard practice time-efficient.
	
	```{r Normalization, warning=FALSE, message=FALSE}
	tryCatch(
	  {
	    if (file.exists(params$metadata_file)) {
	      Metadata=read.table(params$metadata_file, sep="\\t", check.names=FALSE, header=TRUE, row.names = NULL)
	      if ("Sample" %in% colnames(Metadata)) {
	        AttributeList = colnames(Metadata)[colnames(Metadata) != "Sample"]
	        if (params$output_prefix %in% Metadata$Sample) {
	          for (i in 1:length(AttributeList)) {
	            data[[AttributeList[i]]] = as.character(Metadata[Metadata$Sample == params$output_prefix, AttributeList[i]])
	          }
	        } else {
	          for (i in 1:length(AttributeList)) {
	            data[[AttributeList[i]]] = ""
	          }
	        }
	      }
	    }
	  },
	  error = function(err) {
	    print("No metadata information is added to dataset")
	  }
	)
	
	if (params$normalization_method == "SCT") {
	  if (all(data[["percent.mt"]] == 0)) {
	    data = SCTransform(data, variable.features.n = params$variable_features)
	  } else {
	    data = SCTransform(data, variable.features.n = params$variable_features, vars.to.regress = "percent.mt")
	  }
	} else {
	  data = NormalizeData(object = data, normalization.method = params$normalization_method)
	}
	
	if (MultiModal) {
	  for (name in names(Raw)[names(Raw) != "Gene Expression"]) {
	    data[[gsub(" ", "", name)]] = CreateAssayObject((Raw[[name]][, colnames(data)]))
	  }
	}
	
	saveRDS(data, paste0(params$output_prefix, '.rds'))
	```'''))

def session_info():

	return(dedent(
	'''
	```{r, session_info, echo=FALSE, results='asis'}
	cat('# Session Info {.tabset .tabset-pills} \\n')
	cat('## Hide\\n')
	cat('## Show\\n')
	sessionInfo()
	```'''))

def write_output(output, prefix):

	with open('%s_filtering_report.Rmd' % prefix, 'w') as out:
		out.write('\n'.join(output)) 

def run_markdown(prefix):

	cmd = "Rscript -e 'rmarkdown::render(\"%s_filtering_report.Rmd\", \"html_document\")'" % (prefix)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print(out.decode())
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="build_QC_report.py", description='', usage='%(prog)s [options]')
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-o', '--output-prefix', required=True, help='Output prefix', metavar='', dest='output_prefix')
	input_args.add_argument('-i', '--input-file', required=True, help='Path to input file.', metavar='', dest='input_file')
	input_args.add_argument('--metadata-file', default="", help='Path to metadata file.', metavar='', dest='metadata_file')
	input_args.add_argument('--min-genes', default=0.01, type=float, help='Cutoff quantile for minimum number of genes in a cell.', metavar='', dest='min_genes')
	input_args.add_argument('--max-genes', default=0.99, type=float, help='Cutoff quantile for maximum number of genes in a cell.', metavar='', dest='max_genes')
	input_args.add_argument('--min-UMIs', default=0.01, type=float, help='Cutoff quantile for minimum number of UMIs in a cell.', metavar='', dest='min_UMIs')
	input_args.add_argument('--max-UMIs', default=0.99, type=float, help='Cutoff quantile for maximum number of UMIs in a cell.', metavar='', dest='max_UMIs')
	input_args.add_argument('--percent-mitochondrial-cutoff', default=50, type=float, help='Mitochondrial read percentage cutoff per cell.', metavar='', dest='percent_mitochondrial_cutoff')
	input_args.add_argument('--percent-ribosomal-cutoff', default=25, type=float, help='Ribosomal read percentage cutoff per cell.', metavar='', dest='percent_ribosomal_cutoff')
	input_args.add_argument('--variable-features', default=3000, type=int, help='Number of variable features to use.', metavar='', dest='variable_features')
	input_args.add_argument('--normalization-method', default='LogNormalize', choices=["LogNormalize","CLR","RC","SCT"], help='Normalization method to use.', metavar='', dest='normalization_method')
	input_args.add_argument('--doublet-removal', default='DEFAULT', choices=["DEFAULT", "TRUE", "FALSE"], help='Runs doublet detection and removal.', metavar='', dest='doublet_removal')
	input_args.add_argument('--doublet-percentage', default=0.01, type=float, help='Doublet percentage.', metavar='', dest='doublet_percentage')
	input_args.add_argument('--remove-mitochondrial-genes', action='store_true', help='Remove mitochondrial genes.', dest='remove_mitochondrial_genes')
	input_args.add_argument('--remove-ribosomal-genes', action='store_true', help='Remove ribosomal genes.', dest='remove_ribosomal_genes')
	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)