#!/usr/bin/env python3

import sys, os, argparse, subprocess, csv
from collections import defaultdict
from textwrap import dedent

def main(args):

	build_report(args)

def build_report(args):
	output = []
	output.append(print_header(args))
	output.append(print_libraries())
	output.append(load_data())
	output.append(sample_statistics())
	output.append(pca())
	output.append(dimensionality_reduction())
	output.append(clustering())
	output.append(cluster_markers())
	output.append(session_info())
	write_output(output)
	run_markdown()

def print_header(args):

	if args.all_resolution_cluster_markers:
		all_resolution_cluster_markers = 'TRUE'
	else:
		all_resolution_cluster_markers = 'FALSE'
	
	return(dedent(
	'''
	---
	title: "Single Cell RNA-Seq Clustering Report"
	author: "Via Scientific"
	output: 
	  html_document:
	    toc: true
	    toc_float:
	      toc_collapsed: true
	      toc_depth: 3
	    fig_caption: yes
	    code_folding: hide
	params:
	  sample_path: "{}"
	  num_pc: {}
	  min_resolution: {}
	  max_resolution: {}
	  filtered_output: ""
	  algorithm: 2
	  find_markers_for_all_resolutions: {}
	---
	''').format(args.sample_path, args.num_pc, args.min_resolution, args.max_resolution, all_resolution_cluster_markers).strip())

def print_libraries():
	
	return(dedent(
	'''
	```{r, setup, message=FALSE, include=FALSE}
	suppressPackageStartupMessages({
	  library(dplyr)
	  library(Seurat)
	  library(ggplot2)
	  library(DoubletFinder)
	  library(data.table)
	  library(DT)
	  library(clustree)
	  library(cluster)
	})
	```'''))

def load_data():
	
	return(dedent(
	'''
	```{r load_data}
	data = readRDS(params$sample_path)

	number_samples = length(unique(data@meta.data$sample))
	number_rows = floor((number_samples+1)/2)
	```'''))

def sample_statistics():
	
	return(dedent(
	'''
	# Sample Statistics {.tabset .tabset-pills}

	Distribution of key properties per cell after filtering:
	
	## Number of Genes
	
	```{r Number of genes per cell by sample, fig.height = 8}
	VlnPlot(data, "nFeature_RNA", group.by= "sample", pt.size = 0) + 
	  ggtitle("Distribution of Genes per Cell") +
	  labs(x="Sample", y="Number of Genes") +
	  NoLegend()
	```
	
	## Number of UMIs
	
	```{r Number of UMIs per cell by sample, fig.height = 8}
	VlnPlot(data, "nCount_RNA", group.by = "sample", pt.size = 0) +
	  ggtitle("Distribution of UMIs per Cell") +
	  labs(x="Sample", y="Number of UMIs") +
	  NoLegend()
	```
	
	## % Mitochondrial
	
	```{r mitochondrial percentage per cell by sample, fig.height = 8}
	VlnPlot(data, "percent.mt", group.by="sample", pt.size = 0) +
	  ggtitle("Distribution of Mitochondrial Read Percentage per Cell") +
	  labs(x="Sample", y="Percent Mitochondrial Reads") +
	  NoLegend()
	```
	
	## % Ribosomal
	
	```{r ribosomal percentage per cell by sample, fig.height = 8}
	VlnPlot(data, "percent.ribo", group.by="sample", pt.size = 0) +
	  ggtitle("Distribution of Ribosomal Read Percentage per Cell") +
	  labs(x="Sample", y="Percent Ribosomal Reads") +
	  NoLegend()
	```'''))

def pca():
	
	return(dedent(
	'''
	# Principal Component Analysis {.tabset .tabset-pills}

	## Dimension Heatmap

	Primary sources of heterogeneity can be observed in this heatmap, which can be useful when trying to decide which principal components should be used in downstream analysis.
	
	```{r Dimension Reduction heatmap, fig.width = 10, fig.height = 10}
	DimHeatmap(data, dims = 1:10, nfeatures = 9, balanced = TRUE, cells = 500, ncol = 2)
	```
	
	## Elbow Plot

	The standard deviation of each principal component is proportional to how much varation the principal component explains.
	
	```{r Elbow Plot}
	reduction = ifelse("harmony" %in% names(data@reductions), yes = "harmony", no = "pca")
	
	elbow_plot_data = data.frame(
	  stdev = data@reductions[[reduction]]@stdev, 
	  PCs = seq(1, length(data@reductions[[reduction]]@stdev))
	)
	pc.changes=(diff((elbow_plot_data$stdev)))*(-1)
	pc.changes.raw=pc.changes
	pc.changes=which(pc.changes >= mean(pc.changes.raw[pc.changes.raw>0]))
	
	ggplot(elbow_plot_data, aes(x=PCs,y=stdev,label=PCs)) +
	  theme_classic()+ 
	  geom_point() + 
	  geom_vline(xintercept = max(pc.changes)+2, color = "firebrick") +
	  geom_vline(xintercept = params$num_pc, color = "steelblue") 
	```'''))

def dimensionality_reduction():
	
	return(dedent(
	'''

	# Dimensionality Reduction {.tabset .tabset-pills}

	Both visualizations are widely used. One consideration is that it is faster to reduce using UMAP than tSNE.

	## UMAP {.tabset .tabset-pills}

	```{r umap, fig.width = 10,fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	dimensions = 1:ifelse(params$num_pc == 0, yes = max(pc.changes)+1, no = params$num_pc)

	data=RunUMAP(data, dims = dimensions, reduction = reduction)
	```
	
	### Colored by Sample {.tabset}
	
	#### Batch Corrected

	```{r UMAP colored by sample, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	if (length(unique(data$sample)) == 1) {
	  DimPlot(data, reduction = "umap") + NoLegend()
	} else {
	  DimPlot(data, reduction = "umap", group.by = "sample")
	}
	```
	
	#### Without Batch Correction
	
	```{r UMAP without batch effect correction, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	temp=data
	temp=RunUMAP(temp, dims = dimensions, reduction = "pca")
	DimPlot(temp,reduction = "umap", group.by = "sample")
	```
	
	### Split by Sample {.tabset}
	
	#### Batch Corrected
	
	```{r UMAP split by sample, fig.width = 10, fig.height = number_rows * 4, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	if (length(unique(data$sample))==1) {
	  DimPlot(data, reduction = "umap") + NoLegend()
	} else {
	  DimPlot(data, reduction = "umap", split.by = "sample", ncol = 2) + NoLegend()
	}
	```
	
	#### Without Batch Correction
	
	```{r UMAP split by sample without batch effect correction, fig.width = 10, fig.height = number_rows * 4, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	DimPlot(temp, reduction = "umap", split.by="sample", ncol=2) + NoLegend()
	rm(temp)
	```

	## tSNE {.tabset .tabset-pills}

	```{r tSNE, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	data=RunTSNE(data, check_duplicates = FALSE, dims = dimensions, reduction = reduction)
	```
	
	### Colored by Sample {.tabset}

	#### Batch Corrected
	
	```{r tSNE colored by sample, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	if (length(unique(data$sample)) == 1){
	  DimPlot(data, reduction = "tsne") + NoLegend()
	} else {
	  DimPlot(data, reduction = "tsne", group.by = "sample")
	}
	```
	
	#### Without Batch Correction
	
	```{r tSNE without batch effect correction, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	temp=data
	temp=RunTSNE(temp, check_duplicates = FALSE, dims = dimensions, reduction = 'pca')
	DimPlot(temp, reduction = "tsne", group.by = "sample")
	```
	
	### Split by Sample {.tabset}
	
	#### Batch Corrected

	```{r tSNE split by sample, fig.width = 10, fig.height = number_rows * 4, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	if (length(unique(data$sample)) == 1) { 
	  DimPlot(data, reduction = "tsne") + NoLegend()
	} else {
	  DimPlot(data, reduction = "tsne", split.by = "sample", ncol=2) + NoLegend()
	}
	```
	
	#### Without Batch Correction
	
	```{r tSNE split by sample without batch effect correction, fig.width = 10, fig.height = number_rows * 4, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	DimPlot(temp, reduction = "tsne", split.by = "sample", ncol=2) + NoLegend()
	rm(temp)	
	```'''))

def clustering():
	
	return(dedent(
	'''
	# Clustering

	## Assessment and Quality Control {.tabset .tabset-pills}
	
	### Processing

	In order to cluster the cells, the shared nearest neighborhood graph of cells are constructed using the top principle components (default is 25).

	And then Graph Based Community Detection Algorithm is used to cluster the cells.
	
	In order to select the best clustering resolution, the sc3 stability index is calculate for each resolution. The resolution with the highest mean sc3 stability index (marked by red line in the figure below).
	
	```{r Build snn graph, error=FALSE, fig.height =10, fig.width=10, message = FALSE, warning = FALSE, results = FALSE}
	
	if ("wpca" %in% names(data@reductions)) {
	  w_reduction = ifelse("wharmony"%in%names(data@reductions),yes="wharmony",no="wpca")
	  
	  elbow_plot_data = data.frame(
	    stdev = data@reductions[[w_reduction]]@stdev,
	    PCs = seq(1, length(data@reductions[[w_reduction]]@stdev))
	  )
	
	  wpc.changes=(diff((elbow_plot_data$stdev)))*(-1)
	  wpc.changes.raw = wpc.changes
	  wpc.changes = which(wpc.changes >= mean(wpc.changes.raw[wpc.changes.raw>0]))
		
	  data = FindMultiModalNeighbors(data, 
	    reduction.list = list(reduction, w_reduction),
	    dims.list = list(1:max(pc.changes+1), 1:max(wpc.changes+1)),
	    modality.weight.name = "multi_modal_weight"
	  )

	  data = RunUMAP(data, nn.name = "weighted.nn")
	  data = RunTSNE(data, check_duplicates = FALSE, nn.name = "weighted.nn")
	
	} else {
	  data=FindNeighbors(data, dims = dimensions, reduction = reduction)
	}
	```
	
	```{r Clustering, fig.width=10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	
	if ("wsnn" %in% names(data@graphs)) {
	  data <- FindClusters(data, graph.name = "wsnn", algorithm = params$algorithm, resolution = seq(params$max_resolution, params$min_resolution, -0.1))
	  names = paste0("wsnn", "_res.")
	} else {
	  data = FindClusters(data, algorithm = params$algorithm, resolution = seq(params$max_resolution, params$min_resolution, -0.1))
	  names = paste0(DefaultAssay(data), "_snn_res.")
	}
	
	SC3_Stability = clustree(data, prefix = names)
	SC3_Stability.results = SC3_Stability$data
	SC3_Stability.results = SC3_Stability.results[, c(names, "sc3_stability")]
	colnames(SC3_Stability.results)[1] = "resolution"
	SC3_Stability.results.mean = aggregate(sc3_stability ~ resolution, SC3_Stability.results, mean)
	colnames(SC3_Stability.results.mean)[2] = "sc3_stability_mean"
	
	if ("wsnn" %in% names(data@graphs)) {
	  Idents(data) = paste0("wsnn","_res.", max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean == max(SC3_Stability.results.mean$sc3_stability_mean)]))
	} else {
	  Idents(data) = paste0(DefaultAssay(data), "_snn_res.", max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean == max(SC3_Stability.results.mean$sc3_stability_mean)]))
	}
	
	stability_selected_resolution = max(as.numeric(as.character(SC3_Stability.results.mean$resolution))[SC3_Stability.results.mean$sc3_stability_mean == max(SC3_Stability.results.mean$sc3_stability_mean)])
	data$seurat_clusters = Idents(data)
	Cluster.distribution = data.frame(table(data$seurat_clusters, data$sample))
	colnames(Cluster.distribution) = c("Cluster", "Sample", "CellNumber")
	```

	### Stability

	```{r Cluster stability assessment, fig.width=10, fig.height = 10, error=FALSE, warning=FALSE, message=FALSE, results=FALSE, echo=FALSE}
	ggplot(SC3_Stability.results, aes(x = resolution,y = sc3_stability)) + 
	  theme_classic() +
	  geom_boxplot(fill='grey90') + 
	  geom_line(data = SC3_Stability.results.mean, aes(x = resolution,y = sc3_stability_mean, group = 1)) + 
	  geom_vline(xintercept = as.factor(stability_selected_resolution), color='firebrick')
	```
	
	### Clustree evaluation
	
	This figure shows how the cells are assigned as the resolution changes. The color of the arrow shows the amount of cells going into the cluster in the next level and the direction of the arrow shows the identity of cluster that the cells are going to.
	
	As the resolution increases, the arrows will start to appear "messy". This means that the clustering algorithm is having trouble assigning cells.
	
	```{r Clustree assessment, fig.width=10, fig.height=10, error=FALSE, warning=FALSE, message=FALSE, results=FALSE, echo=FALSE}
	SC3_Stability
	```

	### Sample Distribution by Cluster

	```{r Sample distribution, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	ggplot(Cluster.distribution, aes(x = CellNumber, y = Cluster, fill = Sample)) +
	  theme_classic() +
	  labs(x = "Fraction of Cells") +
	  geom_bar(stat="identity", position="fill")
	```
	
	### Cluster Distribution by Sample

	```{r Cluster distribution, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	Cluster.distribution = data.frame(table(data$seurat_clusters, data$sample))
	colnames(Cluster.distribution) = c("Cluster", "Sample", "CellNumber")
	ggplot(Cluster.distribution, aes(x=CellNumber, y = Sample, fill = Cluster)) + 
	  theme_classic() +
	  labs(x = "Fraction of Cells") +
	  geom_bar(stat="identity", position="fill")
	```
	
	### Gene Count Distribution by Cluster {.tabset}
	
	#### Violin Plot

	```{r number of genes per cluster v, fig.width = 10, fig.height = 10}
	VlnPlot(data, "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0) +
	  ggtitle("Distribution of Gene Count per Cell by Cluster") +
	  labs(x="Cluster", y="Number of Genes") +
	  NoLegend()	
	```

	#### Ridge Plot

	```{r number of genes per cluster r, fig.width=10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	RidgePlot(data, "nFeature_RNA", group.by = "seurat_clusters") +
	  ggtitle("Distribution of Gene Count per Cell by Cluster") +
	  NoLegend()
	```
	
	### UMI Count Distribution by Cluster {.tabset}
	
	#### Violin Plot

	```{r number of UMIs per cluster v, fig.width=10,fig.height = 10}
	VlnPlot(data,"nCount_RNA", group.by = "seurat_clusters", pt.size = 0) +
	  ggtitle("Distribution of UMI Count per Cell by Cluster") +
	  labs(x="Cluster", y="Number of UMIs") +
	  NoLegend()	
	```
	
	#### Ridge Plot

	```{r number of UMIs per cluster r, fig.width=10,fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	RidgePlot(data, "nCount_RNA", group.by = "seurat_clusters") +
	  ggtitle("Distribution of UMI Count per Cell by Cluster") +
	  NoLegend()
	```

	### % Mitochondrial by Cluster {.tabset}
	
	#### Violin Plot

	```{r mitochondrial percentages per cluster v, fig.width = 10, fig.height = 10}
	VlnPlot(data, "percent.mt", group.by = "seurat_clusters", pt.size = 0) +
	  ggtitle("Distribution of Mitochonrial Reads per Cell by Cluster") +
	  labs(x="Cluster", y="Percent Mitochondrial Reads") +
	  NoLegend()
	```
	
	#### Ridge Plot

	```{r mitochondrial percentages per cluster r, fig.width= 10,fig.height =10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	RidgePlot(data, "percent.mt", group.by = "seurat_clusters") +
	  ggtitle("Distribution of Mitochonrial Reads per Cell by Cluster") +
	  NoLegend()
	```

	### % Ribosomal by Cluster {.tabset}
	
	#### Violin Plot

	```{r ribosomal percentages per cluster v, fig.width = 10, fig.height = 10}
	VlnPlot(data, "percent.ribo", group.by = "seurat_clusters", pt.size = 0) +
	  ggtitle("Distribution of % Ribosomal Reads per Cell by Cluster") +
	  labs(x="Cluster", y="Percent Ribosomal Reads") +
	  NoLegend()
	```

	#### Ridge Plot

	```{r ribosomal percentages per cluster r, fig.width=10,fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = FALSE, echo = FALSE}
	RidgePlot(data, "percent.ribo", group.by = "seurat_clusters") + 
	  ggtitle("Distribution of % Ribosomal Reads per Cell by Cluster") +
	  NoLegend()
	```

	## Visualization {.tabset .tabset-pills}

	### UMAP {.tabset .tabset-pills}

	#### Full Dataset {.tabset}
	
	Astrix indicates the resolution with maximal S3C stability.

	```{r UMAP Visualization of the Cluster, fig.width = 10, fig.height = 10, error = FALSE, warning = FALSE, message = FALSE, results = 'asis', echo = FALSE}
	for (resolution in seq(params$min_resolution, params$max_resolution, .1)) {
	  
	  if (stability_selected_resolution == resolution) {
	    tab_title = paste0('\\n##### ', resolution, '* {.active}\\n')
	  } else {
	    tab_title = paste0('\\n##### ', resolution, '\\n')
	  }
	  
	  cat(tab_title)
	  print(DimPlot(data, reduction = "umap", label = TRUE, group.by=paste0(names, resolution)))
	  cat('\\n')
	}
	```

	```{r UMAP distributed by sample, eval= length(unique(data$sample)) > 1, fig.width = 10, fig.height = number_rows * 4, results='asis', error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}	
	
	cat("\\n#### Split by Sample {.tabset} \\n")
	cat("\\nAstrix indicates the resolution with maximal S3C stability.\\n")

	for (resolution in seq(params$min_resolution, params$max_resolution, .1)) {
	  
	  if (stability_selected_resolution == resolution) {
	    tab_title = paste0('\\n##### ', resolution, '* {.active}\\n')
	  } else {
	    tab_title = paste0('\\n##### ', resolution, '\\n')
	  }
	  cat(tab_title)
	  print(DimPlot(data, reduction = "umap", group.by=paste0(names, resolution), split.by = "sample", ncol = 2))
	  cat('\\n')
	}
	```

	### tSNE {.tabset .tabset-pills}

	#### Full Dataset {.tabset}

	Astrix indicates the resolution with maximal S3C stability.
	
	```{r tSNE Visualization of the Cluster, fig.width = 10, fig.height = 10, results='asis', error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	for (resolution in seq(params$min_resolution, params$max_resolution, .1)) {
	  
	  if (stability_selected_resolution == resolution) {
	    tab_title = paste0('\\n##### ', resolution, '* {.active}\\n')
	  } else {
	    tab_title = paste0('\\n##### ', resolution, '\\n')
	  }
	  
	  cat(tab_title)
	  print(DimPlot(data, reduction = "tsne", label = TRUE, group.by=paste0(names, resolution)))
	  cat('\\n')
	}
	```
	
	```{r tSNE distributed by sample, eval=length(unique(data$sample)) > 1, results='asis', fig.width = 10, fig.height = number_rows * 4, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}

	cat("\\n#### Split by Sample {.tabset}\\n")
	cat("\\nAstrix indicates the resolution with maximal S3C stability.\\n")

	for (resolution in seq(params$min_resolution, params$max_resolution, .1)) {
	  
	  if (stability_selected_resolution == resolution) {
	    tab_title = paste0('\\n##### ', resolution, '* {.active}\\n')
	  } else {
	    tab_title = paste0('\\n##### ', resolution, '\\n')
	  }
	  
	  cat(tab_title)
	  print(DimPlot(data, reduction = "tsne", split.by = "sample", group.by=paste0(names, resolution), ncol = 2))
	  cat('\\n')
	}
	```'''))

def cluster_markers():
	
	return(dedent(
	'''
	## Cluster Markers {.tabset}
	
	Use differential expression analysis to find gene markers for each cluster which can then be used to help identify cell types.
	
	```{r Find Cluster markers, fig.width = 15, fig.height = 20, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	data=NormalizeData(data,assay = "RNA")
	data.markers=FindAllMarkers(data, only.pos = TRUE, assay = "RNA")
	write.table(data.markers, "Cluster.Markers.tsv", quote=FALSE, sep="\t")

	if (params$find_markers_for_all_resolutions) {
	  for (i in colnames(data@meta.data)[grepl("snn_res.", colnames(data@meta.data))]) {
	    temp = data
	    Idents(temp) = i
	    temp = FindAllMarkers(temp, only.pos = TRUE)
	    write.table(temp, paste0(i, ".Cluster.Markers.tsv"), quote=FALSE, sep="\t")
	  }
	}

	if ("cluster" %in% colnames(data.markers)) {
	  top10 = data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
	} else {
	  top10 = NULL
	}

	saveRDS(data, "Final_Analysis.rds")
	```

	### Top gene markers for clusters

	```{r Top gene markers for clusters, fig.width = 15, fig.height = 20, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	DT::datatable(top10 %>% select(Cluster=cluster, Gene=gene, pct.1, pct.2, p_value=p_val, log2FoldChange=avg_log2FC, padj=p_val_adj),
	  rownames = FALSE,
	  extensions = 'Buttons',
	  options=list(
	    dom = 'lftBipr',
	    buttons = list(
	      list(extend = 'csvHtml5', text='Download', filename = "cluster_markers", extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	    )
	  )
	) %>%
	formatRound(columns=c('log2FoldChange'), digits=4) %>%
	formatSignif(columns=c("p_value", "padj"), digits=4)
	```
	
	### Heatmap of top gene markers for clusters
	
	```{r heatmap of Top gene markers for clusters, fig.width = 15, fig.height = 20, error = FALSE, warning = FALSE, message = FALSE, echo = FALSE}
	if (!is.null(top10)) {
	  Vis=ScaleData(data, features = top10$gene, assay = "RNA")
	  DoHeatmap(Vis, features = top10$gene, assay = "RNA") + NoLegend()
	}
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

def write_output(output):

	with open('final_report.Rmd', 'w') as out:
		out.write('\n'.join(output)) 

def run_markdown():

	cmd = "Rscript -e 'rmarkdown::render(\"final_report.Rmd\", \"html_document\")'"
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print(out.decode())
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="build_clustering_and_find_markers.py", description='', usage='%(prog)s [options]')
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-s', '--sample-path', required=True, help='Path to input file.', metavar='', dest='sample_path')
	input_args.add_argument('-m', '--min-resolution', default=0.1, type=float, help='Minimal resolution.', metavar='', dest='min_resolution')
	input_args.add_argument('-n', '--max-resolution', default=2.0, type=float, help='Maximal resolution.', metavar='', dest='max_resolution')
	input_args.add_argument('-p', '--num-pc', default=0, type=int, help='Number of principal components.', metavar='', dest='num_pc')
	input_args.add_argument('-a', '--all-resolution-cluster-markers', action='store_true', help='Find cluster markers for all resolutions.', dest='all_resolution_cluster_markers')

	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)