#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(openxlsx)
    library(biomaRt)
})

# Parse arguments manually
args_raw <- commandArgs(trailingOnly = TRUE)
args <- list()

# Defaults
args$tissue <- "Immune system"
args$db <- "./ScTypeDB_full.xlsx"

# Helper to parse args
i <- 1
while (i <= length(args_raw)) {
    arg <- args_raw[i]
    if (arg == "--input") {
        args$input <- args_raw[i + 1]
        i <- i + 2
    } else if (arg == "--output") {
        args$output <- args_raw[i + 1]
        i <- i + 2
    } else if (arg == "--organism") {
        args$organism <- args_raw[i + 1]
        if (!args$organism %in% c("human", "mouse", "zebrafish", "d_melanogaster")) {
            stop("Invalid organism. Must be one of: human, mouse, zebrafish, d_melanogaster")
        }
        i <- i + 2
    } else if (arg == "--tissue") {
        args$tissue <- args_raw[i + 1]
        i <- i + 2
    } else if (arg == "--db") {
        args$db <- args_raw[i + 1]
        i <- i + 2
    } else {
        stop(paste("Unknown argument:", arg))
    }
}

# Get script directory
get_script_dir <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(sub("--file=", "", file_arg[1])))
    } else {
        return(".")
    }
}
script_dir <- get_script_dir()

# Validation
if (is.null(args$input)) stop("--input is required")
if (is.null(args$output)) stop("--output is required")
if (is.null(args$organism)) stop("--organism is required")

# Check input file
if (!file.exists(args$input)) {
    stop(paste("Input file not found:", args$input))
}

# Resolve DB path
if (args$db == "./ScTypeDB_full.xlsx" && !file.exists(args$db)) {
    # Check common locations
    possible_paths <- c(
        file.path(script_dir, "ScTypeDB_full.xlsx"),
        "/ScTypeDB_full.xlsx"
    )
    for (p in possible_paths) {
        if (file.exists(p)) {
            args$db <- p
            message(paste("Found DB file at:", p))
            break
        }
    }
}

if (!file.exists(args$db)) {
    stop(paste("DB file not found:", args$db))
}

# Source sctype_score function
sctype_score_script <- file.path(script_dir, "sctype_score_.R")
if (!file.exists(sctype_score_script)) {
     # Try current dir as fallback
     sctype_score_script <- "sctype_score_.R"
}

if (file.exists(sctype_score_script)) {
    source(sctype_score_script)
} else {
    stop(paste("sctype_score_.R not found. looked in:", script_dir, "and ."))
}

# Helper function to parse markers
parse_markers <- function(marker_str) {
    if (is.na(marker_str) || marker_str == "") return(character(0))
    # Replace /// with ,
    marker_str <- gsub("///", ",", marker_str)
    # Remove spaces
    marker_str <- gsub(" ", "", marker_str)
    # Split by comma
    markers <- unlist(strsplit(marker_str, ","))
    # Remove empty
    markers <- markers[markers != ""]
    return(unique(markers))
}

# Main logic
tryCatch({
    
    # Load Seurat Object
    message(paste("Loading Seurat object from", args$input))
    seurat_obj <- readRDS(args$input)
    
    # Read DB
    message(paste("Reading DB from", args$db))
    db <- read.xlsx(args$db)
    
    # Filter by tissue
    message(paste("Filtering for tissue:", args$tissue))
    db_filtered <- db[db$tissueType == args$tissue, ]
    
    if (nrow(db_filtered) == 0) {
        stop(paste("No entries found for tissue:", args$tissue))
    }
    
    # Extract all markers
    all_positive_markers <- unique(unlist(lapply(db_filtered$geneSymbolmore1, parse_markers)))
    all_negative_markers <- unique(unlist(lapply(db_filtered$geneSymbolmore2, parse_markers)))
    all_markers <- unique(c(all_positive_markers, all_negative_markers))
    
    # Check feature presence in Seurat object to decide on conversion
    features <- rownames(seurat_obj)
    
    overlap_human <- sum(all_markers %in% features)
    overlap_target <- 0
    mapped_markers_flat <- c()
    
    use_target <- FALSE
    
    if (args$organism != "human") {
        message(paste("Converting gene symbols from Human to", args$organism))
        
        # ... logic to get gene_map_list ...
        # Copied from previous logic (Step 1-3)
        # Note: I need to reconstruct the logic here carefully to avoid duplication or confusion
        # But I will first perform the mapping to CHECK overlap
        
        message("Step 1: Retrieving Human Ensembl IDs...")
        human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
        human_genes <- getBM(mart = human, filters = "hgnc_symbol", values = all_markers, attributes = c("ensembl_gene_id", "hgnc_symbol"))
        valid_human_ids <- human_genes$ensembl_gene_id
        
        message("Step 2: Retrieving Homologs...")
        species_to <- switch(args$organism, "mouse" = "mouse", "zebrafish" = "zebrafish", "d_melanogaster" = "fly")
        homologs <- getHomologs(ensembl_gene_ids = valid_human_ids, species_from = "human", species_to = species_to)
        target_ensembl_ids <- homologs[, 2]
        
        message("Step 3: Retrieving Target Gene Symbols...")
        target_dataset <- switch(args$organism, "mouse" = "mmusculus_gene_ensembl", "zebrafish" = "drerio_gene_ensembl", "d_melanogaster" = "dmelanogaster_gene_ensembl")
        target_mart <- useEnsembl("ensembl", dataset = target_dataset)
        target_genes <- getBM(mart = target_mart, filters = "ensembl_gene_id", values = target_ensembl_ids, attributes = c("ensembl_gene_id", "external_gene_name"))
        
        colnames(homologs)[1] <- "ensembl_gene_id"; colnames(homologs)[2] <- "target_ensembl_id"
        step1_2 <- merge(human_genes, homologs, by = "ensembl_gene_id")
        colnames(target_genes)[1] <- "target_ensembl_id"; colnames(target_genes)[2] <- "target_symbol"
        final_map <- merge(step1_2, target_genes, by = "target_ensembl_id")
        final_map <- final_map[final_map$target_symbol != "", ]
        
        gene_map_list <- split(final_map$target_symbol, final_map$hgnc_symbol)
        
        # Check overlap
        mapped_markers_flat <- unique(final_map$target_symbol)
        overlap_target <- sum(mapped_markers_flat %in% features)
        
        message(paste("Overlap with Human markers:", overlap_human))
        message(paste("Overlap with Mapped", args$organism, "markers:", overlap_target))
        
        if (overlap_human > overlap_target) {
            message("WARNING: Higher overlap with Human markers detected. The Input data likely uses Human gene symbols (e.g. All Caps) or orthologs mapped to Human symbols.")
            message("Switching to use Human markers directly (ignoring organism conversion).")
            use_target <- FALSE
        } else {
            use_target <- TRUE
        }
        
    } else {
        message("Organism is human. Using DB genes as is.")
        use_target <- FALSE
    }
    
    # Prepare gs_list (gs_positive and gs_negative)
    gs_positive <- list()
    gs_negative <- list()
    
    message("Constructing gene sets...")
    for (i in 1:nrow(db_filtered)) {
        cell_type <- db_filtered$cellName[i]
        pos_markers <- parse_markers(db_filtered$geneSymbolmore1[i])
        neg_markers <- parse_markers(db_filtered$geneSymbolmore2[i])
        
        if (use_target) {
            # Map markers
            map_markers_func <- function(markers, map_list) {
                mapped <- c()
                for (m in markers) {
                    if (m %in% names(map_list)) {
                        mapped <- c(mapped, map_list[[m]])
                    }
                }
                return(unique(mapped))
            }
            pos_markers_final <- map_markers_func(pos_markers, gene_map_list)
            neg_markers_final <- map_markers_func(neg_markers, gene_map_list)
        } else {
            # Human (or Human-like symbols)
            # Ensure proper casing behavior. 
            # If using human markers, we pass them as is (normally All Caps).
            # We will set gene_names_to_uppercase = TRUE later if use_target is FALSE, usually.
            # But wait, Input might be Human symbols but NOT valid.
            # Actually, sc-type default behavior is toupper.
            pos_markers_final <- toupper(pos_markers)
            neg_markers_final <- toupper(neg_markers)
        }
        
        gs_positive[[cell_type]] <- pos_markers_final
        gs_negative[[cell_type]] <- neg_markers_final
    }
    
    # ...
    # Set uppercase flag logic
    # If use_target is TRUE (conversion used), assume Input matches Target (TitleCase usually), so toupper=FALSE.
    # If use_target is FALSE (Human used), assume Input is Human (AllCaps), so toupper=TRUE (or match input?)
    # sc-type default is toupper=TRUE.
    
    sctype_toupper_flag <- !use_target    
    # Run sc-type scoring
    message("Running sc-type scoring...")
    
    assay_name <- "RNA"
    
    # Check if 'scale.data' is present and populated
    has_scale_data <- FALSE
    if (packageVersion("Seurat") >= "5.0.0") {
        if ("scale.data" %in% Layers(seurat_obj, assay = assay_name)) {
             has_scale_data <- TRUE
        }
    } else {
        if (nrow(seurat_obj[[assay_name]]@scale.data) > 0) {
            has_scale_data <- TRUE
        }
    }
    
    if (has_scale_data) {
        message("Using scaled data...")
        if (packageVersion("Seurat") >= "5.0.0") {
            scRNAseqData <- GetAssayData(seurat_obj, assay = assay_name, layer = "scale.data")
        } else {
            scRNAseqData <- GetAssayData(seurat_obj, assay = assay_name, slot = "scale.data")
        }
        scaled_flag <- TRUE
    } else {
        message("Scaled data not found. Using counts...")
        if (packageVersion("Seurat") >= "5.0.0") {
            scRNAseqData <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
        } else {
            scRNAseqData <- GetAssayData(seurat_obj, assay = assay_name, slot = "counts")
        }
        scaled_flag <- FALSE
    }
    
    scRNAseqData <- as.matrix(scRNAseqData)

    es.max <- sctype_score(scRNAseqData = scRNAseqData, 
                           scaled = scaled_flag, 
                           gs = gs_positive, 
                           gs2 = gs_negative,
                           gene_names_to_uppercase = sctype_toupper_flag)
                           
    # Merge by cluster and assign
    message("Assigning cell types to clusters...")
    
    if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
        stop("seurat_clusters column not found in meta.data. Please run clustering (FindClusters) first.")
    }
    
    clusters <- unique(seurat_obj@meta.data$seurat_clusters)
    cL_resutls <- do.call("rbind", lapply(clusters, function(cl){
        cells_in_cluster <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == cl, ])
        es.max.cl <- sort(rowSums(es.max[, cells_in_cluster, drop=FALSE]), decreasing = TRUE)
        
        head(data.frame(cluster = cl, 
                        type = names(es.max.cl), 
                        scores = es.max.cl, 
                        ncells = length(cells_in_cluster)), 10)
    }))
    
    sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    
    # Assign to Seurat object
    seurat_obj@meta.data$sctype_classification <- ""
    for (j in unique(sctype_scores$cluster)) {
        cl_type <- sctype_scores %>% filter(cluster == j)
        seurat_obj@meta.data$sctype_classification[seurat_obj@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
    }
    Idents(seurat_obj) <- "sctype_classification"
    # Save output
    message(paste("Saving results to", args$output))
    saveRDS(seurat_obj, args$output)
    message("Done.")
    
}, error = function(e) {
    message("An error occurred:")
    message(e)
    quit(status = 1)
})
