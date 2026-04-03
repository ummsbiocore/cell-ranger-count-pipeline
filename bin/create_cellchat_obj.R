#!/usr/bin/env Rscript
library(CellChat)
library(Seurat)

args<-commandArgs(trailingOnly=TRUE)

obj <- args[[1]]
labels_col <- args[[2]]
db <- args[[3]]
subset_db <- args[[4]]
threads <- as.numeric(args[[5]])
smooth <- args[[6]]
cell_groups <- args[[7]]
cell_groups_sources <- args[[8]]
cell_groups_targets <- args[[9]]
trim <- args[[10]]
min_cells <- args[[11]]

trim <- ifelse(trim == "", NULL, as.numeric(trim))

CellChatDB <- ifelse (db == "human", "CellChatDB.human", "CellChatDB.mouse")
CellChatDB <- get(CellChatDB)

seurat_obj <- readRDS(obj)
obj.list <- SplitObject(seurat_obj, split.by = "Condition")

cellchat.list <- list()
for (i in seq_along(obj.list)) {
    
    # 1. Extract the specific object and its name
    current_obj <- obj.list[[i]]
    condition_name <- names(obj.list)[i]
    
    message(paste("Processing condition:", condition_name))
    
    # 2. Extract data and metadata
    data.input <- GetAssayData(current_obj, assay = 'RNA', layer = 'data')
    
    labels_df <- current_obj[[labels_col]]
    labels <- labels_df[[labels_col]]
    names(labels) <- rownames(labels_df)
    
    meta <- data.frame(labels = labels, row.names = names(labels))
    
    # 3. Create CellChat object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    
    # 4. Set Database
    if (subset_db == "except_nonprotein") {
        CellChatDB.use <- subsetDB(CellChatDB)
    } else if (subset_db == "all") {
        CellChatDB.use <- CellChatDB
    } else {
        CellChatDB.use <- subsetDB(CellChatDB, search = subset_db)
    }
    
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    
    # 5. Overexpressed genes/interactions
    future::plan("multisession", workers = threads)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # 6. Compute Communication Probability
    use_raw <- ifelse(smooth == "true", FALSE, TRUE)
    if (smooth == "true") { cellchat <- smoothData(cellchat, adj = PPI.human) }
    
    cellchat <- computeCommunProb(cellchat, type = "triMean", trim = trim, raw.use = use_raw)
    cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  
    # 7. Pathway Level Inference
    cellchat <- computeCommunProbPathway(cellchat)
    
    # 8. Aggregate Network
    if (cell_groups == "no"){
        cellchat <- aggregateNet(cellchat)
    } else {
        sources_vec <- strsplit(cell_groups_sources, ",", fixed = TRUE)[[1]]
        targets_vec <- strsplit(cell_groups_targets, ",", fixed = TRUE)[[1]]
        cellchat <- aggregateNet(cellchat, sources.use = sources_vec, targets.use = targets_vec)
    }
    
    # 9. Save file using the condition name
    saveRDS(cellchat, file = paste0(condition_name, "_cellchat.rds"))
    
    cellchat.list[[condition_name]] <- cellchat
}

saveRDS(cellchat.list, file = paste0("cellchat_list.rds"))

