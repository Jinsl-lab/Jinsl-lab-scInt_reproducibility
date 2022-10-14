# All integration methods to compare, containing MNN, fastMNN, LIGER, online iNMF, Seurat, Harmony, Scanorama, RPCI, Conos, and scMC. 
library(reticulate)
use_condaenv("/home/jsl/anaconda3/envs/r-4.0.2/bin/python3")
Scanorama = import("scanorama")
library(batchelor)
library(conos)
library(harmony)
library(rliger)
library(RISC)
library(Seurat)
library(scMC)
library(irlba)

###################################################################################
### MNN ###
###################################################################################
run_MNN <- function(batches, meta,
                    nfeatures = 2000,
                    npcs = 20,
                    is.normalize = TRUE,
                    vargenes = NULL)
{
  # preprocessing
  b_seurat <- CreateSeuratObject(do.call(cbind, batches), 
                                 min.cells = 0, min.genes = 0)
  cell.names <- colnames(b_seurat)
  rownames(meta) <- cell.names
  b_seurat <- CreateSeuratObject(do.call(cbind, batches), meta.data = meta, 
                                 min.cells = 0, min.genes = 0)
  if (is.normalize) {
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  b_seurat <- ScaleData(object = b_seurat)
  b_seurat <- FindVariableFeatures(object = b_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargenes) == TRUE) {
    VariableFeatures(b_seurat) = vargenes
  }
  # highly variable genes
  hvg_all <- b_seurat@assays[["RNA"]]@var.features
  meta_data <- b_seurat@meta.data
  data_filtered <- b_seurat@assays[["RNA"]]@scale.data[hvg_all,]
  data_list <- list()
  meta_list <- list()
  list_batch <- unique(meta_data[, "batchlb"])
  for (i in 1:length(list_batch)) {
    selected <- row.names(meta_data)[which (meta_data[, "batchlb"] == list_batch[i])]
    data_list[[i]] <- data_filtered[, selected]
    meta_list[[i]] <- meta_data[selected, ]
  }   
  meta <- do.call(rbind, meta_list)
  
  ##########################################################
  # run MNN
  t1 = Sys.time()
  out_mnn_total <- do.call(mnnCorrect, data_list)
  t2 = Sys.time()
  print(t2-t1)
  
  mnn_correct <- out_mnn_total@assays@data@listData[["corrected"]]
  dataAll_pca <- prcomp_irlba(t(mnn_correct), n = npcs)
  mnn_res <- t(dataAll_pca$x)
  return(mnn_res)
}


###################################################################################
### fastMNN ###
###################################################################################
run_fastMNN <- function(batches, meta,
                        nfeatures = 2000,
                        npcs = 20,
                        is.normalize = TRUE,
                        vargenes = NULL)
{
  # preprocessing
  b_seurat <- CreateSeuratObject(do.call(cbind, batches),
                                 min.cells = 0, min.genes = 0)
  cell.names <- colnames(b_seurat)
  rownames(meta) <- cell.names
  b_seurat <- CreateSeuratObject(do.call(cbind, batches), meta.data = meta,
                                 min.cells = 0, min.genes = 0)
  if (is.normalize) {
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  b_seurat <- ScaleData(object = b_seurat)
  b_seurat <- FindVariableFeatures(object = b_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargenes) == TRUE) {
    VariableFeatures(b_seurat) = vargenes
  }
  # highly variable genes
  hvg_all <- b_seurat@assays[["RNA"]]@var.features
  meta_data <- b_seurat@meta.data
  data_filtered <- b_seurat@assays[["RNA"]]@scale.data[hvg_all,]
  data_list <- list()
  meta_list <- list()
  list_batch <- unique(meta_data[, "batchlb"])
  for (i in 1:length(list_batch)) {
    selected <- row.names(meta_data)[which (meta_data[, "batchlb"] == list_batch[i])]
    data_list[[i]] <- data_filtered[, selected]
    meta_list[[i]] <- meta_data[selected, ]
  }   
  meta <- do.call(rbind, meta_list)
  
  ##########################################################
  # run fastMNN
  t1 = Sys.time()
  out_mnn_total <- do.call(batchelor::fastMNN, data_list)
  t2 = Sys.time()
  print(t2-t1)
  
  fastmnn_res <- out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components
  fastmnn_res <- t(fastmnn_res)
  return(fastmnn_res)
}


###################################################################################
### LIGER ###
###################################################################################
run_LIGER <- function(batches, meta,
                      is.normalize = TRUE,
                      nfeatures = 2000,
                      k = 20,
                      lambda = 5,
                      nrep = 1) 
{
  # preprocessing
  liger_object = createLiger(batches, remove.missing = F) 
  if (is.normalize == TRUE) {
    liger_object = liger::normalize(liger_object)
  } else {
    liger_object = liger::normalize(liger_object)
    liger_object@norm.data = liger_object@raw.data
  }
  liger_object = selectGenes(liger_object, var.thresh = 0.1)
  liger_object = scaleNotCenter(liger_object, remove.missing = F)
  
  ##########################################################
  # run LIGER
  t1 = Sys.time()
  liger_object = optimizeALS(liger_object, k = k, lambda = lambda, nrep = nrep)
  liger_object = quantile_norm(liger_object)
  t2 = Sys.time()
  print(t2-t1)
  
  liger_res <- t(liger_object@H.norm)
  return(liger_res)
}


###################################################################################
### online iNMF ###
###################################################################################
run_online_iNMF <- function(batches, meta,
                            is.normalize = TRUE,
                            k = 20,
                            lambda = 5,
                            max.epochs = 1,
                            vargenes = NULL)
{
  # preprocessing
  liger_obj = rliger::createLiger(batches, remove.missing = F)
  if (is.normalize == TRUE) {
    liger_obj = rliger::normalize(liger_obj, remove.missing = F)
  } else {
    liger_obj = rliger::normalize(liger_obj, remove.missing = F)
    liger_obj@norm.data = liger_obj@raw.data
  }
  liger_obj = rliger::selectGenes(liger_obj)
  liger_obj = rliger::scaleNotCenter(liger_obj, remove.missing = F)
  mini = min(floor(min(sapply(batches, ncol))), 5000)
  
  ##########################################################
  # run online iNMF
  t1 = Sys.time()
  liger_obj = rliger::online_iNMF(liger_obj, k = k, max.epochs = max.epochs, lambda = lambda, miniBatch_size = mini)
  liger_obj = rliger::quantile_norm(liger_obj)
  t2 = Sys.time()
  print(t2-t1)
  
  liger_res <- t(liger_obj@H.norm)
  return(liger_res)
}


###################################################################################
### Seurat ###
###################################################################################
run_Seurat <- function(batches, meta,
                       is.normalize = TRUE,
                       nfeatures = 2000,
                       npcs = 20,
                       vargenes = NULL)
{
  # preprocessing
  batches <- do.call(cbind, batches)
  rownames(meta) <- colnames(batches)
  batch_seurat <- CreateSeuratObject(batches, meta.data = meta, 
                                     min.cells = 0, min.features = 0)
  batch_list <- SplitObject(batch_seurat, split.by = "batchlb")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, 
                                            verbose = FALSE)
    print(dim(batch_list[[i]]))
  }
  cell_anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:npcs)
  
  # run Seurat 
  t1 = Sys.time()
  batch_correct <- IntegrateData(anchorset = cell_anchors, dims = 1:npcs)
  t2 = Sys.time()
  print(t2-t1)
  
  DefaultAssay(batch_correct) <- "integrated"
  batch_correct <- ScaleData(object = batch_correct)
  batch_correct <- RunPCA(object = batch_correct, npcs = npcs, verbose = FALSE)
  seurat_res <- t(as.data.frame(batch_correct@reductions$pca@cell.embeddings))
  return(seurat_res)
}


###################################################################################
### Harmony ###
###################################################################################
run_Harmony <- function(batches, meta,
                        is.normalize = TRUE,
                        nfeatures = 2000,
                        npcs = 20,
                        vargenes = NULL)
{
  # preprocessing
  if (is.data.frame(batches[[1]])) {
    batches <- lapply(batches, as.matrix)
  }
  batches <- do.call(cbind, batches)
  rownames(meta) <- colnames(batches)
  batch_seurat <- CreateSeuratObject(batches, meta.data = meta, 
                                     min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(object = batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }            
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", 
                                       nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargenes) == TRUE) {
    VariableFeatures(batch_seurat) = vargenes
  }
  batch_seurat <- ScaleData(object = batch_seurat)
  batch_seurat <- RunPCA(object = batch_seurat, npcs = npcs, features = VariableFeatures(object = batch_seurat))
  
  #run Harmony
  t1 = Sys.time()
  batch_seurat <- RunHarmony(object = batch_seurat, "batchlb", theta = 2, plot_convergence = TRUE, 
                             nclust = 50, max.iter.cluster = 100)
  t2 = Sys.time()
  print(t2-t1)
  
  harmony_res <- t(as.data.frame(batch_seurat@reductions[["harmony"]]@cell.embeddings))
  return(harmony_res)
}

###################################################################################
### Scanorama ###
###################################################################################
run_Scanorama <- function(batches, meta,
                          is.normalize = TRUE,
                          nfeatures = 2000,
                          npcs = 20,
                          vargenes = NULL)
{
  # preprocessing
  batches <- lapply(batches, as.matrix)
  b_seurat <- CreateSeuratObject(do.call(cbind, batches), 
                                 min.cells = 0, min.genes = 0)
  cell.names <- colnames(b_seurat)
  rownames(meta) <- cell.names
  meta$cell <- cell.names
  batch_seurat <- CreateSeuratObject(do.call(cbind, batches), meta.data = meta, 
                                     min.cells = 0, min.genes = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(object = batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }            
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", 
                                       nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargenes) == TRUE) {
    VariableFeatures(batch_seurat) = vargenes
  }
  datasets <- setNames(lapply(SplitObject(batch_seurat, split.by = "batchlb"), function(x) as.matrix(t(x@assays$RNA@data[VariableFeatures(batch_seurat),]))), NULL)
  gene_list <- lapply(unique(batch_seurat$batchlb), function(x) VariableFeatures(batch_seurat))
  
  # run Scanorama
  t1 = Sys.time()
  batch_integrated = Scanorama$integrate(datasets, gene_list)
  t2 = Sys.time()
  print(t2-t1)
  
  scanorama_res <- do.call(rbind, batch_integrated[[1]])
  rownames(scanorama_res) <- do.call(c, lapply(datasets, rownames))
  scanorama_res <- t(x = scanorama_res[colnames(batch_seurat),])
  return(scanorama_res)
}


###################################################################################
### RPCI ###
###################################################################################
run_RPCI <- function(batches, meta,
                     is.normalize = TRUE,
                     nfeatures = 2000,
                     npcs = 20,
                     eigens = 8,
                     method = 'RPCI',
                     align = 'OLS',
                     vargenes = NULL)
{
  # preprocessing
  N <- length(batches)
  metas <- lapply(batches, function(x) meta[colnames(x), ])
  sces <- lapply(X = 1:N, 
                 FUN = function(i) readscdata(count = batches[[i]], cell = metas[[i]], 
                                              gene = data.frame(Symbol = rownames(batches[[i]]), row.names = rownames(batches[[i]])), is.filter = FALSE))
  process0 <- function(obj0){
    obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, is.filter = FALSE)
    if (is.normalize == TRUE) {
      # normalize data
      obj0 = scNormalize(obj0)
    } else {
      # normalize data
      obj0 = scNormalize(obj0)
      obj0@assay[["logcount"]] <- obj0@assay[["count"]]
    }
    # highly variable genes
    obj0 = scDisperse(obj0)
    obj0@vargene = rownames(obj0@assay$count)
    return(obj0)
  }
  sces <- lapply(sces, process0)
  set.seed(1)
  var.genes = rownames(sces[[1]]@assay$count)
  
  # run RPCI
  t1 = Sys.time()
  pcr = scMultiIntegrate(sces, eigens = eigens, var.gene = var.genes, method = 'RPCI', align = 'OLS', npc = npcs, do.fast = TRUE)
  t2 = Sys.time()
  print(t2-t1)
  
  rpci_res <- t(x = pcr@DimReduction[["cell.pls"]])
  return(rpci_res)
}


###################################################################################
### Conos ###
###################################################################################
run_Conos <- function(batches, meta,
                      is.normalize = TRUE,
                      nfeatures = 2000,
                      npcs = 20,
                      vargenes = NULL)
{
  # preprocessing
  b_seurat <- CreateSeuratObject(do.call(cbind, batches), 
                                 min.cells = 0, min.genes = 0)
  cell.names <- colnames(b_seurat)
  rownames(meta) <- cell.names
  meta$cell <- cell.names
  batch_seurat <- CreateSeuratObject(do.call(cbind, batches), meta.data = meta, 
                                     min.cells = 0, min.genes = 0)
  batch_list <- SplitObject(batch_seurat, split.by = "batchlb")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, 
                                            verbose = FALSE)
    batch_list[[i]] <- ScaleData(object = batch_list[[i]])
    batch_list[[i]] <- RunPCA(object = batch_list[[i]], npcs = npcs)
    print(dim(batch_list[[i]]))
  }
  
  # Construct Conos object
  con <- Conos$new(batch_list, n.cores=1)
  
  # run Conos
  t1 = Sys.time()
  # Build joint graph
  con$buildGraph()
  t2 = Sys.time()
  print(t2-t1)
  
  # Find communities
  con$findCommunities()
  # Generate embedding
  conos_res <- t(con$embedGraph(target.dims = npcs, method = "largeVis", verbose = FALSE)[cell.names, ])
  return(conos_res)
}


###################################################################################
### scMC ###
###################################################################################
run_scMC <- function(batches, meta,
                     is.normalize = TRUE,
                     nfeatures = 2000,
                     npcs = 20,
                     vargenes = NULL)
{
  # preprocessing
  batches <- do.call(cbind, batches)
  rownames(meta) <- colnames(batches)
  batch_seurat <- CreateSeuratObject(batches, meta.data = meta,
                                     min.cells = 0, min.features = 0)
  batch_list <- SplitObject(batch_seurat, split.by = "batchlb")
  for(i in 1:length(batch_list)) {
    if (is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, 
                                            verbose = FALSE)
    batch_list[[i]] <- ScaleData(object = batch_list[[i]])
    print(dim(batch_list[[i]]))
  }
  
  # run scMC
  t1 = Sys.time()
  combined <- RunscMC(batch_list) # If an error is reported, set resolution = 0.8
  combined <- RunUMAP(combined, reduction = "scMC", dim = 1:npcs)
  t2 = Sys.time()
  print(t2-t1)
  
  scmc_res <- t(as.data.frame(Embeddings(combined@reductions[["scMC"]])))
  return(scmc_res)
}