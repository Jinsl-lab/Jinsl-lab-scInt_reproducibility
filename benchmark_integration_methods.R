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
library(scater)
library(irlba)
library(BiocNeighbors)
library(Matrix)
library(umap)
library(dplyr)
library(Rcpp)

# input data: (1) list of filtered gene by cell expression matrix; (2) data frame of meta information.

# batches
# meta

###################################################################################
### Raw ###
###################################################################################
run_Raw = function(batches, meta, is.normalize = TRUE, out.npcs = 20) {
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  VariableFeatures(batch_seurat) = rownames(batch_seurat)
  batch_seurat = ScaleData(batch_seurat)
  batch_seurat = RunPCA(object = batch_seurat, features = VariableFeatures(batch_seurat), npcs = out.npcs, verbose = FALSE)
  raw_res = t(as.data.frame(batch_seurat@reductions$pca@cell.embeddings))
  colnames(raw_res) = colnames(batch_seurat)
  return(raw_res)
}

###################################################################################
### MNN ###
###################################################################################
# key parameter: k (20 by default)
run_MNN = function(batches, meta, is.normalize = TRUE, nfeatures = 2000,
                    k = 20,
                    out.npcs = 20)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  vargenes <- batch_seurat@assays[["RNA"]]@var.features
  
  # log-normalized data matrices as input
  data_list = lapply(1:length(batches), function(i) batch_seurat@assays[["RNA"]]@data[vargenes, colnames(batches[[i]])])
  ##########################################################
  # run MNN
  t1 = Sys.time()
  out_mnn_total = do.call(mnnCorrect, c(data_list, k = k))
  t2 = Sys.time()
  print(t2-t1)
  
  mnn_correct = out_mnn_total@assays@data@listData[["corrected"]]
  dataAll_pca = prcomp_irlba(t(mnn_correct), n = out.npcs)
  mnn_res = t(dataAll_pca$x)
  colnames(mnn_res) = colnames(batch_seurat)
  return(mnn_res)
}

###################################################################################
### fastMNN ###
###################################################################################
# key parameter: k (20 by default)
run_fastMNN = function(batches, meta, is.normalize = TRUE, nfeatures = 2000,
                    k = 20,
                    out.npcs = 20)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  vargenes <- batch_seurat@assays[["RNA"]]@var.features

  # log-normalized data matrices as input
  data_list = lapply(1:length(batches), function(i) batch_seurat@assays[["RNA"]]@data[vargenes, colnames(batches[[i]])])
  ##########################################################
  # run fastMNN
  t1 = Sys.time()
  out_mnn_total = do.call(batchelor::fastMNN, c(data_list, k = k, d = out.npcs))
  t2 = Sys.time()
  print(t2-t1)
  
  fastmnn_res = t(out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components)
  colnames(fastmnn_res) = colnames(batch_seurat)
  return(fastmnn_res)
}


###################################################################################
### LIGER ###
###################################################################################
# key parameters: k (20 by default), lambda (5 by default)
run_LIGER = function(batches, meta, is.normalize = TRUE,
                      k = 20, lambda = 5, nrep = 1) 
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat = NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_list = SplitObject(batch_seurat, split.by = "batchlb")

  liger_object = createLiger(batches, remove.missing = F)
  liger_object = liger::normalize(liger_object)
  if (is.normalize == TRUE) {
    liger_object@norm.data = lapply(batch_list, function(x) x@assays[["RNA"]]@data)
  } else {
    liger_object@norm.data = liger_object@raw.data
  }
  liger_object = rliger::selectGenes(liger_object, var.thresh = 0.1)
  liger_object = rliger::scaleNotCenter(liger_object, remove.missing = F)
  
  ##########################################################
  # run LIGER
  t1 = Sys.time()
  liger_object = rliger::optimizeALS(liger_object, k = k, lambda = lambda, nrep = nrep)
  liger_object = rliger::quantile_norm(liger_object)
  t2 = Sys.time()
  print(t2-t1)
  
  liger_res = t(liger_object@H.norm)
  colnames(liger_res) = colnames(batch_seurat)
  return(liger_res)
}

###################################################################################
### online iNMF ###
###################################################################################
# key parameters: k (20 by default), lambda (5 by default)
run_online_iNMF = function(batches, meta, is.normalize = TRUE,
                            k = 20, lambda = 5, mini = NULL, max.epochs = 5)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat = NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_list = SplitObject(batch_seurat, split.by = "batchlb")

  liger_object = createLiger(batches, remove.missing = F)
  liger_object = liger::normalize(liger_object)
  if (is.normalize == TRUE) {
    liger_object@norm.data = lapply(batch_list, function(x) x@assays[["RNA"]]@data)
  } else {
    liger_object@norm.data = liger_object@raw.data
  }
  liger_object = rliger::selectGenes(liger_object, var.thresh = 0.1)
  liger_object = rliger::scaleNotCenter(liger_object, remove.missing = F)
  if (is.null(mini)) {
    mini = min(floor(min(sapply(batches, ncol))), 5000) # miniBatch_size be less than batch size
  }
  
  ##########################################################
  # run online iNMF
  t1 = Sys.time()
  liger_object = rliger::online_iNMF(liger_object, k = k, max.epochs = max.epochs, lambda = lambda, miniBatch_size = mini)
  liger_object = rliger::quantile_norm(liger_object)
  t2 = Sys.time()
  print(t2-t1)
  
  onlineiNMF_res = t(liger_object@H.norm)
  colnames(onlineiNMF_res) = colnames(batch_seurat)
  return(onlineiNMF_res)
}


###################################################################################
### Seurat ###
###################################################################################
# k.filter (200 by default)
run_Seurat = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, 
                      k.filter = 200, out.npcs = 20)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  batch_list = SplitObject(batch_seurat, split.by = "batchlb")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  }

  # run Seurat 
  t1 = Sys.time()
  cell_anchors = FindIntegrationAnchors(object.list = batch_list, k.filter = k.filter)
  batch_correct = IntegrateData(anchorset = cell_anchors)
  t2 = Sys.time()
  print(t2-t1)
  
  DefaultAssay(batch_correct) = "integrated"
  batch_correct = ScaleData(object = batch_correct)
  batch_correct = RunPCA(object = batch_correct, npcs = out.npcs, verbose = FALSE)
  seurat_res = t(as.data.frame(batch_correct@reductions$pca@cell.embeddings))
  colnames(seurat_res) = Reduce('c', lapply(batch_list, colnames))
  return(seurat_res)
}

###################################################################################
### Harmony ###
###################################################################################
# key parameters: group.by.vars ("batchlb" by default), theta (2 by default)
run_Harmony = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, group.by.vars = "batchlb", theta = 2,
                        out.npcs = 20)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(object = batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  batch_seurat <- ScaleData(object = batch_seurat)
  batch_seurat = RunPCA(object = batch_seurat, npcs = out.npcs, features = VariableFeatures(object = batch_seurat))
  
  #run Harmony
  t1 = Sys.time()
  batch_seurat = RunHarmony(object = batch_seurat, group.by.vars = group.by.vars, theta = theta, plot_convergence = TRUE, 
                             nclust = 50, max.iter.cluster = 100)
  t2 = Sys.time()
  print(t2-t1)
  
  harmony_res = t(as.data.frame(batch_seurat@reductions[["harmony"]]@cell.embeddings))
  colnames(harmony_res) = colnames(batch_seurat)
  return(harmony_res)
}

###################################################################################
### Scanorama ###
###################################################################################
# key parameter: knn (20 by default)
run_Scanorama = function(batches, meta, is.normalize = TRUE, nfeatures = 2000)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  vargenes <- batch_seurat@assays[["RNA"]]@var.features
  batch_list = SplitObject(batch_seurat, split.by = "batchlb")
  data_list = setNames(lapply(batch_list, function(x) t(as.matrix(x@assays[["RNA"]]@data[vargenes, ]))), NULL)
  gene_list = lapply(1:length(data_list), function(x) vargenes)
  
  # run Scanorama
  t1 = Sys.time()
  batch_integrated = Scanorama$integrate(data_list, gene_list)
  t2 = Sys.time()
  print(t2-t1)
  
  scanorama_res = t(do.call(rbind, batch_integrated[[1]]))
  colnames(scanorama_res) = colnames(batch_seurat)
  return(scanorama_res)
}


###################################################################################
### RPCI ###
###################################################################################
# key parameters: eigens (12 by default), refer (reference batch, the batch containing the most cell types by default)
run_RPCI = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, 
                    eigens = 12, out.npcs = 12, refer = 1)
{
  # preprocessing
  N = length(batches)
  orderr = c(refer, (1:N)[-refer])
  metas <- lapply(batches, function(x) meta[colnames(x), ])
  sces = lapply(X = 1:N, 
                FUN = function(i) readscdata(count = batches[[i]], cell = metas[[i]],
                                             gene = data.frame(Symbol = rownames(batches[[i]]), 
                                                               row.names = rownames(batches[[i]])), 
                                             is.filter = FALSE))
  process0 = function(obj0){
    obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, is.filter = FALSE)
    if (is.normalize == TRUE) {
      obj0 = scNormalize(obj0)
    } else {
      obj0 = scNormalize(obj0)
      obj0@assay[["logcount"]] <- obj0@assay[["count"]]
    }
    obj0 = scDisperse(obj0)
    return(obj0)
  }
  sces = lapply(sces, process0)
  set.seed(1)
  
  # run RPCI
  t1 = Sys.time()
  pcr = scMultiIntegrate(sces[orderr], eigens = eigens, method = 'RPCI', align = 'OLS', npc = out.npcs, do.fast = TRUE)
  t2 = Sys.time()
  print(t2-t1)
  
  rpci_res = t(x = pcr@DimReduction[["cell.pls"]])
  colnames(rpci_res) = Reduce('c', lapply(metas[orderr], rownames))
  rpci_res = rpci_res[, Reduce('c', lapply(metas, rownames))]
  return(rpci_res)
}


###################################################################################
### Conos ###
###################################################################################
# key parameter: k (20 by default)
run_Conos = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, 
                     k = 20,
                     out.npcs = 20)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  batch_list = SplitObject(batch_seurat, split.by = "batchlb")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  }
  batch_list = lapply(batch_list, function(x) ScaleData(x) %>% RunPCA())
  
  # Construct Conos object
  con = Conos$new(batch_list, n.cores=1)
  
  # run Conos
  t1 = Sys.time()
  # Build joint graph
  con$buildGraph(k = k)
  t2 = Sys.time()
  print(t2-t1)
  
  # Find communities
  con$findCommunities()
  # Generate embedding
  conos_res = t(con$embedGraph(target.dims = out.npcs, method = "largeVis", verbose = FALSE)[Reduce('c', lapply(batch_list, colnames)),])
  return(conos_res)
}


###################################################################################
### scMC ###
###################################################################################
# key parameter: TT (0.6 by default), lambda (1 by default)
# If an error is reported, set resolution = 0.8
run_scMC = function(batches, meta, is.normalize = TRUE, nfeatures = 2000,
                    TT = 0.6, lambda = 1, out.npcs = 20, resolution = NULL)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  batch_list = SplitObject(batch_seurat, split.by = "batchlb")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  }
  batch_list = lapply(batch_list, ScaleData)
  
  # run scMC
  t1 = Sys.time()
  combined = RunscMC(batch_list, similarity.cutoff = TT, lambda = lambda, nDims.scMC = out.npcs, resolution = resolution)
  t2 = Sys.time()
  print(t2-t1)
  
  scmc_res = t(as.data.frame(Embeddings(combined@reductions[["scMC"]])))
  colnames(scmc_res) = colnames(batch_seurat)
  return(scmc_res)
}