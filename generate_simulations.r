library("splatter")

simmeta <- function(sim) {
  meta <- data.frame("CellType" = sim$Group,
                    "cell" = sim$Cell,
                    "batch" = as.factor(as.numeric(as.factor(sim$Batch))),
                    "batchlb" = paste0("Batch_", 
                                      as.factor(as.numeric(as.factor(sim$Batch)))))
  rownames(meta) <- meta$"cell"
  return(meta)
}
simdata <- function(sim) {
  meta_sim <- simmeta(sim)
  batch_label <- meta_sim$"batch"
  N <- as.vector(unique(batch_label))
  Nb <- length(N)
  batch_list <- lapply(X = N, 
                       FUN = function(n) assay(sim, "counts")[, which(batch_label == n)])
  names(batch_list) <- as.vector(unique(meta_sim$batchlb))
  return(batch_list)
}

# Simulation 1
Sim1  <- splatSimulate(batchCells  =  c(1000, 1000, 1500), batch.facScale  =  c(0.1, 0.2, 0.3), #c(0.3, 0.3, 0.4), 
                       group.prob  = c(0.03, 0.15, 0.2, 0.2, 0.4, 0.02), method = "groups", verbose = FALSE)
meta_Sim1 <- simmeta(Sim1)
data_Sim1 <- simdata(Sim1)
a51 <- which(meta_Sim1$CellType[which(meta_Sim1$batchlb == "Batch_1")] == "Group6")
a1 <- c(which(meta_Sim1$CellType[which(meta_Sim1$batchlb == "Batch_2")] == "Group5"),
        which(meta_Sim1$CellType[which(meta_Sim1$batchlb == "Batch_2")] == "Group6"))
a6 <- c(which(meta_Sim1$CellType[which(meta_Sim1$batchlb == "Batch_3")] == "Group1"),
        which(meta_Sim1$CellType[which(meta_Sim1$batchlb == "Batch_3")] == "Group5"))
data_Sim1[["Batch_1"]] <- data_Sim1[["Batch_1"]][, -a51]
data_Sim1[["Batch_2"]] <- data_Sim1[["Batch_2"]][, -a1]
data_Sim1[["Batch_3"]] <- data_Sim1[["Batch_3"]][, -a6]
data_Sim1 <- lapply(data_Sim1, function(x) as(x, "sparseMatrix"))
meta_Sim1 <- rbind(meta_Sim1[colnames(data_Sim1[["Batch_1"]]), ], meta_Sim1[colnames(data_Sim1[["Batch_2"]]), ], meta_Sim1[colnames(data_Sim1[["Batch_3"]]), ])
rm(Sim1, a1, a6, a51)

# Simulation 2
Sim2 <- readRDS("/home/jsl/zhouyang/simulations/sim5_qc.Rds")
data_Sim2 <- simdata(Sim2)
meta_Sim2 <- simmeta(Sim2)
rm(Sim2)

# Simulation 3
Sim3 <- readRDS("/home/jsl/zhouyang/simulations/sim6_qc.Rds")
simmeta0 <- function(sim) {
  meta <- data.frame("CellType" = sim$Group,
                     "cell" = sim$Cell,
                     "batch" = as.factor(as.numeric(as.factor(sim$SubBatch))),
                     "batchlb" = as.character(sim$SubBatch))
  rownames(meta) <- meta$"cell"
  
  return(meta)
}
simdata0 <- function(sim) {
  meta_sim <- simmeta0(sim)
  batch_label <- meta_sim$"batchlb"
  N <- as.vector(unique(batch_label))
  Nb <- length(N)
  batch_list <- lapply(X = N, 
                       FUN = function(n) assay(sim, "counts")[, which(batch_label == n)])
  names(batch_list) <- as.vector(unique(meta_sim$batchlb))
  return(batch_list)
}
data_Sim3 <- simdata0(Sim3)
meta_Sim3 <- simmeta0(Sim3)
data_Sim3 <- lapply(data_Sim3, function(x) as(x, "sparseMatrix"))
rm(Sim3)

# Simulation 4
Sim4  <- splatSimulate(batchCells  =  c(1000, 1000, 1000, 1000), batch.facScale  =  c(0.4, 0.5, 0.5, 0.4),
                       group.prob  = c(0.42, 0.42, 0.04, 0.04, 0.04, 0.04), method = "groups", verbose = FALSE)
meta_Sim4 <- simmeta(Sim4)
data_Sim4 <- simdata(Sim4)
rm(Sim4)