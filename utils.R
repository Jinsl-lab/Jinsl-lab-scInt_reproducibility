############################ plot
library(ggplot2)
library(cowplot)
library(patchwork)
require(RcppAnnoy)
library(RColorBrewer)
require(ggplot2)
require(grid)
require(ggpubr)
######################################################## Visualization using R ggplot2 package
visual_plot <- function(res_method, meta, method.use = "", dataset.use = "", output.dir = NULL, reduction = "umap", color_used = NULL, 
                        color_Cell_type = NULL, color_Batch = NULL, celltype_title = NULL) {
  library(ggplot2)
  library(cowplot)
  if (is.null(color_used) == T) {
    color_used = c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076", "#999999", "#8B008B", "#63BAAB", "#FF7F00",
                   "#E19F20", "#BDB76B", "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", "#DCDCDC",
                   "#A65628", "#8B7D6B", "#86115A", "#474747", "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", "#377EB8", "#B3B3B3", "#E5C494", "#FFD92F",
                   "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  }
  if (is.null(color_Cell_type) == T) {
    color_Cell_type <- color_used
  }
  if (is.null(color_Batch) == T) {
    color_Batch <- rev(color_used)
  }
  if (reduction == "umap") {
    library(umap)
    out_umap <- uwot::umap(t(res_method), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
                           metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1,
                           min_dist = .1, n_threads = 4, ret_model = TRUE)
    umap_df<- as.data.frame(out_umap$embedding)
    umap_df$'batchlb' <- as.factor(meta$batchlb)
    if (is.numeric(meta$CellType) == TRUE) {
      umap_df$'CellType' <- meta$CellType
    } else {
      umap_df$'CellType' <- as.factor(meta$CellType)
    }
    rownames(umap_df) <- meta$cell
    colnames(umap_df) <- c('UMAP_1', 'UMAP_2', 'batchlb', 'CellType')
    
    if (is.numeric(umap_df$'CellType') == TRUE) {
      p01 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = batchlb)) + geom_point(size = 1, alpha = 0.6) + theme_bw() 
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      library(RColorBrewer)
      p02 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 1, aes(colour = CellType), alpha = 0.6) + theme_bw() +
        scale_colour_gradientn(colours = topo.colors(10))
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)###########################celltype = NULL
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    } else {
      p01 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = batchlb)) + geom_point(size = 1, alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Batch, name = 'Batch')
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      
      p02 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = CellType)) + geom_point(size = 1, alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Cell_type, name = 'Cell type')
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    }
  } else if (reduction == "tsne") {
    library(Rtsne)
    out_tsne <- Rtsne(t(res_method))
    tsne_df<- as.data.frame(out_tsne$Y)
    tsne_df$'batchlb' <- as.factor(meta$batchlb)
    if (is.numeric(meta$CellType) == TRUE) {
      tsne_df$'CellType' <- meta$CellType
    } else {
      tsne_df$'CellType' <- as.factor(meta$CellType)
    }
    rownames(tsne_df) <- meta$cell
    colnames(tsne_df) <- c('TSNE_1', 'TSNE_2', 'batchlb', 'CellType')
    
    if (is.numeric(tsne_df$'CellType') == TRUE) {
      p01 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) +
        scale_color_manual(values = color_Batch, name = 'Batch')
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      library(RColorBrewer)
      p02 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2)) + geom_point(aes(colour = CellType), alpha = 0.6) + theme_bw()  +
        scale_color_manual(values = color_Cell_type, name = 'Cell type')
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    } else {
      p01 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Batch, name = 'Batch')
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      
      p02 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Cell_type, name = 'Cell type')
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    }
  }
  return(res_p)
}

dim.plot <- function (object, meta, colFactor = NULL, col.rev = F, title = NULL, ord = NULL,
                      genes = NULL, legend = TRUE, Colors = NULL, size = 0.5, Alpha = 0.8, 
                      plot.ncol = NULL, raw.count = FALSE, exp.range = NULL, exp.col = "firebrick2", 
                      dim1 = "UMAP1", dim2 = "UMAP2", dim.name = F,
                      label = FALSE, adjust.label = 0.25, label.font = 5) 
{
  X = Y = label0 = as.character()
  exp.col = as.character(exp.col)
  adjust0 = as.numeric(adjust.label)
  font0 = as.integer(label.font)
  dimReduce0 = object
  size0 = as.numeric(size)
  
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  colFactor = as.character(colFactor)
  colData = meta
  m0$factor = as.factor(colData[, colFactor])
  if (!is.null(ord)) {
    m0$factor <- factor(m0$factor, levels=ord)
  }
  Cols = c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076", "#999999", "#8B008B", "#63BAAB", "#FF7F00",
           "#E19F20", "#BDB76B", "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", "#DCDCDC",
           "#A65628", "#8B7D6B", "#86115A", "#474747", "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", "#377EB8", "#B3B3B3", "#E5C494", "#FFD92F",
           "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  if (col.rev) {
    Cols = rev(Cols)
  }
  Cols0 = Cols[1:length(levels(m0$factor))]
  if (!isTRUE(dim.name)) {
    dim1 = NULL; dim2 = NULL
  }
  
  if (length(levels(m0$factor)) <= 35 | !is.null(Colors)) {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[,2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[,2]
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)), alpha = Alpha, size = size0) + 
        geom_text(data = m1, aes(X, Y, label = label0), nudge_x = adjust0, nudge_y = adjust0, size = font0) + 
        scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)), alpha = Alpha, size = size0) + 
        scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color="black", size=16, hjust = 0.5),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
  }
  else {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[,2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[,2]
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)), 
                   alpha = Alpha, size = size0) + 
        geom_text(data = m1, aes(X, Y, label = label0), 
                  nudge_x = adjust0, nudge_y = adjust0, size = font0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5)) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)),alpha = Alpha, size = size0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5)) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
  }
  if (!isTRUE(legend)) plot = plot + theme(legend.position="none")
  return(plot)
}

feature.plot <- function (object, expr.mat, genes = NULL, legend = TRUE, 
                          Colors = NULL, size = 0.5, Alpha = 0.8, plot.ncol = NULL, 
                          raw.count = FALSE, exp.range = NULL, exp.col = "firebrick2") 
{
  X = Y = as.character()
  exp.col = as.character(exp.col)
  size0 = as.numeric(size)
  
  dimReduce0 = object
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  if (is.null(genes)) {
    stop("Please input gene symbol")
  }
  
  count0 = expr.mat
  genes = intersect(genes, rownames(count0))
  if (length(genes) == 1) {
    count0 = as.matrix(count0[genes, ])
    if (is.null(exp.range)) {
      m0$draw = as.numeric(count0)
    }
    else {
      m0$draw = as.numeric(as.matrix(count0))
      min0 = exp.range[1]
      max0 = exp.range[2]
      m0$draw[m0$draw < min0] = min0
      m0$draw[m0$draw > max0] = max0
    }
    if (isTRUE(legend)) {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = draw), alpha = Alpha, size = size0) + 
        scale_color_gradient(low = "gray90", high = exp.col) +
        labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = genes) +
        theme_bw(base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = draw), alpha = Alpha, size = size0) + 
        scale_color_gradient(low = "gray90", high = exp.col, guide = "none") + 
        labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = genes) + 
        theme_bw(base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5))
    }
    return(plot)
  }
  else {
    count0 = as.matrix(count0[genes, ])
    if (is.null(exp.range)) {
      count0 = count0
    }
    else {
      min0 = exp.range[1]
      max0 = exp.range[2]
      count0[count0 < min0] = min0
      count0[count0 > max0] = max0
    }
    if (isTRUE(legend)) {
      plot0 = lapply(genes, function(z) {
        ggplot(m0, aes(X, Y)) + 
          geom_point(aes(color = count0[z,]), alpha = Alpha, size = size0) + 
          scale_color_gradient(low = "gray90", high = exp.col) + 
          labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = z) +
          theme_bw(base_line_size = 0) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(color="black", size=16, hjust = 0.5))
      })
    }
    else {
      plot0 = lapply(genes, function(z) {
        ggplot(m0, aes(X, Y)) + 
          geom_point(aes(color = count0[z, ]), alpha = Alpha, size = size0) + 
          scale_color_gradient(low = "gray90", high = exp.col, guide = "none") +
          labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = z) +
          theme_bw(base_line_size = 0) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(color="black", size=16, hjust = 0.5))
      })
    }
    names(plot0) = genes
    if (is.null(plot.ncol)) {
      return(do.call(grid.arrange, c(plot0, ncol = ceiling(sqrt(length(genes))))))
    }
    else {
      return(do.call(grid.arrange, c(plot0, ncol = plot.ncol)))
    }
  }
}

Color0 <- function(n, col.rev = F) {
  cc <- c("#E41A1C", "#FF664F", "#00525B", "#377EB8", "#004B80", "#875300", "#E19F20", "#4DAF4A", "#005600", "#2A628F", 
          "#984EA3", "#DE8EE8", "#00C9BF", "#005B58", "#FF7F00", "#006B5E", "#63BAAB", "#FFFF33", "#586500", "#006B5F",
          "#A65628", "#EF9563", "#00C9AA", "#005B46", "#F781BF", "#86115A", "#007664", "#00C9B2", "#999999", "#474747", "#00754B", "#00C896")
  color_used <- c(
    "cyan4",      "skyblue3",   "darkolivegreen3",   "lightpink",   "darkmagenta",   "brown",   "blueviolet", "bisque4",  "deeppink3",       "darkkhaki",      
    "dodgerblue4",     "goldenrod4",            "gainsboro",       "firebrick4",      "cadetblue3",
    "greenyellow",     "gray6",           "coral2",                     "yellow4",         
    "darkgoldenrod3",  "navy",            "deepskyblue3","antiquewhite3"
  )
  col_f <- function(x) {
    cl <- rgb2hsv(col2rgb(x))
    hsv(h = cl[1,1], s = cl[2,1], v = cl[3,1])
  }
  cc <- c(cc, sapply(color_used, col_f))
  #scales::show_col(cc[order(cc, decreasing = T)])
  cc <- cc[order(cc, decreasing = T)][1:36]
  names(cc) <- NULL
  colors <- c()
  for(i in 1:6) {
    for(j in 1:6) {
      colors <- c(colors, cc[6*(j-1) + i])
    }
  }
  colors <- c(colors, rev(brewer.pal(8, "Set2")))
  Col = col2rgb(colors)
  dist.color = as.matrix(dist(t(Col)))
  diag(dist.color) = 1e10
  while(length(colors) > n) {
    minCol = apply(dist.color, 1, FUN = min)
    ids = which(minCol == min(minCol))[1]
    dist.color = dist.color[-ids, -ids]
    colors = colors[-ids]
  }
  if(col.rev) {
    colors <- rev(colors)
  }
  return(colors)
}

FeaturePlots<- function(obj, features, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[['RNA']]@data[feature, ])
  maximal<- max(obj[['RNA']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  return(ps)
}
FeaturePlots<- function(obj, features, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  # the minimal and maximal of the value to make the legend scale the same. 
  ps<- list()
  for (feature in features) {
    minimal<- max(-3, min(obj[['RNA']]@data[feature, ]))
    maximal<- min(3, max(obj[['RNA']]@data[feature, ]))
    subset_indx<- obj@meta.data[, metadata_column] == feature
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(feature) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[feature]]<- p
  }
  return(ps)
}

palette_colors <- function(alpha) {
  color = c('#ff6666', '#a30000', '#ff9d00', '#c3ff00', '#3a470e', '#34a321', 
            '#2198a3', '#001547', '#3400a3', '#9666ff', '#470e28', '#a3416d')
  color_rgb <- t(col2rgb(color)/255)
  grDevices::rgb(color_rgb[, 1], color_rgb[, 2], color_rgb[, 3], alpha = alpha)
}

# test functions
ff <- function(i, j, meta) {
  d1 <- as.matrix(S[[i]][[j]])
  d2 <- d1
  a1 <- which(d2>0, arr.ind = T)
  a2 <- a1
  a2 <- as.data.frame(a2)
  print(nrow(a2))
  a3 <- a2
  batches <- unique(meta$batchlb)
  a3[,1] <- as.character(meta$CellType[which(meta$batchlb == batches[i])][a3[,1]])
  a3[,2] <- as.character(meta$CellType[which(meta$batchlb == batches[j])][a3[,2]])
  print(length(which(a3[,1] == a3[, 2])))
  print(length(which(a3[,1] == a3[, 2]))/nrow(a2))
  print(summary(as.factor(a3[,1])))
  print(sum(summary(as.factor(a3[,1]))))
  print(summary(as.factor(a3[,2])))
  print(sum(summary(as.factor(a3[,2]))))
}

fff <- function(S, meta, r = 0) {
  if (!"dgCMatrix" %in% class(S)){
    S <- as(S, "dgCMatrix")
  }
  uS <- upper_tri(S@x, S@p, S@i, ncol(S))
  S <- sparseMatrix(i = as.vector(uS$i), j = as.vector(uS$j), x = as.vector(uS$x),
                    dims = c(ncol(S), ncol(S)), repr = "C")
  if (!"dgCMatrix" %in% class(S)){
    S <- as(S, "dgCMatrix")
  }
  if (r > 0) {
    S@x <- pmax(0, S@x - r)
  }
  a1 <- sparse_positive(S@x, S@p, S@i, ncol(S), length(which(S@x>0)))
  a1 <- as.data.frame(a1)
  print(nrow(a1))
  a1[,1] <- as.character(meta$CellType[a1[,1]])
  a1[,2] <- as.character(meta$CellType[a1[,2]])
  print(length(which(a1[,1] == a1[, 2])))
  print(length(which(a1[,1] == a1[, 2]))/nrow(a1))
  print(summary(as.factor(a1[,1])))
  print(sum(summary(as.factor(a1[,1]))))
  print(summary(as.factor(a1[,2])))
  print(sum(summary(as.factor(a1[,2]))))
}
ff0 <- function(S, meta) {
  
  for (i in 1:(length(S)-1)) {
    for (j in (i+1):length(S)) {
      d1 <- as.matrix(S[[i]][[j]])
      d2 <- d1
      a1 <- which(d2>0, arr.ind = T)
      
    }
  }
  
  a2 <- a1
  a2 <- as.data.frame(a2)
  print(nrow(a2))
  a3 <- a2
  batches <- unique(meta$batchlb)
  a3[,1] <- as.character(meta$CellType[which(meta$batchlb == batches[i])][a3[,1]])
  a3[,2] <- as.character(meta$CellType[which(meta$batchlb == batches[j])][a3[,2]])
  print(length(which(a3[,1] == a3[, 2])))
  print(length(which(a3[,1] == a3[, 2]))/nrow(a2))
  print(summary(as.factor(a3[,1])))
  print(sum(summary(as.factor(a3[,1]))))
  print(summary(as.factor(a3[,2])))
  print(sum(summary(as.factor(a3[,2]))))
  a1[,1] <- as.character(meta$CellType[a1[,1]])
  a1[,2] <- as.character(meta$CellType[a1[,2]])
  print(length(which(a1[,1] == a1[, 2])))
  print(length(which(a1[,1] == a1[, 2]))/nrow(a1))
  print(summary(as.factor(a1[,1])))
  print(sum(summary(as.factor(a1[,1]))))
  print(summary(as.factor(a1[,2])))
  print(sum(summary(as.factor(a1[,2]))))
}

library(gridExtra)
library(grid)
library(ggplot2)

grid_arrange_shared_legend <- function(plot_list, ncol = length(plot_list), nrow = 1, position = c("bottom", "right")) {
  
  plots <- plot_list
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

grid_arrange_shared_multi_legends <- function(p_list1, p_list2, position = c("bottom", "right")) {
  plots1 <- p_list1
  plots2 <- p_list2
  
  position <- match.arg(position)
  
  g1 <- ggplotGrob(plots1[[1]] + theme(legend.position = position))$grobs
  legend1 <- g1[[which(sapply(g1, function(x) x$name) == "guide-box")]]
  lheight1 <- sum(legend1$height)
  lwidth1 <- sum(legend1$width)
  gl1 <- lapply(plots1, function(x) x + theme(legend.position="none"))
  gl1 <- c(gl1, ncol = length(gl1), nrow = 1)
  
  g2 <- ggplotGrob(plots2[[1]] + theme(legend.position = position))$grobs
  legend2 <- g2[[which(sapply(g2, function(x) x$name) == "guide-box")]]
  lheight2 <- sum(legend2$height)
  lwidth2 <- sum(legend2$width)
  gl2 <- lapply(plots2, function(x) x + theme(legend.position="none"))
  gl2 <- c(gl2, ncol = length(gl2), nrow = 1)
  
  g_blank <- ggplotGrob(p_list1[[1]] + theme(legend.position = NULL))$grobs
  g_blank <- g_blank[[1]]
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(
                       arrangeGrob(do.call(arrangeGrob, gl1),legend1,g_blank,ncol = 1,heights = unit.c(unit(1, "npc") - lheight2, lheight1, lheight2 - lheight1)),
                       arrangeGrob(do.call(arrangeGrob, gl2),legend2,ncol = 2,widths = unit.c(unit(1, "npc") - lheight2, lheight2)),
                       ncol = 2),
                     "right" = arrangeGrob(
                       arrangeGrob(do.call(arrangeGrob, gl1),legend1,g_blank,ncol = 3,widths = unit.c(unit(1, "npc") - lwidth2, lwidth1, lwidth2 - lwidth1)),
                       arrangeGrob(do.call(arrangeGrob, gl2),legend2,ncol = 2,widths = unit.c(unit(1, "npc") - lwidth2, lwidth2)),
                       nrow = 2)) %>% as.ggplot()
  combined          
}

library(pheatmap)
library(grid)
library(gtable)

colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)

data.random <- data.frame(matrix(data=rnorm(150), ncol=15, nrow=10))
rownames(data.random) <- letters[1:10]
colnames(data.random) <- LETTERS[1:15]

# Modified pheatmap:::heatmap_motor
heatmap_motor0 <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
                           tree_row, treeheight_col, treeheight_row, filename, width, 
                           height, breaks, color, legend, annotation_row, annotation_col, 
                           annotation_colors, annotation_legend, annotation_names_row, 
                           annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                           hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
                           gaps_col, gaps_row, labels_row, labels_col, ...) 
{
  lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
                     ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
                     treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
                     legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
                     annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                     annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
                     main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                     fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
                     gaps_col = gaps_col, ...)
  res = lo$gt
  mindim = lo$mindim
  if (!is.na(filename)) {
    if (is.na(height)) {
      height = convertHeight(gtable_height(res), "inches", valueOnly = T)
    }
    if (is.na(width)) {
      width = convertWidth(gtable_width(res), "inches", valueOnly = T)
    }
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) 
      stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
               png = function(x, ...) png(x, units = "in", res = 300, 
                                          ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                                                                             res = 300, ...), jpg = function(x, ...) jpeg(x, 
                                                                                                                          units = "in", res = 300, ...), tiff = function(x, 
                                                                                                                                                                         ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                                                                                                                                                                                   ...), bmp = function(x, ...) bmp(x, units = "in", 
                                                                                                                                                                                                                    res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
    f(filename, height = height, width = width)
    gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                       border_color = border_color, tree_col = tree_col, 
                       tree_row = tree_row, treeheight_col = treeheight_col, 
                       treeheight_row = treeheight_row, breaks = breaks, 
                       color = color, legend = legend, annotation_col = annotation_col, 
                       annotation_row = annotation_row, annotation_colors = annotation_colors, 
                       annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
                       annotation_names_col = annotation_names_col, filename = NA, 
                       main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                       fontsize_col = fontsize_col, hjust_col = hjust_col, 
                       vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
                       fontsize_number = fontsize_number, number_color = number_color, 
                       labels_row = labels_row, labels_col = labels_col, 
                       gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()
    return(gt)
  }
  if (mindim < 3) 
    border_color = NA
  if (!is.na(main)) {
    elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
                          clip = "off")
  }
  if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
    elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
    res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
  }
  if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
    elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
  }
  elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
                                fmat, fontsize_number, number_color)
  res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                        name = "matrix")
  if (length(labels_col) != 0) {
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
                hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
                ...)
    elem = do.call(pheatmap:::draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
                          name = "col_names")
  }
  if (length(labels_row) != 0) {
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
                ...)
    elem = do.call(pheatmap:::draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                          name = "row_names")
  }
  if (!pheatmap:::is.na2(annotation_col)) {
    converted_annotation = convert_annotations(annotation_col, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_col, fontsize, horizontal = T)
    res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
                          name = "col_annotation")
    if (annotation_names_col) {
      elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                                              horizontal = T)
      res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                            name = "col_annotation_names")
    }
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    converted_annotation = convert_annotations(annotation_row, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_row, fontsize, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
                          name = "row_annotation")
    if (annotation_names_row) {
      elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                                              horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                                              angle_col = angle_col)
      res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                            name = "row_annotation_names")
    }
  }
  annotation = c(annotation_col[length(annotation_col):1], 
                 annotation_row[length(annotation_row):1])
  annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
  if (length(annotation) > 0 & annotation_legend) {
    elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
                                             border_color, fontsize = fontsize, ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
                          clip = "off", name = "annotation_legend")
  }
  if (!pheatmap:::is.na2(legend)) {
    elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
                                  ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
                          clip = "off", name = "legend")
  }
  return(res)
}

# Modified pheatmap:::lo    
lo0 <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
                annotation_colors, annotation_legend, annotation_names_row, 
                annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                angle_col, gaps_row, gaps_col, ...) 
{
  if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
    if (!is.null(coln[1])) {
      t = coln
    }
    else {
      t = ""
    }
    tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
    if (annotation_names_row) {
      t = c(t, colnames(annotation_row))
      tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
    }
    longest_coln = which.max(tw)
    gp = list(fontsize = ifelse(longest_coln <= length(coln), 
                                fontsize_col, fontsize), ...)
    coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
                                                 rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
                                                                                                  "bigpts")
  }
  else {
    coln_height = unit(5, "bigpts")
  }
  if (!is.null(rown[1])) {
    t = rown
    tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
    if (annotation_names_col) {
      t = c(t, colnames(annotation_col))
      tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
    }
    longest_rown = which.max(tw)
    gp = list(fontsize = ifelse(longest_rown <= length(rown), 
                                fontsize_row, fontsize), ...)
    rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
                                               rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else {
    rown_width = unit(5, "bigpts")
  }
  gp = list(fontsize = fontsize, ...)
  if (!pheatmap:::is.na2(legend)) {
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", 
                         textGrob(as.character(names(legend))[longest_break], 
                                  gp = do.call(gpar, gp)))
    title_length = unit(1.1, "grobwidth", textGrob("Scale", 
                                                   gp = gpar(fontface = "bold", ...)))
    legend_width = unit(12, "bigpts") + longest_break * 1.2
    legend_width = max(title_length, legend_width)
  }
  else {
    legend_width = unit(0, "bigpts")
  }
  if (is.na(main)) {
    main_height = unit(0, "npc")
  }
  else {
    main_height = unit(1.5, "grobheight", textGrob(main, 
                                                   gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }
  textheight = unit(fontsize, "bigpts")
  if (!pheatmap:::is.na2(annotation_col)) {
    annot_col_height = ncol(annotation_col) * (textheight + 
                                                 unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
    annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_col_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_col_height = unit(0, "bigpts")
    annot_col_legend_width = unit(0, "bigpts")
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    annot_row_width = ncol(annotation_row) * (textheight + 
                                                unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
    annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_row_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_row_width = unit(0, "bigpts")
    annot_row_legend_width = unit(0, "bigpts")
  }
  annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
  treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
                                                         "bigpts")
  treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
                                                         "bigpts")
  if (is.na(cellwidth)) {
    mat_width = unit(1, "npc") - rown_width - legend_width - 
      treeheight_row - annot_row_width - annot_legend_width
  }
  else {
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
      unit(4, "bigpts")
  }
  if (is.na(cellheight)) {
    mat_height = unit(1, "npc") - main_height - coln_height - 
      treeheight_col - annot_col_height
  }
  else {
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
      unit(4, "bigpts")
  }
  gt = gtable(widths = unit.c(treeheight_row, rown_width,  
                              mat_width, treeheight_row, legend_width, annot_legend_width), 
              heights = unit.c(main_height, treeheight_col, annot_col_height, 
                               mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
                                                                                    gp)))
  cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
                                                         "bigpts")), "bigpts", valueOnly = T)/ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
                                                           "bigpts")), "bigpts", valueOnly = T)/nrow
  mindim = min(cw, ch)
  res = list(gt = gt, mindim = mindim)
  return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames0 <- function (rown, gaps, ...) 
{
  coord = pheatmap:::find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
                 hjust = 1, gp = gpar(...))
  return(res)
}

#assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
#assignInNamespace(x="lo", value=lo, ns="pheatmap")
#assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap")