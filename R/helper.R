rand_mut_generator <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

MutsTransform <- function(muts){
  interval_df <- data.frame(start = numeric(),end = numeric())
  dropout_num <- unique(rand_mut_generator(5000))
  k = 1
  for(i in 1:nrow(muts)){
    barcode <- muts[i,]
    for (j in 1:length(barcode)){
      if (j >1){
        if (barcode[j] == '-' & barcode[j-1] != '-'){
          start_reg <- j
        }
        if (barcode[j]!='-' & barcode[j-1] == '-'){
          if (nrow(merge(data.frame(start = start_reg, end = j),interval_df)) == 0){
            temp_df <- data.frame(start = start_reg,end = j)
            interval_df <- rbind(interval_df,temp_df)
            muts[i,start_reg:j] <- dropout_num[k]
            k <- k + 1
          }
          else{
            index <- which(interval_df$start == start_reg & interval_df$end == j)[1]
            muts[i,start_reg:j] <- dropout_num[index]
          }
        }
      } else if (j == 1){
        if (barcode[j] == '-'){
          start_reg <- j
        }
      } 
      if (j == length(barcode)){
        if (barcode[j]=='-' & barcode[j-1] == '-'){
          if (nrow(merge(data.frame(start = start_reg, end = j),interval_df)) == 0){
            temp_df <- data.frame(start = start_reg,end = j)
            interval_df <- rbind(interval_df,temp_df)
            muts[i,start_reg:j] <- dropout_num[k]
            k <- k + 1
          }
          else{
            index <- which(interval_df$start == start_reg & interval_df$end == j)[1]
            muts[i,start_reg:j] <- dropout_num[index]
          }
        }
      }
    }
  }
  return(muts)
}

CellStateCal <- function(true_counts,meta){
  #set.seed(1)
  counts_t <- t(true_counts[[1]])
  rownames(counts_t) <- paste("cell",meta$cellID, sep = "_")
  colnames(counts_t) <- paste("gene",seq(1:ncol(counts_t)),sep = "_")
  
  custom.config <- umap.defaults
  custom.config$min_dist <- 0.3
  pca_res <- prcomp(counts_t,scale. = TRUE)
  
  cluster_id <- kmeans(pca_res$x[,1:20], centers = 7)$cluster
  meta$cluster_kmeans <- cluster_id
  
  
  #umap_res_leaves <- umap(pca_res$x[,1:50],custom.config)
  plot(pca_res$x[,1:2], col = meta$cluster, pch=16, asp = 1)
  #legend(-15, 10,legend=unique(meta$cluster),col=unique(meta$cluster))
  
  #sce <- SingleCellExperiment(assays = List(counts = t(counts_t)))
  #sce <- slingshot(pca_res$x[,1:50], clusterLabels = meta$cluster_kmeans,start.clus = )
      sce <- slingshot(pca_res$x[,1:20], clusterLabels = meta$cluster_kmeans)
  
  plot(pca_res$x[,1:2], col = meta$cluster_kmeans, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
  #legend(-15, 10,legend=unique(meta$cluster_kmeans),col=unique(meta$cluster_kmeans))
  
  returnlist <- list(sce@lineages,meta)
  return(returnlist)
  #dataset <- wrap_expression(counts = counts_t,expression = log2(counts_t+1))
  #groups <- cbind(rownames(counts_t),meta$cluster_kmeans)
  #colnames(groups) <- c("cell_id", "group_id")
  #dataset <- add_prior_information(dataset = dataset, groups_id = as.data.frame(groups), start_id = paste("cell",sample(meta[meta$cluster == 6,5],1),sep = "_"))
  #dataset <- add_dimred(dataset = dataset, umap_res_leaves$layout)
  
  #trajectory <- infer_trajectory(dataset,method,verbose = TRUE, give_priors = c("start_id","group_id"))
  
  #p <- plot_dimred(trajectory, dimred = dataset$dimred, grouping = dataset$prior_information$groups_id)
  #p
}

