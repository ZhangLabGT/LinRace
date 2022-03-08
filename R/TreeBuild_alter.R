DivideMut <- function(X){
  mut_table_str <- c()
  for (i in 1:nrow(X)){
    barcode <- X[i,]
    cs <- paste(X[i,], collapse = '|')
    mut_table_str <- c(mut_table_str,cs)
  }
  dt <- as.data.table(mut_table_str)[, list(list(.I)), by = mut_table_str]
  colnames(dt) <- c("barcode","cellids")

  #obtain tree backbone using nj()
  X_unique <- X[!duplicated(X), ]
  l <- dim(X_unique)[1]
  rownames(X_unique) <- 1:l
  sim_tree <- list()
  for (j in 1:l){
    sim_tree[[j]] <- as.character(X_unique[j,])
  }

  names(sim_tree) <- rownames(X_unique)
  #nmstrings <- c('0','-',c(1:100))
  nmstrings <- unique(unlist(sim_tree))
  sim_data <- phyDat(sim_tree,type = 'USER',levels = nmstrings)

  #dist_wh2 <- WH(sim_data, InfoW, dropout=TRUE)
  dist_h <- dist.hamming(sim_data)
  Treenj <- nj(dist_h)
  tree_backbone <- multi2di(Treenj)
  tree_backbone$edge.length <- rep(1, length(tree_backbone$edge.length))

  return(list(tree_backbone,dt))
}

#' LinRace main function(alter): asymmetric division based Neighbor Joining
#'
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @import data.table
#' @import phangorn
#' @export
NJ_asym_alter <- function(muts,states,state_lineages,max_Iter = 200){
  labels <- rownames(muts)
  returnList <- DivideMut(muts)

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]

  subtree_list <- list()
  for (i in 1:nrow(dt)){
    cellids <- unlist(dt$cellids[i])
    muts_sub <- muts[cellids,]
    states_sub <- states[cellids]
    labels_sub <- labels[cellids]

      res <- FindExpTree(states_sub,labels = labels_sub,state_lineages,muts,maxIter = max_Iter)
    subtree_opt <- res[[1]]
    subtree_opt$name <- as.character(i)
    subtree_list[[length(subtree_list)+1]] <- subtree_opt
  }

  tree_final <- ConstructTree(tree_backbone,subtree_list)
}
