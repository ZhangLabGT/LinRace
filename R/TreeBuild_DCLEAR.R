#' LinRace main function with DCLEAR: asymmetric division based DCLEAR-kmer
#'
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param counts gene expression data of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @param lambda1 hyperparameter for asymmetric division likelihood
#' @param lambda2 hyperparameter for neighbor distance likelihood
#' @import data.table
#' @import phangorn
#' @export
NJ_asym_DCLEAR <- function(muts, states, counts, state_lineages, max_Iter = 200,lambda1 = 10,lambda2 = 1) {
  labels <- rownames(muts)
  X <- muts
  #X <- BarToMut(muts)
  returnList <- DivideMut_DCLEAR(X)

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]
  dt$merged <- FALSE

  subtree_list <- list()

  dm <- DiffusionMap(data = log2(counts + 1), n_pcs = min(nrow(counts),30))
  prob_t <- dm@transitions

  #names_subtree <- c()

  #df_pairs <- FindPairs(tree_backbone, dt)
  #res <- MergeDT(tree_backbone, df_pairs, dt)
  #tree_backbone <- res[[1]]
  #dt <- res[[2]]

  for (i in 1:nrow(dt)) {
    if (!dt$merged[i]) {
      cellids <- unlist(dt$cellids[i])
      muts_sub <- muts[cellids,]
      counts_sub <- counts[cellids,]
      states_sub <- states[cellids]
      labels_sub <- labels[cellids]

      if (length(labels_sub) > 1) {
        #res <- FindExpTree_Dif(states,labels = labels_sub,state_lineages,muts,prob_t,maxIter = max_Iter)
        res <- FindExpTree_Dif_Newick(states, labels = labels_sub, state_lineages, muts, prob_t, maxIter = max_Iter, lambda1 = lambda1,lambda2 = lambda2)

        subtree_opt <- res[[1]]
        subtree_opt$name <- as.character(i)
        #names_subtree <- c(names_subtree,as.character(i))
        subtree_list[[length(subtree_list) + 1]] <- subtree_opt
      }else {
        tree_backbone$tip.label[tree_backbone$tip.label == as.character(i)] <- labels_sub
      }
    }
  }
  #names(subtree_list) <- names_subtree

  tree_final <- ConstructTree(tree_backbone, subtree_list)
}

#' Divide muts using DCLEAR-kmer
#'
#' @param muts lineage barcode matrix
#' @import DCLEAR
#' @import phangorn
#' @export
DivideMut_DCLEAR <- function(X) {
  mut_table_str <- c()
  for (i in 1:nrow(X)) {
    barcode <- X[i,]
    cs <- paste(X[i,], collapse = '|')
    mut_table_str <- c(mut_table_str, cs)
  }
  dt <- as.data.table(mut_table_str)[, list(list(.I)), by = mut_table_str]
  colnames(dt) <- c("barcode", "cellids")


  X_unique <- X[!duplicated(X),]
  l <- dim(X_unique)[1]
  rownames(X_unique) <- 1:l

  sim_tree <- list()
  states_barcode <- unique(unlist(X_unique))
  states_mutated <- states_barcode[states_barcode != "0"]
  nstates <- length(states_mutated)
  #overflow_states <- c()
  #if (nstates > 26){
  #  overflow_states <- states_mutated[27:nstates]
  #}
  for (i in 1:dim(X_unique)[1]){
    vec <- as.character(X_unique[i,])
    #vec[vec %in% overflow_states] <- "0"
    sim_tree[[i]] <- vec
    #names(sim_tree[[i]]) <- rownames(muts_leaves)[i]
  }
  #print(length(sim_tree))
  names(sim_tree) <- rownames(X_unique)
  #nmstrings <- c('0','-',c(1:100))
  nmstrings <- unique(unlist(sim_tree))

  sim_data<- phyDat(sim_tree,type = 'USER',levels = nmstrings)
  sim_data <- subset(sim_data,1:length(sim_tree))

  dist_kmer <- dist_replacement(sim_data,reps = 20L, k = 2L)
  Treenj <- nj(dist_kmer)
  tree_backbone <- multi2di(Treenj)
  tree_backbone$edge.length <- rep(1, length(tree_backbone$edge.length))
  #tree_backbone$tip.label <- 1

  return(list(tree_backbone, dt))
}
