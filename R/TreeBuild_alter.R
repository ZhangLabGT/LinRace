BarToMut <- function(X){
  mutations <- as.character(unique(unlist(X)))
  mutations <- mutations[mutations != "0"]

  Mut <- as.data.frame(matrix(0,nrow = nrow(X), ncol = length(mutations)))
  colnames(Mut) <- mutations
  rownames(Mut) <- rownames(X)
  for (i in 1:nrow(X)){
    barcode <- X[i,]
    for (ch in barcode){
      if (ch != "0"){
        Mut[i,which(mutations==ch)] <- 1
      }
    }
  }

  return(Mut)
}

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
  #tree_backbone$tip.label <- 1

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
    if (length(labels_sub)==1){
      res <- FindExpTree(states_sub,labels = labels_sub,state_lineages,muts,maxIter = max_Iter)
      subtree_opt <- res[[1]]
      subtree_opt$name <- as.character(i)
      subtree_list[[length(subtree_list)+1]] <- subtree_opt
    }
  }

  tree_final <- ConstructTree(tree_backbone,subtree_list)
}


#' LinRace main function(Dif): asymmetric division based Neighbor Joining
#'
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param counts gene expression data of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @import data.table
#' @import phangorn
#' @export
NJ_asym_Dif <- function(muts,states,counts,state_lineages,max_Iter = 200){
  labels <- rownames(muts)
  X <- muts
  #X <- BarToMut(muts)
  returnList <- DivideMut(X)

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]

  subtree_list <- list()

  dm <- DiffusionMap(data=log2(counts+1),n_pcs = 30)
  prob_t <- dm@transitions

  #names_subtree <- c()

  df_pairs <- FindPairs(tree_backbone,dt)
  res <- MergeDT(tree_backbone,df_pairs,dt)
  tree_backbone <- res[[1]]
  dt <- res[[2]]


  for (i in 1:nrow(dt)){
    if (!dt$merged[i]){
      cellids <- unlist(dt$cellids[i])
      muts_sub <- muts[cellids,]
      counts_sub <- counts[cellids,]
      states_sub <- states[cellids]
      labels_sub <- labels[cellids]

      if (length(labels_sub)>1){
        res <- FindExpTree_Dif_Newick(states,labels = labels_sub,state_lineages,muts,prob_t,maxIter = max_Iter)
        subtree_opt <- res[[1]]
        subtree_opt$name <- as.character(i)
        #names_subtree <- c(names_subtree,as.character(i))
        subtree_list[[length(subtree_list)+1]] <- subtree_opt
      }else {
        tree_backbone$tip.label[tree_backbone$tip.label==as.character(i)] <- labels_sub
      }
    }
  }
  #names(subtree_list) <- names_subtree

  tree_final <- ConstructTree(tree_backbone,subtree_list)
}

#' Finding the best tree structure based on local search
#'
#' @param states the states of the gene expressions of cells
#' @param labels the labels of cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param muts lineage barcode data of cells
#' @param prob_f transition probabilities between cells
#' @param newick_lookup lookup table for internal nodes
#' @param maxIter maximum number of iterations of local search
#' @import castor
FindExpTree_Dif <- function(states,labels,state_lineages,muts,prob_f,newick_lookup = NULL, maxIter = 200){
  #browser()
  leafs <- c()
  for (i in 1:length(labels)){
    label <- labels[i]
    n <- nchar(label)
    if(substr(label, n, n) == ';'){
      labels[i] <- substr(label,1,n-1)
    }else{
      leafs <- c(leafs,as.numeric(substr(label,6,n)))
    }
  }

  if (length(labels)>2){
    prob_t_sub <- prob_f[leafs,leafs]
    rownames(prob_t_sub) <- labels
    colnames(prob_t_sub) <- labels

    Treenj <- nj(as.dist(exp(-prob_t_sub)))
    #dist_h <- distH(exp(-prob_t_sub))
    #Treenj <- nj(dist_h)
    tree_init <- multi2di(Treenj)
  }else {
    tree_init <- rtree(length(labels), rooted = TRUE, tip.label = labels)
  }

  tree <- tree_init
  #tree <- rtree(length(labels), rooted = TRUE, tip.label = labels)
  tree$edge.length <- rep(1, length(tree$edge.length))
  edges <- tree$edge
  #maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  max_likelihood <- LikelihoodCal_Dif(tree,prob_f,states,state_lineages,newick_lookup)

  likelihood_curve <-c()
  seeds <- runif(10000,1,99999)
  best_tree <- tree
  best_tree_list <- list()
  maxl_list <- c()
  ptm <- proc.time()
  LikelihoodCal_time <- 0
  TreeLC_time <- 0
  for (i in 1:maxIter){
    c_time <- proc.time()

    tree_new <- TreeLC2(tree)
    #tree_new <-  rNNI(tree,1,1)
    TreeLC_time <- TreeLC_time + proc.time()-c_time
    #cl <- LikelihoodCal(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
    c_time <- proc.time()
    likelihood_new <- LikelihoodCal_Dif(tree_new,prob_f,states,state_lineages,newick_lookup)
    #cl_list <- LikelihoodCal2(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
    LikelihoodCal_time <- LikelihoodCal_time + proc.time()-c_time
    likelihood_curve <- c(likelihood_curve,max_likelihood)
    if (likelihood_new > max_likelihood){
      max_likelihood <- likelihood_new
      tree <- tree_new
      best_tree <- tree
    }
    if (i %% 100 == 0){
      likelihood_check <- sprintf("After %g iterations, the best likelihood is %g.",
                                  i, max_likelihood)
      print(likelihood_check)
    }
    if (i>100){
      if (length(unique(likelihood_curve[(i-20):i])) == 1){
        #record the local optima, and
        best_tree_list[[length(best_tree_list)+1]] <- best_tree
        maxl_list <- c(maxl_list,max_likelihood)

        tree <- tree_init
        #tree <- rtree(length(labels), rooted = TRUE, tip.label = labels)
        tree$edge.length <- rep(1, length(tree$edge.length))
        #maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
        max_likelihood <- LikelihoodCal_Dif(tree,prob_f,states,state_lineages,newick_lookup)
      }
    }
  }
  #browser()
  #return the best tree found so far
  if (length(maxl_list) > 0){
    best_tree <- best_tree_list[[which.max(maxl_list)]]
  }else{
    best_tree <- tree
  }
  ances_res <- AncesInfer_sub(tree,muts,states,state_lineages,newick_lookup)
  root_barcode <- ances_res[[1]][dim(ances_res[[1]])[1],]
  root_state_label <- ances_res[[2]][length(ances_res[[2]])]

  return(list(best_tree,root_state_label,root_barcode))
}


FindExpTree_Dif_Newick <- function(states,labels,state_lineages,muts,prob_f,newick_lookup = NULL, maxIter = 200){
  #browser()
  leafs <- c()
  for (i in 1:length(labels)){
    label <- labels[i]
    n <- nchar(label)
    if(substr(label, n, n) == ';'){
      labels[i] <- substr(label,1,n-1)
    }else{
      leafs <- c(leafs,as.numeric(substr(label,6,n)))
    }
  }
  N <- length(unique(states[leafs]))
  maxIter <- min(c(3 * N * N, 1000))
  sliding <- min(100, N * N)

  if (length(labels)>2){
    prob_t_sub <- prob_f[leafs,leafs]
    rownames(prob_t_sub) <- labels
    colnames(prob_t_sub) <- labels

    Treenj <- nj(as.dist(exp(-prob_t_sub)))
    #dist_h <- distH(exp(-prob_t_sub))
    #Treenj <- nj(dist_h)
    tree_init <- multi2di(Treenj)
  }else {
    tree_init <- rtree(length(labels), rooted = TRUE, tip.label = labels)
  }

  tree <- tree_init
  #tree <- rtree(length(labels), rooted = TRUE, tip.label = labels)
  tree$edge.length <- rep(1, length(tree$edge.length))
  edges <- tree$edge
  #maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  max_likelihood <- LikelihoodCal_Dif(tree,prob_f,states,state_lineages,newick_lookup)

  likelihood_curve <-c()
  seeds <- runif(10000,1,99999)
  best_tree <- tree
  best_tree_list <- list()
  maxl_list <- c()
  ptm <- proc.time()
  LikelihoodCal_time <- 0
  TreeLC_time <- 0

  # newick swapping for faster speeds
  tree_newick <- stringi::stri_replace_all_fixed(NewickTree(tree), " ", "_")
  tree_newick <- stringi::stri_replace_all_fixed(tree_newick, ":1", "")
  num_leaves <- length(tree$tip.label)
  # print("starting iterations")

  for (i in 1:maxIter){
    c_time <- proc.time()

    #tree_new <- TreeLC2(tree)
    tree_new_newick <- TreeLC2Newick(tree_newick, num_leaves)
    #tree_new <-  rNNI(tree,1,1)
    TreeLC_time <- TreeLC_time + proc.time()-c_time
    #cl <- LikelihoodCal(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
    c_time <- proc.time()

    #likelihood_new <- LikelihoodCal_Dif(tree_new,prob_f,states,state_lineages,newick_lookup)
    # do likelihood calculation with new newick tree
    likelihood_new <- LikelihoodCal_Dif(ape::read.tree(text=tree_new_newick),prob_f,states,state_lineages,newick_lookup)

    #cl_list <- LikelihoodCal2(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
    LikelihoodCal_time <- LikelihoodCal_time + proc.time()-c_time
    likelihood_curve <- c(likelihood_curve,max_likelihood)
    if (likelihood_new > max_likelihood){
      max_likelihood <- likelihood_new
      # store tree_newick instead of tree
      tree_newick <- tree_new_newick
      best_tree <- ape::read.tree(text=tree_newick)
      best_tree$edge.length <- rep(1, length(tree$edge.length))
    }
    if (i %% 100 == 0){
      likelihood_check <- sprintf("After %g iterations, the best likelihood is %g.",
                                  i, max_likelihood)
      print(likelihood_check)
    }
    if (i > sliding){
      if (length(unique(likelihood_curve[(i-sliding):i])) == 1){
        #record the local optima, and
        best_tree_list[[length(best_tree_list)+1]] <- best_tree
        maxl_list <- c(maxl_list,max_likelihood)

        tree <- tree_init
	#tree <- rtree(length(labels), rooted = TRUE, tip.label = labels)
        tree$edge.length <- rep(1, length(tree$edge.length))
        #maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
        max_likelihood <- LikelihoodCal_Dif(tree,prob_f,states,state_lineages,newick_lookup)
        tree_newick <- stringi::stri_replace_all_fixed(NewickTree(tree), " ", "_")
        tree_newick <- stringi::stri_replace_all_fixed(tree_newick, ":1", "")
      }
    }
  }
  #browser()
  #return the best tree found so far
  if (length(maxl_list) > 0){
    best_tree <- best_tree_list[[which.max(maxl_list)]]
  }else{
    best_tree <- tree
  }

  ances_res <- AncesInfer_sub(tree,muts,states,state_lineages,newick_lookup)
  root_barcode <- ances_res[[1]][dim(ances_res[[1]])[1],]
  root_state_label <- ances_res[[2]][length(ances_res[[2]])]

  return(list(best_tree,root_state_label,root_barcode))
}

distH <- function(muts_leaves){
  sim_tree <- list()
  for (i in 1:dim(muts_leaves)[1]){
    sim_tree[[i]] <- as.character(muts_leaves[i,])
    #names(sim_tree[[i]]) <- rownames(muts_leaves)[i]
  }
  #print(length(sim_tree))
  names(sim_tree) <- rownames(muts_leaves)
  #nmstrings <- c('0','-',c(1:100))
  nmstrings <- unique(unlist(sim_tree))

  sim_data<- phyDat(sim_tree,type = 'USER',levels = nmstrings)
  sim_data <- subset(sim_data,1:length(sim_tree))
  #print(length(sim_data))
  #dist_wh2 <- WH(sim_data, InfoW, dropout=TRUE)
  dist_h <- dist.hamming(sim_data)
  #print(dim(as.matrix(dist_h)))
  return(dist_h)
}

