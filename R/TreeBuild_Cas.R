DivideMut_Cas <- function(X){
  mut_table_str <- c()
  for (i in 1:nrow(X)){
    barcode <- X[i,]
    cs <- paste(X[i,], collapse = '|')
    mut_table_str <- c(mut_table_str,cs)
  }
  dt <- as.data.table(mut_table_str)[, list(list(.I)), by = mut_table_str]
  colnames(dt) <- c("barcode","cellids")

  #prob_table <- table(unlist(c(X)))
  #prob_table <- prob_table[rownames(prob_table)!="-"]
  #prob_table <- prob_table[rownames(prob_table)!="0"]
  #prob_table <- log2(prob_table+1)
  #prob_table <- prob_table/sum(prob_table)
  #prob_table <- log2(prob_table)

  #obtain tree backbone using nj()
  X_unique <- X[!duplicated(X), ]
  prob_table <- table(unlist(c(X_unique)))
  prob_table <- prob_table[rownames(prob_table)!="-"]
  prob_table <- prob_table[rownames(prob_table)!="0"]
  #prob_table <- log2(prob_table+1)
  prob_table <- prob_table/sum(prob_table)
  prob_table <- log2(prob_table)


  l <- dim(X_unique)[1]
  rownames(X_unique) <- 1:l

  newick <- Cas_greedy(X_unique,prob_table)
  tree_backbone <- read.tree(text=paste(newick,";"))
  tree_backbone$edge.length <- rep(1, length(tree_backbone$edge.length))

  return(list(tree_backbone,dt))
}

#' Generate unique barcode data for cas hybrid runs
#'
#' @param mut_dir directory to lineage barcode matrix

#' @export
Generate_Mut_Unique <- function(mut_dir) {
  X <- read.csv(file = mut_dir,row.names = 1,stringsAsFactors = FALSE)

  X_unique <- X[!duplicated(X),]
  l <- dim(X_unique)[1]
  rownames(X_unique) <- 1:l

  mut_unique_dir <- paste0(gsub(".csv","",mut_dir),"_unique.csv")
  write.csv(X_unique, mut_unique_dir)
}

Cas_greedy <- function(X,prob_table){
  if (nrow(X) == 1){return(rownames(X))}
  #if (nrow(X) == 1){return(paste(X, collapse = '|'))}
  freq_table <- plyr::alply(X, 2, table)
  min_loss <- 0
  min_target <- 0
  min_mut <- '0'

  for (i in 1:length(freq_table)){
    freq_table_target <- freq_table[[i]]
    #if (length(rownames(freq_table_target)) < 1){
    #  browser()
    #}
    rownames(freq_table_target) <- gsub(" ","",rownames(freq_table_target))
    freq_table_target <- freq_table_target[rownames(freq_table_target)!="-"]

    if (length(freq_table_target) > 1){
      freq_table_target <- freq_table_target[rownames(freq_table_target)!="0"]
      for (mutation in rownames(freq_table_target)){
	loss <- prob_table[[mutation]]*freq_table_target[[mutation]]
        if (loss < min_loss){
          min_loss <- loss
          min_target <- i
          min_mut <- mutation
        }

      }

    }
  }
  #browser()
  #print(c(max_target,max_mut))
  if(min_target == 0){
    I <- c(1)
    E <- 2:nrow(X)

  }else {
    I <- c()
    E <- c()
    D <- c()
    for (j in 1:nrow(X)){
      barcode <- X[j,]
      if (barcode[min_target] == min_mut){
        I <- c(I,j)
      }else if (barcode[min_target] == "-"){
        D <- c(D,j)
      } else{
        E <- c(E,j)
      }
    }

    freq_table_I <- apply(X[I,], 2, table)
    freq_table_E <- apply(X[E,], 2, table)
    for (d in D){
      freq_I <- 0
      freq_E <- 0
      barcode <- X[d,]
      for (k in 1:length(barcode)){
        tryCatch({freq_I <- freq_I + freq_table_I[[k]][[barcode[[k]]]]},error = function(e){})
        tryCatch({freq_E <- freq_E + freq_table_E[[k]][[barcode[[k]]]]},error = function(e){})
      }

      if (freq_I > freq_E){
        I <- c(I,d)
      } else{
        E <- c(E,d)
      }
    }
  }
  t_I <- Cas_greedy(X[I,],prob_table)
  t_E <- Cas_greedy(X[E,],prob_table)

  newick = sprintf("(%s:%f, %s:%f)",t_I, 1, t_E, 1)

  return(newick)
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
NJ_asym_Cas_Dif <- function(muts,states,counts,state_lineages,max_Iter = 200){
  labels <- rownames(muts)
  returnList <- DivideMut_Cas(muts)

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]

  subtree_list <- list()

  dm <- DiffusionMap(data=log2(counts+1),n_pcs = 30)
  prob_t <- dm@transitions

  for (i in 1:nrow(dt)){
    cellids <- unlist(dt$cellids[i])
    muts_sub <- muts[cellids,]
    counts_sub <- counts[cellids,]
    states_sub <- states[cellids]
    labels_sub <- labels[cellids]
    if (length(labels_sub)>1){
      #res <- FindExpTree_Dif(states,labels = labels_sub,state_lineages,muts,prob_t,maxIter = max_Iter)
      res <- FindExpTree_Dif_Newick(states, labels = labels_sub, state_lineages, muts, prob_t, maxIter = max_Iter, lambda1 = 1,lambda2 = 1)
      subtree_opt <- res[[1]]
      subtree_opt$name <- as.character(i)
      subtree_list[[length(subtree_list)+1]] <- subtree_opt
    }else {
      tree_backbone$tip.label[tree_backbone$tip.label==as.character(i)] <- labels_sub
    }
  }

  tree_final <- ConstructTree(tree_backbone,subtree_list)
}



#' LinRace main function with Cassiopeia backbone: asymmetric division based Cassiopeia-hybrid
#'
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param counts gene expression data of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @import data.table
#' @import phangorn
#' @export
NJ_asym_Cas_hybrid <- function(muts,states,counts,state_lineages,max_Iter = 200,tree_backbone){
  labels <- rownames(muts)
  returnList <- DivideMut_Cas(muts)

  dt <- returnList[[2]]

  subtree_list <- list()

  dm <- DiffusionMap(data=log2(counts+1),n_pcs = 30)
  prob_t <- dm@transitions

  for (i in 1:nrow(dt)){
    cellids <- unlist(dt$cellids[i])
    muts_sub <- muts[cellids,]
    counts_sub <- counts[cellids,]
    states_sub <- states[cellids]
    labels_sub <- labels[cellids]
    if (length(labels_sub)>1){
      #res <- FindExpTree_Dif(states,labels = labels_sub,state_lineages,muts,prob_t,maxIter = max_Iter)
      res <- FindExpTree_Dif_Newick(states, labels = labels_sub, state_lineages, muts, prob_t, maxIter = max_Iter, lambda1 = 1,lambda2 = 1)
      subtree_opt <- res[[1]]
      subtree_opt$name <- as.character(i)
      subtree_list[[length(subtree_list)+1]] <- subtree_opt
    }else {
      tree_backbone$tip.label[tree_backbone$tip.label==as.character(i)] <- labels_sub
    }
  }

  tree_final <- ConstructTree(tree_backbone,subtree_list)
}
