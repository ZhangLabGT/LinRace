#' Calculate combined likelihood
#'
#' This function computes the combined likelihood of a tree based on the lineage barcode and cell state data
#' @param tree input tree
#' @param muts lineage barcode matrix
#' @param cell_state_labels cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param transfer_prob barcode transition probabilities
#' @param state_shift_prob state shift probabilities
#' @param lambda the parameter that balances the importance of the two likelihoods
#' @import ape
#' @export
LikelihoodCal <- function(tree,muts,cell_state_labels,state_lineages,transfer_prob = NULL,state_shift_prob = NULL, lambda = 1){
  ances_res <- AncesInfer(tree,muts,cell_state_labels,state_lineages)
  muts <- ances_res[[1]]
  cell_state_labels <- ances_res[[2]]
  edges <- tree$edge
  if (is.null(transfer_prob)){
    transfer_prob <- MutationCal(edges,muts)
  }
  if (is.null(state_shift_prob)){
    state_shift_prob <- ShiftprobCal(edges,cell_state_labels,state_lineages)
  }
  cl <- 0
  cl_barcode <- 0
  cl_expression <- 0
  for (i in 1:dim(edges)[1]){
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }
    barcode_par <- muts[edge[1],]
    barcode_child <- muts[cell_child_id,]
    mut_loc <- which(barcode_par == '0')
    for (loc in mut_loc){
      mutation_child <- barcode_child[loc]
      cl_barcode <- cl_barcode + log(transfer_prob[transfer_prob$allele == mutation_child,2])
    }

    state_par <- cell_state_labels[edge[1]]
    state_child <- cell_state_labels[cell_child_id]
    state_dist <- 0
    if (state_par == '0'){
      cl_expression <- cl_expression - 50
    }else{
      for (lineage in state_lineages){
        if((state_par %in% lineage) & (state_child %in% lineage)){
          state_dist <- match(state_child,lineage)-match(state_par,lineage)
          #print(state_dist)
        }
      }
      if ((state_dist >= 0)&(state_dist <= max(state_shift_prob$step))){
        cl_expression <- cl_expression + log(state_shift_prob[state_shift_prob$step == state_dist,2])
      } else if (state_dist < 0){
        cl_expression <- cl_expression - 10000
        #cl_expression <- cl
      }
    }

  }
  cl <- cl_barcode + lambda * cl_expression
  return(list(cl = cl,l_barcode = cl_barcode,l_expression = cl_expression))
}

#' Calculate combined likelihood (alternative)
#'
#' This function computes the combined likelihood of a tree based on the lineage barcode and cell state data using a different way of calculating state shift probabilities
#' @param tree input tree
#' @param muts lineage barcode matrix
#' @param cell_state_labels cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param transfer_prob barcode transition probabilities
#' @param state_shift_prob state shift probabilities
#' @param lambda the parameter that balances the importance of the two likelihoods
#' @import ape
#' @export
LikelihoodCal2 <- function(tree,muts,cell_state_labels,state_lineages,transfer_prob = NULL,state_shift_prob = NULL, lambda = 1){
  ances_res <- AncesInfer(tree,muts,cell_state_labels,state_lineages)
  muts <- ances_res[[1]]
  cell_state_labels <- ances_res[[2]]
  edges <- tree$edge
  if (is.null(transfer_prob)){
    transfer_prob <- MutationCal(edges,muts)
  }
  if (is.null(state_shift_prob)){
    state_shift_prob <- ShiftprobCal2(edges,cell_state_labels,state_lineages)
  }
  cl <- 0
  cl_barcode <- 0
  cl_expression <- 0
  for (i in 1:dim(edges)[1]){
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }
    barcode_par <- muts[edge[1],]
    barcode_child <- muts[cell_child_id,]
    mut_loc <- which(barcode_par == '0')
    for (loc in mut_loc){
      mutation_child <- barcode_child[loc]
      cl_barcode <- cl_barcode + log(transfer_prob[transfer_prob$allele == mutation_child,2])
    }

    state_par <- cell_state_labels[edge[1]]
    state_child <- cell_state_labels[cell_child_id]
    flag <- 0
    if (state_par == '0'){
      cl_expression <- cl_expression - 50
    }else{
      for (lineage in state_lineages){
        if((state_par %in% lineage) & (state_child %in% lineage)){
          flag <- 1
        }
      }
      if (flag == 1){
        if (state_shift_prob[state_par,state_child] == 0){
          cl_expression <- cl_expression - 10000
        }else {
          cl_expression <- cl_expression + log(state_shift_prob[state_par,state_child])
        }
      }
    }

  }
  cl <- cl_barcode + lambda * cl_expression
  return(list(cl = cl,l_barcode = cl_barcode,l_expression = cl_expression))
}

Cal_LC_prob <- function(edges,cell_state_labels,state_shift_prob,state_lineages){
  LC_prob <- rep(0,dim(edges)[1])
  for (i in 1:dim(edges)[1]){
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }

    state_par <- cell_state_labels[edge[1]]
    state_child <- cell_state_labels[cell_child_id]

    flag <- 0
    if (state_par == '0'){
      LC_prob[cell_child_id] <- 1
    }else{
      for (lineage in state_lineages){
        if((state_par %in% lineage) & (state_child %in% lineage)){
          flag <- 1
          state_dist <- match(state_child,lineage)-match(state_par,lineage)
        }
      }
      if ((flag == 1) && (state_dist >= 0)){
        #LC_prob[cell_child_id] <- 1-state_shift_prob[state_par,state_child] + 0.001
        LC_prob[cell_child_id] <- 1-state_shift_prob[state_shift_prob$step == state_dist,2] + 0.001
      }else{
        LC_prob[cell_child_id] <- 1
      }
    }
  }
  LC_prob <- LC_prob/sum(LC_prob)
  return(LC_prob)
}

TreeLC2 <- function(tree,prob_internal = NULL){

  tree <- Preorder(tree)
  tree$node.label <-paste0("node", seq(1:tree$Nnode))

  if (is.null(prob_internal)){
    prob_internal <- dunif(1:(tree$Nnode + length(tree$tip.label)),1,tree$Nnode + length(tree$tip.label))
  }

  subtree_swap <- sample(tree$Nnode + length(tree$tip.label),2,prob = prob_internal)
  subtree1 <- Subtree(tree,subtree_swap[1])
  subtree2 <- Subtree(tree,subtree_swap[2])

  while (length(intersect(subtree1$tip.label, subtree2$tip.label)) > 0){
    subtree_swap <- sample(tree$Nnode + length(tree$tip.label),2)
    subtree1 <- Subtree(tree,subtree_swap[1])
    subtree2 <- Subtree(tree,subtree_swap[2])
  }

  subtree1$name <- subtree_swap[1]
  subtree1$edge.length <- rep(1, nrow(subtree1$edge))
  subtree2$name <- subtree_swap[2]
  subtree2$edge.length <- rep(1, nrow(subtree2$edge))

  root1 <- tree$edge[tree$edge[,2]==subtree1$name,1]
  root2 <- tree$edge[tree$edge[,2]==subtree2$name,1]

  if (root1 > root2) {
    reg <- subtree2
    subtree2 <- subtree1
    subtree1 <- reg
  }
  root1 <- tree$edge[tree$edge[,2]==subtree1$name,1]
  root2 <- tree$edge[tree$edge[,2]==subtree2$name,1]
  root1 <- paste0("node",root1-Ntip(tree))
  root2 <- paste0("node",root2-Ntip(tree))
  subtree1$root.edge <- 1
  subtree2$root.edge <- 1
  Ntip1 <- subtree1$Ntip
  Ntip2 <- subtree2$Ntip

  if (length(subtree1$tip.label)>1){
    tree_prune1 <- drop.tip(tree, subtree1$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  } else{
    tree_prune1 <- tree
    tree_prune1$tip.label[tree_prune1$tip.label == subtree1$tip.label] <- "NA"
  }
  if (length(subtree2$tip.label)>1){
    tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  } else{
    tree_prune2 <- tree_prune1
    tree_prune2$tip.label[tree_prune2$tip.label == subtree2$tip.label] <- "NA"
  }
  #tree_new1 <- bind.tree(tree_prune1, subtree1, where = subtree2$node.label[1]-Nnode1)
  #tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  #print(regraft_loc)
  #N_na <- sum(tree_prune2$tip.label == 'NA')
  tree_new1 <- bind.tree(tree_prune2, subtree2, where = Ntip(tree_prune2)+match(root1,tree_prune2$node.label))
  tree_new2 <- bind.tree(tree_new1, subtree1, where = Ntip(tree_new1)+match(root2,tree_new1$node.label))
  tree_final <- drop.tip(tree_new2, tree_new2$tip.label[grepl("node", tree_new2$tip.label, fixed = TRUE)], trim.internal = TRUE)
  tree_final <- drop.tip(tree_final, "NA")
  tree_final$node.label <- NULL
  tree_final$edge.length <- rep(1, length(tree_final$edge.length))
  #plot(tree_final)
  return(tree_final)
}

MutationCal <- function(edges,muts){
  allele_unique <- unique(c(muts))
  transfer_prob <- data.frame(allele = allele_unique, prob = rep(0,length(allele_unique)))
  total <- 0
  for (i in 1:dim(edges)[1]){
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }
    barcode_par <- muts[edge[1],]
    barcode_child <- muts[cell_child_id,]
    #mut_loc <- which(barcode_par == '0')
    for (loc in 1:length(barcode_par)){
      mutation_par <- barcode_par[loc]
      mutation_child <- barcode_child[loc]
      transfer_prob[transfer_prob$allele == mutation_child,2] <- transfer_prob[transfer_prob$allele == mutation_child,2] + 1
      total <- total + 1
    }
  }
  transfer_prob$prob <- transfer_prob$prob/total
  return(transfer_prob)
}

ShiftprobCal<- function(edges,cell_state_labels,state_lineages){
  max_step <- max(sapply(state_lineages,length))
  state_shift_prob <- data.frame(step = c(0:max_step), prob = rep(0,max_step+1))
  #total = 0
  for (i in 1:dim(edges)[1]){
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }
    state_par <- cell_state_labels[edge[1]]
    state_child <- cell_state_labels[cell_child_id]
    state_dist <- 0
    if (state_par!='0'){
      for (lineage in state_lineages){
        if((state_par %in% lineage) & (state_child %in% lineage)){
          state_dist <- match(state_child,lineage)-match(state_par,lineage)
          state_shift_prob[state_shift_prob$step == state_dist,2] <- state_shift_prob[state_shift_prob$step == state_dist,2] + 1
        }
      }
      #total <- total + 1
    }
  }
  state_shift_prob$prob <- state_shift_prob$prob/sum(state_shift_prob$prob)
  return(state_shift_prob)
}

ShiftprobCal2 <- function(edges,cell_state_labels,state_lineages){
  state_unique <- unique(cell_state_labels)
  state_shift_prob <- matrix(0, length(state_unique), length(state_unique))
  colnames(state_shift_prob) <- state_unique
  rownames(state_shift_prob) <- state_unique
  #transfer_prob <- data.frame(allele = allele_unique, prob = rep(0,length(allele_unique)))
  total <- 0

  for (i in 1:dim(edges)[1]){
    flag <- 0
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }
    state_par <- cell_state_labels[edge[1]]
    state_child <- cell_state_labels[cell_child_id]

    if (state_par!= '0'){
      for (lineage in state_lineages){
        if((state_par %in% lineage) & (state_child %in% lineage)){
          flag <- 1
        }
      }
      if (flag == 1){
        state_shift_prob[state_par,state_child] <- state_shift_prob[state_par,state_child] + 1
      }
      #total <- total + 1
    }
  }
  for (i in 1:ncol(state_shift_prob)){
    rowsum <- rowSums(state_shift_prob)[i]
    if (rowsum != 0){
      state_shift_prob[i,] <- state_shift_prob[i,]/rowsum
    }
  }
  return(state_shift_prob)
}


#Infer ancestor barcodes
AncesInfer <- function(tree,muts,cell_state_labels,state_lineages){
  N_char <- dim(muts)[2]
  #Get the internal node ids
  nodes_internal <- unique(tree$edge[,1])
  #order the internal nodes
  nodes_internal <- sort(nodes_internal, decreasing = TRUE)

  #extend mutation table to include the internal nodes
  muts <- rbind(muts,matrix(0,length(nodes_internal),N_char))
  cell_state_labels <- c(cell_state_labels,rep(0,length(nodes_internal)))

  #determine ancestor barcode based on child barcodes
  for (node in nodes_internal){
    children_index <- tree$edge[tree$edge[,1] == node,2]
    children <- c()
    for (cell_id in children_index){
      if (cell_id > length(tree$tip.label)){
        cell_row <- cell_id
      }else {
        cell_id <- tree$tip.label[cell_id]
        cell_row <- strtoi(substr(cell_id,6,nchar(cell_id)))
      }
      children <- c(children, cell_row)
    }
    barcodes_children <- muts[children,]
    cell_state_children <- cell_state_labels[children]

    for (lineage in state_lineages){
      if((cell_state_children[1] %in% lineage) & (cell_state_children[2] %in% lineage)){
        cell_state_labels[node] <- lineage[min(match(cell_state_children[1],lineage),match(cell_state_children[2],lineage))]
      }
    }
    if (length(children) > 1){
      for (i in 1:N_char){
        characters <- unique(barcodes_children[,i])
        if (length(characters) == 1){
          muts[node,i] <- characters
        }else{
          muts[node,] <- "0"
        }
      }
    }
  }
  return(list(muts,cell_state_labels))
}

BLikelihoodCal <- function(tree,muts,cell_state_labels,state_lineages,transfer_prob = NULL, lambda = 1, pruneQ = FALSE){
  edges <- tree$edge
  ances_res <- AncesInfer(tree,muts,cell_state_labels,state_lineages)
  muts <- ances_res[[1]]

  #Two modes of the pruning algorithm
  if (pruneQ){
    cl_expression <- ShiftprobFelsenstein(edges,cell_state_labels,state_lineages)
  }else{
    cl_expression <- ShiftprobFelsenstein(edges,cell_state_labels)
  }

  if (is.null(transfer_prob)){
    transfer_prob <- MutationCal(edges,muts)
  }

  cl <- 0
  cl_barcode <- 0
  for (i in 1:dim(edges)[1]){
    edge <- edges[i,]
    if (edge[2] > length(tree$tip.label)){
      cell_child_id <- edge[2]
    }else {
      cell_id <- tree$tip.label[edge[2]]
      cell_child_id <- strtoi(substr(cell_id,6,nchar(cell_id)))
    }
    barcode_par <- muts[edge[1],]
    barcode_child <- muts[cell_child_id,]
    mut_loc <- which(barcode_par == '0')
    for (loc in mut_loc){
      mutation_child <- barcode_child[loc]
      cl_barcode <- cl_barcode + log(transfer_prob[transfer_prob$allele == mutation_child,2])
    }

  }
  cl <- cl_barcode + lambda * log(cl_expression)
  return(list(cl = cl,l_barcode = cl_barcode,l_expression = log(cl_expression)))
}

TreeLC3 <- function(tree){

  tree <- Preorder(tree)
  tree$node.label <-paste0("node", seq(1:tree$Nnode))

  subtree_swap <- 1:(tree$Nnode + length(tree$tip.label))
  tree_list <- list()
  for (i in 1:length(subtree_swap)){
    for (j in i:length(subtree_swap)){
      subtree1 <- Subtree(tree,i)
      subtree2 <- Subtree(tree,j)
      if (length(intersect(subtree1$tip.label, subtree2$tip.label)) == 0){
        subtree1$name <- i
        subtree1$edge.length <- rep(1, nrow(subtree1$edge))
        subtree2$name <- j
        subtree2$edge.length <- rep(1, nrow(subtree2$edge))

        root1 <- tree$edge[tree$edge[,2]==subtree1$name,1]
        root2 <- tree$edge[tree$edge[,2]==subtree2$name,1]

        if (root1 > root2) {
          reg <- subtree2
          subtree2 <- subtree1
          subtree1 <- reg
        }
        root1 <- tree$edge[tree$edge[,2]==subtree1$name,1]
        root2 <- tree$edge[tree$edge[,2]==subtree2$name,1]
        root1 <- paste0("node",root1-Ntip(tree))
        root2 <- paste0("node",root2-Ntip(tree))
        subtree1$root.edge <- 1
        subtree2$root.edge <- 1
        Ntip1 <- subtree1$Ntip
        Ntip2 <- subtree2$Ntip

        if (length(subtree1$tip.label)>1){
          tree_prune1 <- drop.tip(tree, subtree1$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
        } else{
          tree_prune1 <- tree
          tree_prune1$tip.label[tree_prune1$tip.label == subtree1$tip.label] <- "NA"
        }
        if (length(subtree2$tip.label)>1){
          tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
        } else{
          tree_prune2 <- tree_prune1
          tree_prune2$tip.label[tree_prune2$tip.label == subtree2$tip.label] <- "NA"
        }
        #tree_new1 <- bind.tree(tree_prune1, subtree1, where = subtree2$node.label[1]-Nnode1)
        #tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
        #print(regraft_loc)
        #N_na <- sum(tree_prune2$tip.label == 'NA')
        tree_new1 <- bind.tree(tree_prune2, subtree2, where = Ntip(tree_prune2)+match(root1,tree_prune2$node.label))
        tree_new2 <- bind.tree(tree_new1, subtree1, where = Ntip(tree_new1)+match(root2,tree_new1$node.label))
        tree_final <- drop.tip(tree_new2, tree_new2$tip.label[grepl("node", tree_new2$tip.label, fixed = TRUE)], trim.internal = TRUE)
        tree_final <- drop.tip(tree_final, "NA")
        tree_final$node.label <- NULL
        tree_final$edge.length <- rep(1, length(tree_final$edge.length))

        tree_list[[length(tree_list)+1]] <- tree_final
      }
    }
  }
  return(tree_list)
}

NeighborTrees <- function(tree,size){
  tree_list <- list()
  for (i in 1:size){
    tree_neighbor <- TreeLC2(tree)
    tree_list[[length(tree_list)+1]] <- tree_neighbor
  }
  return(tree_list)
}
