#' Calculate expression likelihood with diffusion loss
#'
#' This function computes the expression likelihood of a tree based on cell gene expression data
#' @param tree input tree
#' @param counts gene expression data
#' @param cell_state_labels cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param newick_lookup lookup table for internal nodes
#' @import ape
#' @import destiny
#' @export
LikelihoodCal_Dif <- function(tree,counts,cell_state_labels,newick_lookup){

  dm <- DiffusionMap(data=log2(counts+1),n_pcs = min(30,nrow(counts)))
  prob_t <- dm@transitions

  #browser()
  #Get the internal node ids
  nodes_internal <- unique(tree$edge[,1])
  #order the internal nodes
  nodes_internal <- sort(nodes_internal, decreasing = TRUE)

  #extend mutation table to include the internal nodes
  #cell_state_labels <- c(cell_state_labels,rep("0",length(nodes_internal)))
  internal_lookup <- data.frame(node = numeric(),index = numeric())

  #initialize transition frequencies
  max_step <- max(sapply(state_lineages,length))
  state_transitions <- data.frame(step = c(c(0:max_step)), freq = rep(0,max_step+1))
  N_violations <- 0
  total_dist <- 0

  #determine ancestor barcode based on child barcodes
  diffuse_loss <- 0
  for (node in nodes_internal){
    children_index <- tree$edge[tree$edge[,1] == node,2]
    children <- c()
    leaf_counts <- 0
    for (cell_id in children_index){
      if (cell_id > length(tree$tip.label)){
        cell_row <- internal_lookup$index[internal_lookup$node == cell_id]
      }else {
        leaf_counts <- leaf_counts + 1
        cell_id <- tree$tip.label[cell_id]
        if(substr(cell_id,1,1)=="c"){
          cell_row <- strtoi(substr(cell_id,6,nchar(cell_id)))
        }else{
          cell_row <- newick_lookup$index[newick_lookup$newick==cell_id]
        }
      }
      children <- c(children, cell_row)
    }
    cell_state_children <- cell_state_labels[children]

    state_parent <- "0"
    for (lineage in state_lineages){
      match_child1 <- NA
      match_child2 <- NA
      match_child1 <- match(cell_state_children[1],lineage)
      match_child2 <- match(cell_state_children[2],lineage)
      if (is.na(match_child1) || is.na(match_child2)) next
      state_parent <- lineage[min(match_child1,match_child2)]
      state_dist <- abs(match_child1 - match_child2)
      state_transitions[state_transitions$step == state_dist,2] <- state_transitions[state_transitions$step == state_dist,2] + 1
      state_transitions[state_transitions$step == 0,2] <- state_transitions[state_transitions$step == 0,2] + 1
      break
    }
    cell_state_labels <- c(cell_state_labels,state_parent)
    temp <- data.frame(node = node,index = length(cell_state_labels))
    internal_lookup <- rbind(internal_lookup,temp)

    if((cell_state_children[1] == cell_state_children[2])&(leaf_counts == 2)){
      diffuse_loss <- diffuse_loss + log2(prob_t[children[1]],prob_t[children[2]])
    }

    if (cell_state_labels[length(cell_state_labels)] == "0"){
      N_violations <- N_violations + 2
    }
  }
  state_transitions <- state_transitions[state_transitions$freq > 0,]
  l_expression <- sum(state_transitions$freq * log2(state_transitions$freq/sum(state_transitions$freq))) - N_violations*50
  return(l_expression)
}
