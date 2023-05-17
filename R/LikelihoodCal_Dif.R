#' Calculate expression likelihood with diffusion loss
#'
#' This function computes the expression likelihood of a tree based on cell gene expression data
#' @param tree input tree
#' @param muts lineage barcode matrix
#' @param prob_t transition probabilities between cells
#' @param cell_state_labels cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param newick_lookup lookup table for internal nodes
#' @param lambda1 hyperparameter for asymmetric division likelihood
#' @param lambda2 hyperparameter for neighbor distance likelihood
#' @import ape
#' @import destiny
#' @export
LikelihoodCal_Dif <- function(tree,muts,prob_t,cell_state_labels,state_lineages,newick_lookup,lambda1 = 10,lambda2 = 1){
  #browser()
  N_char <- dim(muts)[2]
  #Get the internal node ids
  nodes_internal <- unique(tree$edge[,1])
  #order the internal nodes
  nodes_internal <- sort(nodes_internal, decreasing = TRUE)

  #extend mutation table to include the internal nodes
  muts <- rbind(muts,matrix(0,length(nodes_internal),N_char))
  #cell_state_labels <- c(cell_state_labels,rep("0",length(nodes_internal)))
  internal_lookup <- data.frame(node = numeric(),index = numeric())

  #initialize transition frequencies
  max_step <- max(sapply(state_lineages,length))
  state_transitions <- data.frame(step = c(c(0:max_step)), freq = rep(0,max_step+1))
  N_violations <- 0
  total_dist <- 0
  allele_unique <- unique(unlist(muts))
  mutation_transitions <- data.frame(allele = allele_unique,freq = rep(0,length(allele_unique)))
  total_muts <- 0

  #determine ancestor barcode based on child barcodes
  diffuse_loss <- 0
  ad_loss <- 0
  p_a <- 0.8
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
    barcodes_children <- muts[children,]
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

    if (length(children) > 1){
      for (i in 1:N_char){
        characters <- unique(barcodes_children[,i])
        if (length(characters) == 1){
          muts[node,i] <- characters
        }else{
          muts[node,i] <- "0"
          #browser()
          mutation_transitions[mutation_transitions$allele == characters[1],2] <- mutation_transitions[mutation_transitions$allele == characters[1],2] + 1
          mutation_transitions[mutation_transitions$allele == characters[2],2] <- mutation_transitions[mutation_transitions$allele == characters[2],2] + 1
        }
        total_muts <- total_muts + length(setdiff(characters, muts[node,i]))
      }
    }

    cell_state_labels <- c(cell_state_labels,state_parent)
    temp <- data.frame(node = node,index = length(cell_state_labels))
    internal_lookup <- rbind(internal_lookup,temp)

    #print(children)
    #print(leaf_counts)
    if (length(children) > 1){
      if((cell_state_children[1] == cell_state_children[2])&(leaf_counts == 2)){
        diffuse_loss <- diffuse_loss + log2(prob_t[children[1],children[2]])
      }
    }

    if (length(children) > 1){
      if(cell_state_children[1] == cell_state_children[2]){
        ad_loss <- ad_loss + log2(1-p_a)
      } else {
        ad_loss <- ad_loss + log2(p_a)
      }
    }

    if (cell_state_labels[length(cell_state_labels)] == "0"){
      N_violations <- N_violations + 2
    }
  }

  mutation_transitions <- mutation_transitions[mutation_transitions$freq > 0,]
  state_transitions <- state_transitions[state_transitions$freq > 0,]
  max_step <- length(state_transitions$freq)

  l_barcode <- sum(mutation_transitions$freq * log2(mutation_transitions$freq/sum(mutation_transitions$freq)))
  if (max_step > 0){
    transit_prob <- exp(-1:-max_step)/sum(exp(-1:-max_step))
    l_expression <- sum(state_transitions$freq * transit_prob) - N_violations*50 + lambda1 *  ad_loss + lambda2 * diffuse_loss
    #l_expression <- - N_violations*50 + lambda1 *  ad_loss + lambda2 * diffuse_loss
  }else {
    l_expression <- - N_violations*50 + lambda1 * ad_loss + lambda2 * diffuse_loss
  }
  #browser()
  #l_expression <- sum(state_transitions$freq * log2(state_transitions$freq/sum(state_transitions$freq))) - N_violations*50 + diffuse_loss
  return(l_expression + l_barcode)

}
