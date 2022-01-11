LikelihoodCal_new <- function(tree,muts,cell_state_labels,state_lineages, lambda = 1, Q = 5, P = 5){
  N_char <- dim(muts)[2]
  #Get the internal node ids
  nodes_internal <- unique(tree$edge[,1])
  #order the internal nodes
  nodes_internal <- sort(nodes_internal, decreasing = TRUE)
  
  #extend mutation table to include the internal nodes
  muts <- rbind(muts,matrix(0,length(nodes_internal),N_char))
  cell_state_labels <- c(cell_state_labels,rep("0",length(nodes_internal)))
  
  #initialize transition frequencies
  max_step <- max(sapply(state_lineages,length))
  allele_unique <- unique(c(muts))
  state_transitions <- data.frame(step = c(c(0:max_step)), freq = rep(0,max_step+1))
  N_violations <- 0
  mutation_transitions <- data.frame(allele = allele_unique,freq = rep(0,length(allele_unique)))
  total_muts <- 0
  total_dist <- 0
  
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
      match_child1 <- NA
      match_child2 <- NA
      match_child1 <- match(cell_state_children[1],lineage)
      match_child2 <- match(cell_state_children[2],lineage)
      if (is.na(match_child1) || is.na(match_child2)) next
      cell_state_labels[node] <- lineage[min(match_child1,match_child2)]
      state_dist <- abs(match_child1 - match_child2)
      state_transitions[state_transitions$step == state_dist,2] <- state_transitions[state_transitions$step == state_dist,2] + 1
      state_transitions[state_transitions$step == 0,2] <- state_transitions[state_transitions$step == 0,2] + 1
    }
    if (cell_state_labels[node] == "0"){
      N_violations <- N_violations + 2
    }
    if (length(children) > 1){
      for (i in 1:N_char){
        characters <- unique(barcodes_children[,i])
        if (length(characters) == 1){
          muts[node,i] <- characters
        }else{
          muts[node,i] <- "0"
          mutation_transitions[mutation_transitions$allele == characters[1],2] <- mutation_transitions[mutation_transitions$allele == characters[1],2] + 1
          mutation_transitions[mutation_transitions$allele == characters[2],2] <- mutation_transitions[mutation_transitions$allele == characters[2],2] + 1
        }
        total_muts <- total_muts + length(setdiff(characters, muts[node,i]))
      }
    }
  }
  state_transitions <- state_transitions[state_transitions$freq > 0,]
  mutation_transitions <- mutation_transitions[mutation_transitions$freq > 0,]
  l_expression <- sum(state_transitions$freq * log2(state_transitions$freq/sum(state_transitions$freq))) - N_violations*50
  l_barcode <- sum(mutation_transitions$freq * log2(mutation_transitions$freq/sum(mutation_transitions$freq))) - Q * total_muts - P * total_dist
  
  cl <- l_barcode + lambda * l_expression
  return(list(cl = cl,l_barcode = l_barcode, l_expression = l_expression))
}
