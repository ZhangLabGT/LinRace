AncesInfer_sub <- function(tree,muts,cell_state_labels,state_lineages,newick_lookup){
  N_char <- dim(muts)[2]
  #Get the internal node ids
  nodes_internal <- unique(tree$edge[,1])
  #order the internal nodes
  nodes_internal <- sort(nodes_internal, decreasing = TRUE)

  #extend mutation table to include the internal nodes
  #muts <- rbind(muts,matrix(0,length(nodes_internal),N_char))
  #cell_state_labels <- c(cell_state_labels,rep(0,length(nodes_internal)))
  internal_lookup <- data.frame(node = numeric(),index = numeric())

  #determine ancestor barcode based on child barcodes
  for (node in nodes_internal){
    children_index <- tree$edge[tree$edge[,1] == node,2]
    children <- c()
    for (cell_id in children_index){
      if (cell_id > length(tree$tip.label)){
        cell_row <- internal_lookup$index[internal_lookup$node == node]
      }else {
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
      break
    }
    cell_state_labels <- c(cell_state_labels,state_parent)
    temp <- data.frame(node = node,index = length(cell_state_labels))
    internal_lookup <- rbind(internal_lookup,temp)

    if (length(children) > 1){
      barcode_par <- barcodes_children[1,]
      for (i in 1:N_char){
        characters <- unique(barcodes_children[,i])
        if (length(characters) == 1){
          barcode_par[i] <- characters
        }else{
          barcode_par[i] <- "0"
        }
      }
      muts <- rbind(muts,barcode_par)
    }
  }
  return(list(muts,cell_state_labels))
}
