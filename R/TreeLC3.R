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
