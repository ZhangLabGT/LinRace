FindSubtreeNewick <- function(tree, node, num_leaves) {
  # assumes node is listed through Preorder indexing
  tree_strsplit <- strsplit(tree, "")[[1]]
  if (node <= num_leaves) {
    i <- 1
    while (node > 0 && i < length(tree_strsplit)) {
      if (tree_strsplit[i] != "(" &&
          tree_strsplit[i] != ")" &&
          tree_strsplit[i] != ",") {
        node <- node - 1

        j <- i

        while (tree_strsplit[j] != "(" &&
               tree_strsplit[j] != ")" &&
               tree_strsplit[j] != ",") {
          j <- j + 1
        }
        if (node == 0) {
          return(c(i, j - 1))
        }
        i <- j
      }
      i <- i + 1
    }
  }
  node <- node - num_leaves
  start <- 1
  for (i in seq_along(tree_strsplit)) {
    if (tree_strsplit[i] == "(") {
      node <- node - 1
    }
    if (node == 0) {
      # traversed down enough of tree
      start <- i
      break
    }
  }

  # find end of subtree in newick string
  num_open <- 0
  i <- start
  while ((num_open > 0 || i == start) && i < length(tree_strsplit)) {
    ch <- tree_strsplit[i]
    if (ch == "(") {
      num_open <- num_open + 1
    } else if (ch == ")") {
      num_open <- num_open - 1
    }

    if (num_open == 0) {
      break
    }

    i <- i + 1
  }
  end <- i
  return(c(start, end))
}

# checks if two ranges overlap
OverlapNewick <- function(r1, r2) {
  if (r1[1] <= r2[1] && r1[2] >= r2[1]) {
    return(TRUE)
  } else if (r1[1] >= r2[1] && r1[1] <= r2[2]) {
    return(TRUE)
  } else if (r1[2] >= r2[1] && r1[2] <= r2[2]) {
    return(TRUE)
  } else if (r2[1] <= r1[1] && r2[2] >= r1[1]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

SwapSubtrees <- function(tree, r1, r2, num_leaves) {
  # assumes subtrees provided are disjoint

  subtree1_loc <- FindSubtreeNewick(tree, r1, num_leaves)
  subtree2_loc <- FindSubtreeNewick(tree, r2, num_leaves)

  subtree1 <- substr(tree, subtree1_loc[1], subtree1_loc[2])
  subtree2 <- substr(tree, subtree2_loc[1], subtree2_loc[2])

  # replace subtrees with placeholders in newick tree string
  tree <- stringi::stri_replace_first_fixed(tree, subtree1, "_TMP_")
  tree <- stringi::stri_replace_first_fixed(tree, subtree2, subtree1)
  tree <- stringi::stri_replace_first_fixed(tree, "_TMP_", subtree2)

  return(tree)
}

SwapTreeTools <- function(tree, r1, r2) {

  tree <- Preorder(tree)
  tree$node.label <- paste0("node", seq(1:tree$Nnode))

  subtree_swap <- c(r1, r2)
  subtree1 <- Subtree(tree, subtree_swap[1])
  subtree2 <- Subtree(tree, subtree_swap[2])

  while (length(intersect(subtree1$tip.label, subtree2$tip.label)) > 0) {
    print("intersection ocurred!")
    subtree_swap <- sample(tree$Nnode + length(tree$tip.label), 2)
    subtree1 <- Subtree(tree, subtree_swap[1])
    subtree2 <- Subtree(tree, subtree_swap[2])
  }

  subtree1$name <- subtree_swap[1]
  subtree1$edge.length <- rep(1, nrow(subtree1$edge))
  subtree2$name <- subtree_swap[2]
  subtree2$edge.length <- rep(1, nrow(subtree2$edge))

  root1 <- tree$edge[tree$edge[, 2] == subtree1$name, 1]
  root2 <- tree$edge[tree$edge[, 2] == subtree2$name, 1]

  if (root1 > root2) {
    reg <- subtree2
    subtree2 <- subtree1
    subtree1 <- reg
  }
  root1 <- tree$edge[tree$edge[, 2] == subtree1$name, 1]
  root2 <- tree$edge[tree$edge[, 2] == subtree2$name, 1]
  root1 <- paste0("node", root1 - Ntip(tree))
  root2 <- paste0("node", root2 - Ntip(tree))
  subtree1$root.edge <- 1
  subtree2$root.edge <- 1
  Ntip1 <- subtree1$Ntip
  Ntip2 <- subtree2$Ntip

  if (length(subtree1$tip.label) > 1) {
    tree_prune1 <- drop.tip(tree, subtree1$tip.label, trim.internal = FALSE, subtree = FALSE, root.edge = 1)
  } else {
    tree_prune1 <- tree
    tree_prune1$tip.label[tree_prune1$tip.label == subtree1$tip.label] <- "NA"
  }
  if (length(subtree2$tip.label) > 1) {
    tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE, root.edge = 1)
  } else {
    tree_prune2 <- tree_prune1
    tree_prune2$tip.label[tree_prune2$tip.label == subtree2$tip.label] <- "NA"
  }
  #tree_new1 <- bind.tree(tree_prune1, subtree1, where = subtree2$node.label[1]-Nnode1)
  #tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  #print(regraft_loc)
  #N_na <- sum(tree_prune2$tip.label == 'NA')
  tree_new1 <- bind.tree(tree_prune2, subtree2, where = Ntip(tree_prune2) + match(root1, tree_prune2$node.label))
  tree_new2 <- bind.tree(tree_new1, subtree1, where = Ntip(tree_new1) + match(root2, tree_new1$node.label))
  tree_final <- drop.tip(tree_new2, tree_new2$tip.label[grepl("node", tree_new2$tip.label, fixed = TRUE)], trim.internal = TRUE)
  tree_final <- drop.tip(tree_final, "NA")
  tree_final$node.label <- NULL
  tree_final$edge.length <- rep(1, length(tree_final$edge.length))
  #plot(tree_final)
  return(tree_final)
}

TestLoopNewick <- function(tree, r1, r2) {
  tree_str <- stringi::stri_replace_all_fixed(TreeTools::NewickTree(tree), ":1", "")

  for (i in 1:1000) {
    SwapSubtrees(tree_str, r1, r2, length(tree[['tip.label']]))
  }
}

TestLoopTreeTools <- function(tree, r1, r2) {
  #tree_str <- stringi::stri_replace_all_fixed(TreeTools::NewickTree(tree), ":1", "")

  for (i in 1:1000) {
    SwapTreeTools(tree, r1, r2)
  }
}

Check <- function(tree, r1, r2, should_plot = FALSE) {
  tree_str <- stringi::stri_replace_all_fixed(TreeTools::NewickTree(tree), ":1", "")
  swap1 <- SwapSubtrees(tree_str, r1, r2, length(tree[['tip.label']]))
  swap2 <- stringi::stri_replace_all_fixed(TreeTools::NewickTree(SwapTreeTools(tree, r1, r2)), ":1", "")

  if (swap1 == swap2) {
    print("Passed Check.")
  } else {
    print("Failed Check.")
  }
  if (should_plot) {

    plot(ape::read.tree(text = swap1), main = "Newick Swapping")
    plot(ape::read.tree(text = swap2), main = "TreeTools Swapping")
  }
}

TestSwaps <- function() {
  #file_name <- "sim/tree_2048_mu_0.05_pd_1_run_10.tree"
  file_name <- "sim/tree_256_mu_0.1_pd_0_run_1.tree"

  tree <- ape::read.tree(file_name)
  tree <- TreeTools::Preorder(tree)

  SwapSubtrees(NewickTree(Subtree(tree, 262)), 1, 2, 8)

  #r1 <- 2060
  #r2 <- 2059
  r1 <- 266
  r2 <- 264

  print(file_name)
  print(paste(replicate(25, "-"), collapse = ""))
  print("Profiling Newick Swap")
  system.time(TestLoopNewick(tree, r1, r2))

  print("Profiling TreeTools Swap")
  system.time(TestLoopTreeTools(tree, r1, r2))
  #check(tree, r1, r2, should_plot=FALSE)

}

TreeLC2Newick <- function(tree_newick, num_leaves, prob_internal = NULL) {
  num_nodes <- 2 * num_leaves - 1
  tree_newick <- stringi::stri_replace_all_fixed(tree_newick, ":1", "")

  if (is.null(prob_internal)) {
    prob_internal <- dunif(1:num_nodes, 1, num_nodes)
  }

  subtree_swap <- sample(num_nodes, 2, prob = prob_internal)
  subtree1 <- FindSubtreeNewick(tree_newick, subtree_swap[1], num_leaves)
  subtree2 <- FindSubtreeNewick(tree_newick, subtree_swap[2], num_leaves)

  while (OverlapNewick(subtree1, subtree2)) {
    subtree_swap <- sample(num_nodes, 2, prob = prob_internal)
    subtree1 <- FindSubtreeNewick(tree_newick, subtree_swap[1], num_leaves)
    subtree2 <- FindSubtreeNewick(tree_newick, subtree_swap[2], num_leaves)
  }

  tree_new <- SwapSubtrees(tree_newick, subtree_swap[1], subtree_swap[2], num_leaves)
  return(tree_new)
}
