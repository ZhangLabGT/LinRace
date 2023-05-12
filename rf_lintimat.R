main <- function() {

  library("devtools")
  load_all()
  ncells <- 256

  tree_dir <- sprintf("results/lintimat_%d_bin_tree.newick", ncells)

  tree <- ape::read.tree(file = tree_dir)
  tree <- drop.tip(tree, "normal")

  tree_gt <- stree(ncells, type = "balanced")
  tip_label <- c()
  for (node in tree_gt$tip.label) {
    tip <- paste('cell_', substr(node, 2, nchar(node)), sep = '')
    tip_label <- c(tip_label, tip)
  }
  tree_gt$tip.label <- tip_label
  score <- RF.dist(tree, tree_gt, normalize = TRUE)
  print(score)
}

main()