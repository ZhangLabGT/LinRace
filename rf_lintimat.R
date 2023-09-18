main <- function() {

  library("devtools")
  load_all()

  args <- commandArgs(trailingOnly = TRUE)
  ncells <- args[1]
  run <- args[2]
  time <- args[3]
  mu <- 0.1
  pd <- 1

  tree_dir <- sprintf("results/lintimat_%s_%s_bin_tree.newick", ncells, run)

  ncells <- strtoi(ncells)

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

  out_dir <- sprintf('results/lintimat_%s_mu_%s_pd_%s.run', ncells, mu, pd)

  write(sprintf("%f\t%f", score, time), file = out_dir, append = TRUE)

}

main()