Preprocessing <- function(muts_leaves, phyla, states_leaves) {
  state_counter <- 1
  states_intermediate <- states_leaves[, 'cluster']
  for (cluster_id in sort(unique(states_leaves[, 'cluster']))) {
    depth_range <- sort(unique(states_leaves[states_leaves[, 'cluster'] == cluster_id, 'depth']))
    state_inter <- state_counter:(state_counter + length(depth_range) - 1)
    state_counter <- state_counter + length(depth_range)
    for (i in 1:length(depth_range)) {
      states_intermediate[states_leaves[, 'cluster'] == cluster_id & states_leaves[, 'depth'] == depth_range[i]] <- state_inter[i]
    }
  }

  state_lineages <- list()
  for (j in 1:6) {
    lineage <- c(sort(unique(states_intermediate[states_leaves[, 'cluster'] == 8])), sort(unique(states_intermediate[states_leaves[, 'cluster'] == j])))
    state_lineages[[length(state_lineages) + 1]] <- lineage
  }

  cell_meta <- as.data.frame(states_leaves)
  cell_meta$cluster_kmeans <- states_intermediate
  #muts_leaves <- MutsTransform(muts_leaves)

  return(list(muts_leaves, state_lineages, cell_meta))
}


main <- function() {
  #.libPaths("/nethome/xpan78/LinRace/Rlib")
  .libPaths("/nethome/pputta7/LinRace/Rlib")
  library("devtools")
  load_all()
  #load_all("/project/xpan78/LinRace-temp/")


  args <- commandArgs(trailingOnly = TRUE)
  ncells <- args[1]
  mu <- args[2]
  pd <- args[3]
  run <- args[4]

  print(args)


  cm_dir <- sprintf('sim/mut_%s_mu_%s_pd_%s_run_%s.csv', ncells, mu, pd, run)
  meta_dir <- sprintf('sim/meta_%s_mu_%s_pd_%s_run_%s.csv', ncells, mu, pd, run)
  count_dir <- sprintf('sim/expr_%s_mu_%s_pd_%s_run_%s.csv', ncells, mu, pd, run)
  out_dir <- sprintf('results/linrace_%s_mu_%s_pd_%s.run', ncells, mu, pd)
  ncells <- strtoi(ncells)

  cm <- read.csv(file = cm_dir, row.names = 1, stringsAsFactors = FALSE)
  cell_meta <- read.csv(file = meta_dir, row.names = 1, stringsAsFactors = FALSE)
  counts <- read.csv(file = count_dir, stringsAsFactors = FALSE)
  counts <- t(counts)
  rownames(counts) <- cell_meta[, "cellID"]

  cell_state_labels <- cell_meta[, "cluster"]

  phyla <- ape::read.tree(text = '((t1:4, t2:4, t3:4, t4:4, t5:4, t6:4):2);')

  returnList <- Preprocessing(cm, phyla, cell_meta)
  muts_leaves <- returnList[[1]]
  state_lineages <- returnList[[2]]
  cell_meta <- returnList[[3]]



  #order <- sample(nrow(muts_leaves))
  #cell_meta <- cell_meta[order,]
  #muts_leaves <- muts_leaves[order,]
  #counts <- counts[order,]

  #dist_h <- distH(muts_leaves);
  #print(as.matrix(dist_h)[1:5,])
  #tree <- NJ_asym(dist_h,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  scores <- c()
  times <- c()
  avg_score <- 0
  avg_time <- 0

  num_trials <- 5

  for (i in 1:num_trials) {
    start_time <- as.numeric(Sys.time())

    tree <- NJ_asym_Dif(muts_leaves, cell_meta$cluster_kmeans, counts, state_lineages, max_Iter = 500, newick = TRUE)
    tree_gt <- stree(ncells, type = "balanced")
    tip_label <- c()
    for (node in tree_gt$tip.label) {
      tip <- paste('cell_', substr(node, 2, nchar(node)), sep = '')
      tip_label <- c(tip_label, tip)
    }
    tree_gt$tip.label <- tip_label
    avg_score <- avg_score + RF.dist(tree, tree_gt, normalize = TRUE)
    end_time <- as.numeric(Sys.time())
    avg_time <- avg_time + (end_time - start_time)
  }

  avg_score <- avg_score / num_trials
  avg_time <- avg_time / num_trials

  write(sprintf("%f\t%f", avg_score, avg_time), file=out_dir, append=TRUE)
}

main()