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
  #.libPaths("/nethome/pputta7/LinRace/Rlib")
  library("devtools")
  #load_all("R")

  start_time <- as.numeric(Sys.time())
  print(start_time)
  args <- commandArgs(trailingOnly = TRUE)
  #cm_dir <- 'sim/mut_0.8_pa_1024_mu_0.35_pd_1_Nchar_16_run_2.csv'
  #meta_dir <- 'sim/meta_0.8_pa_1024_mu_0.35_pd_1_Nchar_16_run_2.csv'
  #count_dir <- 'sim/expr_0.8_pa_1024_mu_0.35_pd_1_Nchar_16_run_2.csv'

  ncells <- args[1]
  print(ncells)

  cm_dir <- sprintf('sim/mut_%d_mu_0.1_pd_0_run_1.csv', ncells)
  meta_dir <- sprintf('sim/meta_%d_mu_0.1_pd_0_run_1.csv', ncells)
  count_dir <- sprintf('sim/expr_%d_mu_0.1_pd_0_run_1.csv', ncells)
  out_dir <- sprintf('results/linrace_%d_bin_tree.newick', ncells)

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


  order <- sample(nrow(muts_leaves))
  cell_meta <- cell_meta[order,]
  muts_leaves <- muts_leaves[order,]
  counts <- counts[order,]
  #dist_h <- distH(muts_leaves)
  #print(as.matrix(dist_h)[1:5,])
  #tree <- NJ_asym(dist_h,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  tree <- NJ_asym_Dif(muts_leaves, cell_meta$cluster_kmeans, counts, state_lineages, max_Iter = 200)
  tree_gt <- stree(ncells,type = "balanced")
  score <- RF.dist(tree, tree_gt, normalize=TRUE)
  print(score)
  end_time <- as.numeric(Sys.time())
  print(end_time)
  print(sprintf('Total Time: %gs', (end_time - start_time)))
  write.tree(tree, file = out_dir)
}

main()
