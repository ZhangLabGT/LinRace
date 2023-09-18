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

distH <- function(muts_leaves) {
  sim_tree <- list()
  states_barcode <- unique(unlist(muts_leaves))
  states_mutated <- states_barcode[states_barcode != "0"]
  nstates <- length(states_mutated)
  #overflow_states <- c()
  #if (nstates > 26){
  #  overflow_states <- states_mutated[27:nstates]
  #}
  for (i in 1:dim(muts_leaves)[1]) {
    vec <- as.character(muts_leaves[i,])
    #vec[vec %in% overflow_states] <- "0"
    sim_tree[[i]] <- vec
    #names(sim_tree[[i]]) <- rownames(muts_leaves)[i]
  }
  #print(length(sim_tree))
  names(sim_tree) <- rownames(muts_leaves)
  #nmstrings <- c('0','-',c(1:100))
  nmstrings <- unique(unlist(sim_tree))

  sim_data <- phyDat(sim_tree, type = 'USER', levels = nmstrings)
  sim_data <- subset(sim_data, 1:length(sim_tree))
  #print(length(sim_data))

  #transform to make k-mer work
  # new_levels <- c()
  # i <- 1
  # for (state in levels(sim_data)){
  #   if (!state %in% c("0","-")){
  #       new_levels <- c(new_levels,LETTERS[i])
  #       i <- i + 1
  #   }
  # }

  #levels(sim_data) <- new_levels

  dist_kmer <- dist_replacement(sim_data, reps = 20L, k = 2L)
  #dist_h <- dist.hamming(sim_data)
  #print(dim(as.matrix(dist_h)))
  return(dist_kmer)
}

main <- function() {
  #.libPaths("/nethome/pputta7/LinRace/Rlib")
  library("devtools")
  load_all()
  #load_all("/project/pputta7/LinRace-temp/")
  library("DCLEAR")


  args <- commandArgs(trailingOnly = TRUE)

  args <- commandArgs(trailingOnly = TRUE)
  ncells <- args[1]
  mu <- args[2]
  pd <- args[3]
  run <- args[4]

  print(args)


  cm_dir <- sprintf('sim/mut_%s_mu_%s_pd_%s_run_%s.csv', ncells, mu, pd, run)
  meta_dir <- sprintf('sim/meta_%s_mu_%s_pd_%s_run_%s.csv', ncells, mu, pd, run)
  count_dir <- sprintf('sim/expr_%s_mu_%s_pd_%s_run_%s.csv', ncells, mu, pd, run)
  out_dir <- sprintf('results/dclear_%s_mu_%s_pd_%s.run', ncells, mu, pd)
  ncells <- strtoi(ncells)

  cm <- read.csv(file = cm_dir, row.names = 1, stringsAsFactors = FALSE)
  cell_meta <- read.csv(file = meta_dir, row.names = 1, stringsAsFactors = FALSE)

  cell_state_labels <- cell_meta[, "cluster"]

  phyla <- read.tree(text = '((t1:4, t2:4, t3:4, t4:4, t5:4, t6:4):2);')

  returnList <- Preprocessing(cm, phyla, cell_meta)
  muts_leaves <- returnList[[1]]
  state_lineages <- returnList[[2]]
  cell_meta <- returnList[[3]]

  #order <- sample(nrow(muts_leaves))
  #cell_meta <- cell_meta[order,]
  #muts_leaves <- muts_leaves[order,]
  start_time <- as.numeric(Sys.time())

  dist_kmer <- distH(muts_leaves)
  #print(as.matrix(dist_h)[1:5,])
  #tree <- NJ_asym(dist_h,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  #tree <- NJ_asym_alter(muts_leaves,cell_meta$cluster_kmeans,state_lineages,max_Iter = 200)
  tree <- nj(dist_kmer)


  tree_gt <- stree(ncells, type = "balanced")
  tip_label <- c()
  for (node in tree_gt$tip.label) {
    tip <- paste('cell_', substr(node, 2, nchar(node)), sep = '')
    tip_label <- c(tip_label, tip)
  }
  tree_gt$tip.label <- tip_label
  score <- RF.dist(tree, tree_gt, normalize = TRUE)
  end_time <- as.numeric(Sys.time())
  time <- end_time - start_time

  write(sprintf("%f\t%f", score, time), file=out_dir, append=TRUE)

}

main()
