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

Generate_profile_multi <- function(observed_counts, muts, states){
  ncells <- nrow(states)
  allele <- c()
  m_l <- c('D','I')
  for (k in 1:40){
    mutator <- m_l[sample(c(1,2),1)]
    num1 <- sample(1:100,1)
    num2 <- sample((2*num1):300,1)
    allele_temp <- sprintf("%g%s+%g",
                           num1, mutator,num2)
    allele <- c(allele,allele_temp)
  }

  #counts <- observed_counts[[1]]
  counts <- observed_counts
  exdata <- preprocessSC(counts)
  sf <- apply(exdata, 2, mean)
  npX <- t(t(exdata) / sf )
  lnpX <- log(npX+1)
  lnpX_imp <- DrImpute(lnpX)

  counts_t <- t(lnpX_imp)
  colnames(counts_t) <- paste("gene", seq(1, ncol(counts_t)), sep = "_")
  profiles <- data.frame(cell_id=character(), ClusterIdent = numeric(), HMID = character(),
                         stringsAsFactors = F)
  for (k in 1:ncells){
    barcode <- muts[k,]
    barcode <- replace(barcode, barcode=='-', '100D+100')
    barcode <- replace(barcode, barcode==0, 'None')
    for (j in 1:20){
      barcode <- replace(barcode, barcode==j, allele[j])
    }
    barcode <- paste(barcode,collapse = "-")
    cell_id <- paste("cell",states[k, 4],sep = "_")
    #cell_id <- paste("f6_DEW",states_leaves[i,4],"bcAAVQ",sep = "_")
    temp_df <- data.frame(cell_id=cell_id,ClusterIdent = states[k,2], HMID = barcode,
                          stringsAsFactors = F)
    profiles <- rbind(profiles,temp_df)
  }

  profiles_out <- cbind(profiles,counts_t)

  var_gene <- c()
  for (gene in colnames(counts_t)){
    var_gene <- c(var_gene,var(counts_t[,gene]))
  }
  rank_gene <- data.frame(name = colnames(counts_t), var = var_gene)
  rank_gene <- rank_gene[order(-rank_gene$var),]
  colnames(rank_gene)<- c()
  returnList <- list(profiles_out,rank_gene)

  return(returnList)
}

main <- function() {
  library("devtools")
  library(DrImpute)
  load_all()

  cm_dir <- 'sim/mut_2048_mu_0.1_pd_0_run_1.csv'
  meta_dir <- 'sim/meta_2048_mu_0.1_pd_0_run_1.csv'
  count_dir <- 'sim/expr_2048_mu_0.1_pd_0_run_1.csv'

  combined_profile_dir <- 'sim/profile_2048_mu_0.1_pd_0_run_1.txt'
  top_genes_dir <- 'sim/top_genes_2048_mu_0.1_pd_0_run_1.txt'

  cm <- read.csv(file = cm_dir, row.names = 1, stringsAsFactors = FALSE)
  cell_meta <- read.csv(file = meta_dir, row.names = 1, stringsAsFactors = FALSE)
  counts <- read.csv(file = count_dir, stringsAsFactors = FALSE)
  counts <- t(counts)
  rownames(counts) <- cell_meta[, "cellID"]

  cell_state_labels <- cell_meta[, "cluster"]
  returnList <- Preprocessing(cm, phyla, cell_meta)
  muts_leaves <- returnList[[1]]
  state_lineages <- returnList[[2]]
  cell_meta <- returnList[[3]]
  counts <- t(counts)

  profile_res <- Generate_profile_multi(counts,muts_leaves,cell_meta)
  profile_out <- profile_res[[1]]
  top_genes <- profile_res[[2]]
  order <- sample(nrow(profile_out))
  profile_out <- profile_out[order,]

  write.table(profile_out,file = combined_profile_dir,row.names = FALSE,sep = "\t", quote = FALSE)
  write.table(subset(top_genes,select = 1),top_genes_dir,row.names = FALSE,sep = "\t", quote = FALSE)
}

main()