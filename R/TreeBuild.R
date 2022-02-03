rowcol<-function(ix,n) { #given index, return row and column
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  cbind(nr,nc)
}

Diver_cal <- function(X,N){
  r <- rep(0,N)
  for (i in 1:length(X)){
    rc <- rowcol(i,N)
    r[rc] <- r[rc] + X[i]
  }
  return(r)
}

compute_Q <- function(X,r,N){
  Q <- X
  for (i in 1:length(X)){
    coord <- rowcol(i,N)
    Q[i] <- X[i] - sum(r[coord])/(N-2)
  }
  return(Q)
}

update_dist <- function(X,index1,index2,labels,new_node_id){
  X_m <- as.matrix(X)
  #X_new <- X_new[-c(index1,index2),]
  #X_new <- X_new[,-c(index1,index2)]
  N <- length(labels)

  labels_new = c(new_node_id,labels)
  results <- matrix(0, N+1, N+1)
  results[2:(N+1),2:(N+1)] <- X_m
  for (i in 1:N){
    dist_temp <- max(0.5*(X_m[i,index1] + X_m[i,index2] - X_m[index1,index2]),0)
    results[1,i+1] <- dist_temp
    results[i+1,1] <- dist_temp
  }
  rownames(results) <- labels_new
  colnames(results) <- labels_new
  
  results <- results[-c(index1+1,index2+1),]
  results <- results[,-c(index1+1,index2+1)]

  return(as.dist(results))
}


NJ <- function(X, muts = NA, states = NA, lineages = NA, it = 10000){
  if (is.matrix(X)) X <- as.dist(X)
  if (anyNA(X))
    stop("missing values are not allowed in the distance matrix\nConsider using njs()")
  if (any(is.infinite(X)))
    stop("infinite values are not allowed in the distance matrix")
  N <- as.integer(attr(X, "Size"))
  if (N < 3) stop("cannot build an NJ tree with less than 3 observations")
  labels <- attr(X, "Labels")
  if (is.null(labels)) labels <- as.character(1:N)

  node_definition <- NULL
  X_Q_min <- 1
  
  .merge <- function(merge_member_1, merge_member_2, X_val) {
    #compute distances between the merged node and the new node
    dist_merged_1 <- X_val / 2 - (r[merge_member_1] - r[merge_member_2])/(2*(N-2))
    if (dist_merged_1 < 0){
      dist_merged_1 <- 0
    } else if (dist_merged_1 > X_val){
      dist_merged_1 <- X_val
    }
    dist_merged_2 <- X_val - dist_merged_1

    node_definition <- sprintf('(%s:%f, %s:%f)',
                               labels[merge_member_1], dist_merged_1,labels[merge_member_2], dist_merged_2)

    X_Q_min <<- X_val
    X <<- update_dist(X,merge_member_1,merge_member_2,labels,node_definition)
    if (anyNA(X)) browser()

    N <<- as.integer(attr(X, "Size"))
    labels <<- attr(X, "Labels")
    
    node_definition
  }
  
  .find_node_idx <- function(id, tree) {
    if (id <= length(tree$tip.label)) {
      name <- tree$tip.label[id]
      lbl <- substr(name, 6, nchar(name))
    } else {
      lbl <- internal_nodemap[which(internal_nodemap[,1] == id), 2]
    }
    which(labels == lbl)
  }

  while(N > 3){
    r <- Diver_cal(X,N)
    Q <- compute_Q(X,r,N)
    Q_min_all <- which(Q == min(Q))
    if (!is.na(lineages) && length(Q_min_all) > 1) {
      candidates <- rowcol(Q_min_all, N)
      cand_labels <- labels[unique(as.vector(candidates))]
      cand_labels <- cand_labels[!is.na(strtoi(cand_labels))]
      if (length(cand_labels) >= 3) {
        # hill climb
        tip_labels <- paste0('cell_', cand_labels)
        res <- hill_climb(tip_labels, muts, states, lineages, it)
        
        internal_nodemap <- matrix(0, nrow = 0, ncol = 2)
        porder <- reorder(res, "postorder")$edge
        for (i in 1:(nrow(porder)/ 2)) {
          c1 <- porder[i * 2 - 1, 2]
          c1_idx <- .find_node_idx(c1, res)
          c2 <- porder[i * 2, 2]
          c2_idx <- .find_node_idx(c2, res)
          pr <- porder[i * 2, 1]
          
          ndef <- .merge(c1_idx, c2_idx, as.matrix(X)[c1_idx, c2_idx])
          internal_nodemap <- rbind(internal_nodemap, c(pr, ndef))
          if (N <= 3) break
        }
      } else {
        Q_min <- which.min(Q)[1]
        coord_merge <- rowcol(Q_min,N)
        .merge(coord_merge[1], coord_merge[2], X[Q_min])
      }
    } else {
      Q_min <- which.min(Q)[1]
      coord_merge <- rowcol(Q_min,N)
      .merge(coord_merge[1], coord_merge[2], X[Q_min])
    }
  }
  merge_member_1 = labels[1]
  merge_member_2 = labels[2]
  r <- Diver_cal(X,N)

  dist_merged_1 <- X[1]/2 - (r[1] - r[2])/(2*(N-2))
  if (dist_merged_1 < 0){
    dist_merged_1 <- 0
  } else if (dist_merged_1 > X_Q_min){
    dist_merged_1 <- X_Q_min
  }
  dist_merged_2 <- X_Q_min - dist_merged_1

  # ...then determine their distance to the other remaining node
  node_definition = labels[3]
  internal_len = max(0.5*(X[2]+X[3]-X[1]),0)
  # ...and finally create the newick string describing the whole tree.
  newick = sprintf("(%s:%f, %s:%f, %s:%f);",merge_member_1, dist_merged_1,
                                       node_definition, internal_len,
                                       merge_member_2, dist_merged_2)

  # return the phylo object transformed from newick.
  tree_nj <- read.tree(text = newick)
  return(tree_nj)
}

NJ_asym <- function(X,states){
  if (is.matrix(X)) X <- as.dist(X)
  if (anyNA(X))
    stop("missing values are not allowed in the distance matrix\nConsider using njs()")
  if (any(is.infinite(X)))
    stop("infinite values are not allowed in the distance matrix")
  N <- as.integer(attr(X, "Size"))
  if (N < 3) stop("cannot build an NJ tree with less than 3 observations")
  labels <- attr(X, "Labels")
  if (is.null(labels)) labels <- as.character(1:N)

  node_definition <- NULL

  while(N > 3){
    r <- Diver_cal(X,N)
    Q <- compute_Q(X,r,N)
    Q_min <- which.min(Q)
  if (length(Q_min) > 1){
    coord_merge <- c()
    for (entry in Q_min){
    coord_merge <- c(coord_merge,rowcol(entry,N))
    }
    states_merge <- states[coord_merge]
    node_definition <- FindExpTree(states_merge,labels = labels[coord_merge])
  }

    coord_merge <- rowcol(Q_min,N)
    merge_member_1 = coord_merge[1]
    merge_member_2 = coord_merge[2]

    #compute distances between the merged node and the new node
    dist_merged_1 <- X[Q_min]/2 - (r[merge_member_1] - r[merge_member_2])/(2*(N-2))
    if (dist_merged_1 < 0){
      dist_merged_1 <- 0
    } else if (dist_merged_1 > X[Q_min]){
      dist_merged_1 <- X[Q_min]
    }
    dist_merged_2 <- X[Q_min] - dist_merged_1

    node_definition <- sprintf('(%s:%f, %s:%f)',
                               labels[merge_member_1], dist_merged_1,labels[merge_member_2], dist_merged_2)

    X <- update_dist(X,merge_member_1,merge_member_2,labels,node_definition)

    N <- as.integer(attr(X, "Size"))
    labels <- attr(X, "Labels")
  }
  merge_member_1 = labels[1]
  merge_member_2 = labels[2]
  r <- Diver_cal(X,N)

  dist_merged_1 <- X[1]/2 - (r[1] - r[2])/(2*(N-2))
  if (dist_merged_1 < 0){
    dist_merged_1 <- 0
  } else if (dist_merged_1 > X[Q_min]){
    dist_merged_1 <- X[Q_min]
  }
  dist_merged_2 <- X[Q_min] - dist_merged_1

  # ...then determine their distance to the other remaining node
  node_definition = labels[3]
  internal_len = max(0.5*(X[2]+X[3]-X[1]),0)
  # ...and finally create the newick string describing the whole tree.
  newick = sprintf("(%s:%f, %s:%f, %s:%f);",merge_member_1, dist_merged_1,
                                       node_definition, internal_len,
                                       merge_member_2, dist_merged_2)

  # return the phylo object transformed from newick.
  tree_nj <- read.tree(text = newick)
  return(tree_nj)
}

FindExpTree <- function(states,labels,maxIter){
  tree <- rtree(length(states), rooted = TRUE, tip.label = labels)
  tree$edge.length <- rep(1, length(tree$edge.length))
  edges <- tree$edge
  #maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  max_likelihood <- LikelihoodCal_Exp(tree,states,state_lineages)

  likelihood_curve <-c()
  seeds <- runif(10000,1,99999)
  best_tree <- tree
  best_tree_list <- list()
  tree_list <- list()
  maxl_list <- c()
  ptm <- proc.time()
  LikelihoodCal_time <- 0
  TreeLC_time <- 0
  for (i in 1:maxIter){
    c_time <- proc.time()

    tree_new <- TreeLC2(tree)
    #tree_new <-  rNNI(tree,1,1)
    TreeLC_time <- TreeLC_time + proc.time()-c_time
    #cl <- LikelihoodCal(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
    c_time <- proc.time()
    likelihood_new <- LikelihoodCal_Exp(tree_new,states,state_lineages)
    #cl_list <- LikelihoodCal2(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
    LikelihoodCal_time <- LikelihoodCal_time + proc.time()-c_time
    likelihood_curve <- c(likelihood_curve,maxcl)
    if (likelihood_new > max_likelihood){
    max_likelihood <- likelihood_new
    tree <- tree_new
    best_tree <- tree
    tree_list[[length(tree_list)+1]] <- tree
    }
    if (i %% 10 == 0){
    likelihood_check <- sprintf("After %g iterations, the best likelihood is %g.", 
                  i, maxcl)
    print(likelihood_check)
    }
    if (i%%100 == 0){
    if (length(unique(likelihood_curve[(i-100):i])) == 1){
      #return the local optimal tree
    }
    }
  }
  #return the best tree found so far
  return()
}



DivideTree <- function(tree, depth_split){
  tree <- Preorder(tree)
  tree$edge.length <- rep(1,length(tree$edge.length))
  #bound_l <- floor(log2(N_split))
  #bound_h <- ceiling(log2(N_split))

  node_depth <- get_all_distances_to_root(tree)
  leaf_depth_min <- min(node_depth[1:length(tree$tip.label)])
  if (depth_split >= leaf_depth_min){
    stop("The dividing depth is below the depth of certain leaves.")
  }
  #depth_df <- data.frame(depth = numeric(),N_node = numeric())
  #for (depth in sort(unique(node_depth))){
  #  node_depth_list <- which(node_depth == depth)
  #  N_node <- length(node_depth_list)
  #  temp <- data.frame(depth = depth, N_node = N_node)
  #  depth_df <- rbind(depth_df,temp)
  #}
  node_split_list <- which(node_depth == depth_split)

  subtree_list <- list()
  tree_backbone  <- tree
  N_node_dropped <- 0
  for (node_split in sort(node_split_list,decreasing=TRUE)){
    subtree <- Subtree(tree,node_split)
    subtree$root.edge <- 1
    subtree$name <- paste("subtree",toString(node_split))
    subtree$edge.length <- rep(1, nrow(subtree$edge))
    subtree_list[[length(subtree_list)+1]] <- subtree
    tree_backbone <- bind.tip(tree_backbone, paste("subtree",toString(node_split)), where=node_split-N_node_dropped,edge.length = 1)
    tree_backbone <- drop.tip(tree_backbone,tip = subtree$tip.label, root.edge = 1)
    N_node_dropped <- N_node_dropped + length(subtree$tip.label)-1
  }
  return(list(tree_backbone,subtree_list))
}

ConstructTree <- function(tree_backbone,subtrees){
  tree <- tree_backbone
  for (subtree in subtrees){
    bind_tip <- subtree$name
    tree <- bind.tree(tree, subtree, where = which(tree_backbone$tip.label==bind_tip))
  }
  return(tree)
}

