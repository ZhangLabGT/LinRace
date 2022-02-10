distdex<-function(i,j,n){
  n*(i-1) - i*(i-1)/2 + j-i
} #given row, column, and n, return index

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
  X <- as.matrix(X)
  N <- length(labels)

  labels_new = c(new_node_id,labels)
  results <- matrix(0, N+1, N+1)
  results[2:(N+1),2:(N+1)] <- X

  for (i in 2:(N+1)){
    dist_temp <- 0.5*(X[i-1,index1]+X[i-1,index2]-X[index1,index2])
    results[1,i] <- dist_temp
    results[i,1] <- dist_temp
  }
  rownames(results) <- labels_new
  colnames(results) <- labels_new

  results <- tryCatch({
    results <- results[-c(index1+1,index2+1),]
    results <- results[,-c(index1+1,index2+1)]
    return(as.dist(results))
  }, error = function(e) {
    print(paste("Only one node left."))
    newick <- labels_new[1]
    return(newick)
  })
}

update_dist_multi <- function(X,index_list,labels,new_node_id){
  X <- as.matrix(X)
  #X_new <- X_new[-c(index1,index2),]
  #X_new <- X_new[,-c(index1,index2)]
  N <- length(labels)

  labels_new = c(new_node_id,labels)
  results <- matrix(0, N+1, N+1)
  results[2:(N+1),2:(N+1)] <- X
  innerdist_sum <- 0
  for (i in 1:(length(index_list)-1)){
    innerdist_sum <- innerdist_sum + X[index_list[i],index_list[i+1]]
  }
  for (i in 2:(N+1)){
    out_dist_sum <- 0
    for (index in index_list){
      out_dist_sum <- out_dist_sum + X[i-1,index]
    }

    dist_temp <- 0.5*(out_dist_sum-innerdist_sum)

    results[1,i] <- dist_temp
    results[i,1] <- dist_temp
  }
  rownames(results) <- labels_new
  colnames(results) <- labels_new


  #browser()
  results <- tryCatch({
    results <- results[-c(index_list+1),]
    results <- results[,-c(index_list+1)]
    return(as.dist(results))
  }, error = function(e) {
    print(paste("Only one node left."))
    newick <- labels_new[1]
    return(newick)
  })
}


NJ <- function(X){
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
    Q_min <- which.min(Q)[1]

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

#' LinRace main function: asymmetric division based Neighbor Joining
#'
#' @param X input distance matrix
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @import data.table
#' @export
NJ_asym <- function(X,muts,states,state_lineages){
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

  newick_lookup <- data.frame(newick = character(),index = numeric())
  while(N > 3){
    r <- Diver_cal(X,N)
    Q <- compute_Q(X,r,N)
    Q_min <- which(Q == min(Q))
	  coord_merge <- c()
	  for (entry in Q_min){
		  coord_merge <- c(coord_merge,rowcol(entry,N))
	  }
	  coord_merge <- unique(coord_merge)

	  muts_merge <- muts[coord_merge,]
	  DT <- as.data.table(muts_merge)
	  group_count <- lapply(1:nrow(DT), function(i){
	    DT_group_count <- apply(DT, 1, function(x) DT_compare <- sum(which(sum(x != DT[i,]) == 0)))
	    #DT_group_count <- sum(DT_group_count)
	    })
	  #sorting
	  for (i in 1:length(group_count)){
	    if (sum(group_count[[i]]) >1){
	      coord_merge <- coord_merge[which(group_count[[i]]==1)]
	      break
	    }
	  }
	  #browser()
	  if(length(coord_merge)>2){
	    states_merge <- states[coord_merge]
	    #consider barcode when all barcodes are different
	    tree_recon <- FindExpTree(states,labels = labels[coord_merge],state_lineages,muts,newick_lookup,maxIter = 200)
	    best_tree <- tree_recon[[1]]

	    node_definition <- write_tree(best_tree,quoting = 0)
	    if (length(coord_merge)==N){
	      print("All nodes are merged.")
	      tree_nj <- read.tree(text = node_definition)
	      return(tree_nj)
	    }


	    n <- nchar(node_definition)
	    if(substr(node_definition, n, n) == ';'){
	      node_definition <- substr(node_definition,1,n-1)
	    }
	    merged_state <- tree_recon[[2]]
	    merged_barcode <- tree_recon[[3]]
	    states <- c(states,merged_state)
	    muts <- rbind(muts,merged_barcode)

	    temp <- data.frame(newick = node_definition,index = length(states))
	    newick_lookup <- rbind(newick_lookup,temp)

	    #browser()
	    nodes_internal <- unique(best_tree$edge[,1])
	    #order the internal nodes
	    nodes_internal <- sort(nodes_internal, decreasing = TRUE)

	    node_newick_df <- data.frame(newick = character(),tree_index = numeric())
	    for(node in nodes_internal){
	      children_index <- best_tree$edge[best_tree$edge[,1] == node,2]
	      coord_merge <- c()
	      for (cell_id in children_index){
	        if (cell_id > length(best_tree$tip.label)){
	          node_name <- node_newick_df$newick[node_newick_df$tree_index==cell_id]
	          node_index <- which(labels==node_name)
	        }else {
	          node_name <- best_tree$tip.label[cell_id]
	          node_index <- which(labels==node_name)
	        }
	        coord_merge <- c(coord_merge, node_index)
	      }

	      merge_member_1 <- coord_merge[1]
	      merge_member_2 <- coord_merge[2]

	      index <- distdex(merge_member_1,merge_member_2,N)
	      if(N >2){
	        dist_merged_1 <- X[index]/2 - (r[merge_member_1] - r[merge_member_2])/(2*(N-2))
	      }else{
	        dist_merged_1 <- X[index]/2
	      }

	      if (dist_merged_1 < 0){
	        dist_merged_1 <- 0
	      } else if (dist_merged_1 > X[index]){
	        dist_merged_1 <- X[index]
	      }
	      dist_merged_2 <- X[index] - dist_merged_1

	      for (i in coord_merge){
	        label <- labels[i]
	        n <- nchar(label)
	        if(substr(label, n, n) == ';'){
	          labels[i] <- substr(label,1,n-1)
	        }
	      }

	      node_definition <- sprintf('(%s:%f, %s:%f)',
	                                 labels[merge_member_1], dist_merged_1,labels[merge_member_2], dist_merged_2)

	      X <- update_dist(X,coord_merge[1],coord_merge[2],labels,node_definition)
	      #N <- as.integer(attr(X, "Size"))
	      #labels <- attr(X, "Labels")

	      if (is.character(X)){
	        #browser()
	        print(X)
	        if(substr(X,nchar(X),nchar(X))!=";"){
	          X <- paste(X,";",sep = "")
	        }
	        tree_nj <- read.tree(text = X)
	        return(tree_nj)
	      }else{
	        N <- as.integer(attr(X, "Size"))
	        labels <- attr(X, "Labels")
	      }

	      temp <- data.frame(newick = node_definition,tree_index = node)
	      node_newick_df <- rbind(node_newick_df,temp)
      }

	    #X <- update_dist_multi(X,coord_merge,labels,node_definition)

	    # if (is.character(X)){
	    #   #browser()
	    #   print(X)
	    #   if(substr(X,nchar(X),nchar(X))!=";"){
	    #     X <- paste(X,";",sep = "")
	    #   }
	    #   tree_nj <- read.tree(text = X)
	    #   return(tree_nj)
	    # }else{
	    #   N <- as.integer(attr(X, "Size"))
	    #   labels <- attr(X, "Labels")
	    # }

	  }else{
  	  merge_member_1 = coord_merge[1]
  	  merge_member_2 = coord_merge[2]

  	  #compute distances between the merged node and the new node
  	  index <- distdex(merge_member_1,merge_member_2,N)
  	  dist_merged_1 <- X[index]/2 - (r[merge_member_1] - r[merge_member_2])/(2*(N-2))
  	  if (dist_merged_1 < 0){
  	    dist_merged_1 <- 0
  	  } else if (dist_merged_1 > X[index]){
  	    dist_merged_1 <- X[index]
  	  }
  	  dist_merged_2 <- X[index] - dist_merged_1

  	  for (i in coord_merge){
  	    label <- labels[i]
  	    n <- nchar(label)
  	    if(substr(label, n, n) == ';'){
  	      labels[i] <- substr(label,1,n-1)
  	    }
  	  }
      #browser()
  	  node_definition <- sprintf('(%s:%f, %s:%f)',
  	                             labels[merge_member_1], dist_merged_1,labels[merge_member_2], dist_merged_2)

  	  X <- update_dist(X,merge_member_1,merge_member_2,labels,node_definition)

  	  N <- as.integer(attr(X, "Size"))
  	  labels <- attr(X, "Labels")
  	}
  }

  for (i in 1:length(labels)){
    label <- labels[i]
    n <- nchar(label)
    if(substr(label, n, n) == ';'){
      labels[i] <- substr(label,1,n-1)
    }
  }

  if(length(labels)==2){
    # ...and finally create the newick string describing the whole tree.
    newick = sprintf("(%s:%f, %s:%f);",labels[1], 1,
                     labels[2], 1)
    #print(newick)
    #browser()
    # return the phylo object transformed from newick.
    tree_nj <- read.tree(text = newick)
  }else{
    newick = sprintf("(%s:%f, %s:%f, %s:%f);",labels[1], 1,
                     labels[2], 1, labels[3], 1)
    #print(newick)
    #browser()
    # return the phylo object transformed from newick.
    tree_nj <- read.tree(text = newick)
  }
  return(tree_nj)
}


#' LinRace main function: asymmetric division based Neighbor Joining
#'
#' @param X input distance matrix
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @import castor
FindExpTree <- function(states,labels,state_lineages,muts,newick_lookup,maxIter){
  #browser()
  for (i in 1:length(labels)){
    label <- labels[i]
    n <- nchar(label)
    if(substr(label, n, n) == ';'){
      labels[i] <- substr(label,1,n-1)
    }
  }

	tree <- rtree(length(labels), rooted = TRUE, tip.label = labels)
	tree$edge.length <- rep(1, length(tree$edge.length))
	edges <- tree$edge
	#maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
	max_likelihood <- LikelihoodCal_Exp(tree,states,state_lineages,newick_lookup)

	likelihood_curve <-c()
	seeds <- runif(10000,1,99999)
	best_tree <- tree
	best_tree_list <- list()
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
	  likelihood_new <- LikelihoodCal_Exp(tree_new,states,state_lineages,newick_lookup)
	  #cl_list <- LikelihoodCal2(tree_new,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
	  LikelihoodCal_time <- LikelihoodCal_time + proc.time()-c_time
	  likelihood_curve <- c(likelihood_curve,max_likelihood)
	  if (likelihood_new > max_likelihood){
  		max_likelihood <- likelihood_new
  		tree <- tree_new
  		best_tree <- tree
	  }
	  if (i %% 100 == 0){
  		likelihood_check <- sprintf("After %g iterations, the best likelihood is %g.",
  									i, max_likelihood)
  		print(likelihood_check)
	  }
	  if (i%%20 == 0){
  		if (length(unique(likelihood_curve[(i-20):i])) == 1){
  		  #record the local optima, and
  		  best_tree_list[[length(best_tree_list)+1]] <- best_tree
  		  maxl_list <- c(maxl_list,max_likelihood)

  		  tree <- rtree(length(labels), rooted = TRUE, tip.label = labels)
  		  tree$edge.length <- rep(1, length(tree$edge.length))
  		  #maxcl <- LikelihoodCal(tree,muts_leaves,cell_meta$cluster_kmeans,state_lineages)
  		  max_likelihood <- LikelihoodCal_Exp(tree,states,state_lineages,newick_lookup)
  		}
	  }
	}
	#browser()
	#return the best tree found so far
	best_tree <- best_tree_list[[which.max(maxl_list)]]
	ances_res <- AncesInfer_sub(tree,muts,states,state_lineages,newick_lookup)
	root_barcode <- ances_res[[1]][dim(ances_res[[1]])[1],]
	root_state_label <- ances_res[[2]][length(ances_res[[2]])]

	return(list(best_tree,root_state_label,root_barcode))
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

