FindPairs <- function(tree_backbone,dt){
  edges_backbone <- as.data.frame(tree_backbone$edge)
  colnames(edges_backbone) <- c("ances","desc")
  edges_backbone_groups <- edges_backbone %>% dplyr::filter(desc < length(tree_backbone$tip.label)) %>% dplyr::group_split(ances)
  edges_backbone_groups <- edges_backbone_groups[unlist(lapply(edges_backbone_groups,function(x){nrow(x)})) == 2]

  df <- data.frame(parent = numeric(),left = numeric(),right = numeric(),desc_side = character(),asc_side = character())
  for(branch in edges_backbone_groups){
    branch %>% arrange(desc)
    barcode_left <- dt$barcode[branch$desc[1]]
    barcode_right <- dt$barcode[branch$desc[2]]

    left_bucket <- stringr::str_split(barcode_left,"\\|",simplify = TRUE)
    right_bucket <- stringr::str_split(barcode_right,"\\|",simplify = TRUE)

    #replace dropout by 99999
    left_bucket[left_bucket == "-"] <- "99999"
    right_bucket[right_bucket == "-"] <- "99999"
    left_bucket <- as.numeric(left_bucket)
    right_bucket <- as.numeric(right_bucket)


    left_diff_bucket <- left_bucket[left_bucket != right_bucket]
    right_diff_bucket <- right_bucket[left_bucket != right_bucket]
    left_diff_bucket[left_diff_bucket != 0] <- 1
    right_diff_bucket[right_diff_bucket != 0] <- 1
    if (all(left_diff_bucket > right_diff_bucket)){
      df <- df %>% add_row(parent = branch$ances[1],left = branch$desc[1],right = branch$desc[2],desc_side = "left",asc_side = "right")
    } else if (all(left_diff_bucket < right_diff_bucket)){
      df <- df %>% add_row(parent = branch$ances[1],left = branch$desc[1],right = branch$desc[2],desc_side = "right",asc_side = "left")
    }
  }

  return(df)
}

PostSearch <- function(tree_backbone,df_pairs,subtree_list,counts,states,state_lineages,labels){
  for (i in 1:nrow(df_pairs)){
    subtree_child_label <- tree_backbone$tip.label[df_pairs[[df_pairs$desc_side[i]]][i]]
    subtree_parent_label <- tree_backbone$tip.label[df_pairs[[df_pairs$asc_side[i]]][i]]

    if (!grepl("cell", subtree_child_label, fixed = TRUE)){
      subtree_child <- subtree_list[[which(names(subtree_list) == subtree_parent_label)]]
      counts_child <- counts[labels %in% subtree_child$tip.label,]
    } else{
      counts_child <- counts[labels==subtree_child_label,]
    }

    if (!grepl("cell", subtree_parent_label, fixed = TRUE)){
      subtree_parent <- subtree_list[[which(names(subtree_list) == subtree_parent_label)]]
      counts_parent <- counts[labels %in% subtree_parent$tip.label,]
    } else{
      counts_parent <- counts[labels==subtree_parent_label,]
    }
  }
}


MergeDT <- function(tree_backbone,df_pairs,dt){
  dt_merge <- dt
  dt_merge$merged <- FALSE
  tree <- makeNodeLabel(tree_backbone)
  drop_tips <- c()
  for (i in 1:nrow(df_pairs)){
    dt_child <- dt[df_pairs[[df_pairs$desc_side[i]]][i],]
    dt_parent <- dt[df_pairs[[df_pairs$asc_side[i]]][i],]

    dt_merge <- dt_merge%>% add_row(barcode = dt_parent$barcode,cellids = list(c(dt_parent$cellids[[1]],dt_child$cellids[[1]])),merged = FALSE)
    dt_merge$merged[c(df_pairs$left[i],df_pairs$right[i])] <- TRUE

    tree <- tree %>% phytools::bind.tip(as.character(nrow(dt_merge)),where = df_pairs$parent[i] + i - 1,edge.length = 1)
    drop_tips <- c(drop_tips,as.character(df_pairs[i,2:3]))
  }

  tree <- tree %>% drop.tip(tip = drop_tips)
  tree <- multi2di(tree)
  tree$edge.length <- rep(1,length(tree$edge.length))

  return(list(tree,dt_merge))
}
