get_pair_label_table <- function(nComp, vine_type) {
  tmp <- VCO_Initialize(
    nComp = nComp,
    prev = rep(0.2, nComp),
    vine_type = vine_type,
    method_list = rep("Gaussian", nComp - 1),
    param_list = rep(0, nComp - 1)
  )
  
  out <- list()
  
  for (tt in seq_len(tmp$nTree)) {
    pair_mat <- tmp$pair_list[[tt]]
    n_pair <- ncol(pair_mat)
    
    pair_labels <- character(n_pair)
    
    for (j in seq_len(n_pair)) {
      a <- pair_mat[1, j]
      b <- pair_mat[2, j]
      
      if (tt == 1) {
        pair_labels[j] <- paste0("(", a, ",", b, ")")
      } else {
        cond_node <- gsub("^V", "", tmp$conditioned_on_nodes[[tt - 1]][j])
        cond_digits <- strsplit(cond_node, "")[[1]]
        pair_labels[j] <- paste0("(", a, ",", b, "|", paste(cond_digits, collapse = ","), ")")
      }
    }
    
    out[[tt]] <- data.frame(
      Tree = tt,
      Pair = pair_labels,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}
