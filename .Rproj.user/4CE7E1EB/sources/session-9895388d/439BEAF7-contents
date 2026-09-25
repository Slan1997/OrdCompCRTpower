generate_dvine_pairs <- function(n_comp) {
  out <- vector("list", n_comp - 1)
  
  for (tree in 1:(n_comp - 1)) {
    n_pairs <- n_comp - tree
    pair_labels <- character(n_pairs)
    
    for (i in 1:n_pairs) {
      j <- i + tree
      
      if (tree == 1) {
        pair_labels[i] <- paste0("(", i, ",", j, ")")
      } else {
        cond <- paste((i + 1):(j - 1), collapse = ",")
        pair_labels[i] <- paste0("(", i, ",", j, "|", cond, ")")
      }
    }
    
    out[[tree]] <- data.frame(
      Tree = tree,
      Pair = pair_labels,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}


generate_cvine_pairs <- function(n_comp) {
  out <- vector("list", n_comp - 1)
  
  for (tree in 1:(n_comp - 1)) {
    root <- tree
    others <- (tree + 1):n_comp
    pair_labels <- character(length(others))
    
    for (k in seq_along(others)) {
      j <- others[k]
      
      if (tree == 1) {
        pair_labels[k] <- paste0("(", root, ",", j, ")")
      } else {
        cond <- paste(1:(tree - 1), collapse = ",")
        pair_labels[k] <- paste0("(", root, ",", j, "|", cond, ")")
      }
    }
    
    out[[tree]] <- data.frame(
      Tree = tree,
      Pair = pair_labels,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

generate_vine_pairs <- function(n_comp, vine_type = "D") {
  if (!is.numeric(n_comp) || length(n_comp) != 1 || is.na(n_comp) || n_comp < 2) {
    stop("n_comp must be a single integer >= 2.")
  }
  
  n_comp <- as.integer(n_comp)
  vine_type <- toupper(trimws(vine_type))
  
  if (vine_type == "D") {
    return(generate_dvine_pairs(n_comp))
  } else if (vine_type == "C") {
    return(generate_cvine_pairs(n_comp))
  } else {
    stop("vine_type must be either 'D' or 'C'.")
  }
}
