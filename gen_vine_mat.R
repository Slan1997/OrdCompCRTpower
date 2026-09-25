generate_d_vine_matrix <- function(d) {
  # Initialize an empty d x d matrix
  mat <- matrix(0, nrow = d, ncol = d)
  
  for (i in 1:d) {
    # Set the diagonal element
    mat[i, i] <- i
    
    # Fill the lower triangle (if not the first row)
    if (i > 1) {
      for (j in 1:(i - 1)) {
        # The pattern: M[i, j] = i - j
        mat[i, j] <- i - j
      }
    }
  }
  return(t(mat))
}

generate_c_vine_matrix <- function(d) {
  # Initialize a d x d matrix with zeros
  mat <- matrix(0, nrow = d, ncol = d)
  
  for (i in 1:d) {
    # Each row i is filled with the sequence 1 to i
    # This matches your logic: Row 1 is [1], Row 2 is [1, 2], etc.
    mat[i, 1:i] <- 1:i
  }
  
  return(t(mat))
}
