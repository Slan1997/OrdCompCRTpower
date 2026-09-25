format_prob_vec_for_copy_5 <- function(probs, digits = 5) {
  probs <- as.numeric(probs)
  
  if (length(probs) < 2 || any(!is.finite(probs)) || any(probs < 0)) {
    stop("probs must be a finite nonnegative vector with at least 2 entries.")
  }
  
  probs <- probs / sum(probs)
  
  scale <- 10^digits
  
  rounded <- numeric(length(probs))
  rounded[-length(probs)] <- round(probs[-length(probs)], digits)
  rounded[length(probs)] <- 1 - sum(rounded[-length(probs)])
  
  # avoid negative residual caused by unlucky rounding
  if (rounded[length(probs)] < 0) {
    floored <- floor(probs * scale) / scale
    remainder_units <- scale - sum(floored * scale)
    
    # distribute remaining units to the largest fractional remainders
    raw_units <- probs * scale
    frac <- raw_units - floor(raw_units)
    add_idx <- order(frac, decreasing = TRUE)[seq_len(remainder_units)]
    
    rounded <- floored
    rounded[add_idx] <- rounded[add_idx] + 1 / scale
  }
  
  paste(
    sprintf(paste0("%.", digits, "f"), rounded),
    collapse = ", "
  )
}