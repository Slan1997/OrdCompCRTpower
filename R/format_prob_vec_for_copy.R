format_prob_vec_for_copy <- function(probs) {
  probs <- as.numeric(probs)
  
  paste(
    format(
      probs,
      digits = 16,
      scientific = FALSE,
      trim = TRUE
    ),
    collapse = ", "
  )
}