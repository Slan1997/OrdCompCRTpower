classify_ordinal_distribution <- function(p) {
  if (length(p) != 3) stop("Input distribution must have length 3.")
  if (any(!is.finite(p)) || any(p < 0)) stop("Probabilities must be nonnegative and finite.")
  if (abs(sum(p) - 1) > 1e-8) stop("Probabilities must sum to 1.")
  
  templates <- rbind(
    Bell         = c(0.10, 0.80, 0.10),
    Common       = c(0.10, 0.10, 0.80),
    Rare         = c(0.80, 0.10, 0.10),
    U_shape      = c(0.45, 0.10, 0.45),
    High01low2   = c(0.45, 0.45, 0.10),
    Even         = c(0.33, 0.33, 0.34)
  )
  
  dists <- apply(templates, 1, function(tpl) {
    sqrt(sum((p - tpl)^2))
  })
  
  best_name <- names(which.min(dists))
  best_dist <- min(dists)
  
  # similarity score on a 0-100 scale
  # max possible Euclidean distance in a 3-category probability simplex is sqrt(2)
  similarity <- 100 * (1 - best_dist / sqrt(2))
  similarity <- max(0, min(100, similarity))
  
  # prettier label for UI
  pretty_name <- c(
    Bell = "Bell",
    Common = "Common",
    Rare = "Rare",
    U_shape = "U-shape",
    High01low2 = "High01low2",
    Even = "Even"
  )[best_name]
  
  list(
    assigned_type = unname(pretty_name),
    assigned_code = best_name,
    distance = best_dist,
    similarity = similarity,
    distances = sort(dists)
  )
}