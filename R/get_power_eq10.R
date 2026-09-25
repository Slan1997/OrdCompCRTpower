### power calculation based on equation 10 (ordinal)  
# Rutterford et al. 2015 PMC Methods for sample size determination in cluster randomized trials.pdf
# m: number of subjects per arm = cluster_size*number_clusters/2
get_power_eq10 <- function(log_OR, prop_categories, ICC, cluster_size, m, alpha = 0.05){
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  z_alpha_2 <- qnorm(1 - alpha/2)
  pnorm( sqrt(m * x1 / x2 / x3) - z_alpha_2 ) # log_OR should not be smaller than 0?
}


solve_number_clusters_eq10 <- function(log_OR,
                                       prop_categories,
                                       ICC,
                                       cluster_size,
                                       target_power = 0.80,
                                       alpha = 0.05) {
  if (!is.numeric(log_OR) || length(log_OR) != 1 || !is.finite(log_OR) || log_OR == 0) {
    stop("log_OR must be a finite nonzero scalar.")
  }
  if (!is.numeric(prop_categories) || any(prop_categories < 0) ||
      abs(sum(prop_categories) - 1) > 1e-8) {
    stop("prop_categories must be a probability vector summing to 1.")
  }
  if (!is.numeric(ICC) || length(ICC) != 1 || !is.finite(ICC)) {
    stop("ICC must be a finite scalar.")
  }
  if (!is.numeric(cluster_size) || length(cluster_size) != 1 ||
      !is.finite(cluster_size) || cluster_size <= 0) {
    stop("cluster_size must be a positive scalar.")
  }
  if (!is.numeric(target_power) || length(target_power) != 1 ||
      target_power <= 0 || target_power >= 1) {
    stop("target_power must be in (0, 1).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).")
  }
  
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  
  z_alpha_2 <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(target_power)
  
  m_per_arm <- (x2 * x3 / x1) * (z_alpha_2 + z_beta)^2
  
  number_clusters_cont <- 2 * m_per_arm / cluster_size
  number_clusters_int <- ceiling(number_clusters_cont)
  
  achieved_power <- get_power_eq10(
    log_OR = log_OR,
    prop_categories = prop_categories,
    ICC = ICC,
    cluster_size = cluster_size,
    m = number_clusters_int * cluster_size / 2,
    alpha = alpha
  )
  
  list(
    m_per_arm_continuous = m_per_arm,
    N_required_continuous = number_clusters_cont,
    N_required_integer = number_clusters_int,
    achieved_power = achieved_power
  )
}

# solve_number_clusters_eq10(
#   log_OR = log(1.5),
#   prop_categories = c(0.1, 0.8, 0.1),
#   ICC = 0.01,
#   cluster_size = 50,
#   target_power = 0.80,
#   alpha = 0.05
# )
