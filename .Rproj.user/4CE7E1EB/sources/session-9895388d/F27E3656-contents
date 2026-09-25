power_binary = function(P0, P1, theta_R, rho,alpha = 0.05,
                        cluster_size,number_cluster,allocate_prop_C=.5){
  # basic checks
  if (any(c(P0, P1) <= 0) || any(c(P0, P1) >= 1)) stop("P0 and P1 must be in (0,1).")
  if (!is.numeric(theta_R) || length(theta_R) != 1 || is.na(theta_R)) stop("theta_R must be scalar numeric.")
  if (!is.numeric(rho) || length(rho) != 1 || is.na(rho)) stop("rho must be scalar numeric.")
  if (cluster_size <= 1) stop("cluster_size must be > 1.")
  if (number_cluster <= 1) stop("number_cluster must be > 1.")
  if (allocate_prop_C <= 0 || allocate_prop_C >= 1) stop("allocate_prop_C must be in (0,1).")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0,1).")
  DE = 1+ (cluster_size-1)*rho
  var_term = 1 / ( allocate_prop_C*P0*(1-P0) ) + 1 / ( (1-allocate_prop_C)*P1*(1-P1) )
  ex_V = DE/cluster_size * var_term / number_cluster
  z_alpha_2 <- qnorm(1 - alpha/2)
  final_term1_numerator = cluster_size*number_cluster*theta_R^2
  final_term1_denominator = DE*var_term
  final_term1 = sqrt(final_term1_numerator/final_term1_denominator)
  pnorm(final_term1 - z_alpha_2  )
}