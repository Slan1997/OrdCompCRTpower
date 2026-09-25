# assumed 1:1 allocation rate for now.
# assumed proportional odds.
# number_cluster: can be a vector, total number of clusters.
# cluster_size
# treatment_effect: marginal log odds ratio, can be null if both p_C_vec and p_T_vec are not null?
# p_C_vec: control group distributions, 
# p_T_vec: can be null and calculated with treatment_effect and p_C_vec using function pi_T_from_PO().
# rho11
# rho12
# rho22
# alpha = 0.05,

calculate_power_3cat <- function(number_cluster,
                                 cluster_size,
                                 p_C_vec = NULL,
                                 p_T_vec = NULL,
                                 pi_vec = NULL,
                                 treatment_effect = NULL,
                                 rho11, rho12, rho22,
                                 add_bin = TRUE,
                                 bin_rho1 = NULL,
                                 bin_rho2 = NULL,
                                 add_WH_DE = TRUE,
                                 WH_DE_ICC = NULL,
                                 add_indep = TRUE,
                                 alpha = 0.05,
                                 tol = 1e-4) {
  
  # Resolve whichever two inputs the user supplied
  resolved <- resolve_3cat_inputs(
    p_C_vec = p_C_vec,
    p_T_vec = p_T_vec,
    pi_vec = pi_vec,
    treatment_effect = treatment_effect,
    tol = tol
  )
  
  p_C_vec <- resolved$p_C_vec
  p_T_vec <- resolved$p_T_vec
  pi_vec  <- resolved$pi_vec
  treatment_effect <- resolved$treatment_effect
  
  # Basic checks
  if (!is.numeric(number_cluster) || length(number_cluster) != 1 || 
      !is.finite(number_cluster) || number_cluster <= 1) {
    stop("number_cluster must be a single numeric value > 1.")
  }
  
  if (!is.numeric(cluster_size) || length(cluster_size) != 1 || 
      !is.finite(cluster_size) || cluster_size <= 1) {
    stop("cluster_size must be a single numeric value > 1.")
  }
  
  if (!all(sapply(list(rho11, rho12, rho22), function(x) {
    is.numeric(x) && length(x) == 1 && is.finite(x)
  }))) {
    stop("rho11, rho12, rho22 must be finite scalar values.")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }
  
  rho_matrix <- matrix(c(rho11, rho12,
                         rho12, rho22), nrow = 2, byrow = TRUE)
  
  if (is.null(bin_rho1)) bin_rho1 <- rho11
  if (is.null(bin_rho2)) {
    bin_rho2 <- rho_binary2_from_p(
      p = pi_vec,
      rho11 = rho11,
      rho12 = rho12,
      rho22 = rho22
    )
  }
  
  # category probabilities
  p0_C <- p_C_vec[1]
  p1_C <- p_C_vec[2]
  p0_T <- p_T_vec[1]
  p1_T <- p_T_vec[2]
  
  ## GEE ordinal exchangeable
  ex_V <- exch_V_theta(
    pi_T1 = p0_T, pi_T2 = p1_T,
    pi_C1 = p0_C, pi_C2 = p1_C,
    rho_mat = rho_matrix,
    cluster_size = cluster_size,
    number_cluster = number_cluster
  )
  
  power_GEE_exch <- power_gee_est(
    V_theta = ex_V,
    theta_star = treatment_effect,
    alpha = alpha
  )
  
  ## Binary collapses
  power_GEE_bin1 <- NULL
  power_GEE_bin2 <- NULL
  
  logit <- function(p) log(p / (1 - p))
  
  theta_bin1 <- logit(p0_T) - logit(p0_C)
  theta_bin2 <- logit(p0_T + p1_T) - logit(p0_C + p1_C)
  
  if (add_bin) {
    # Binary1: first category vs remaining
    power_GEE_bin1 <- power_binary(
      P0 = p0_C,
      P1 = p0_T,
      theta_R = theta_bin1,
      rho = bin_rho1,
      cluster_size = cluster_size,
      number_cluster = number_cluster,
      allocate_prop_C = 0.5,
      alpha = alpha
    )
    
    # Binary2: first two categories vs last
    power_GEE_bin2 <- power_binary(
      P0 = p0_C + p1_C,
      P1 = p0_T + p1_T,
      theta_R = theta_bin2,
      rho = bin_rho2,
      cluster_size = cluster_size,
      number_cluster = number_cluster,
      allocate_prop_C = 0.5,
      alpha = alpha
    )
  }
  
  
  if (is.null(WH_DE_ICC)) WH_DE_ICC <- rho11
  power_WH_DE <- NA
  
  if (add_WH_DE) {
    power_WH_DE <- get_power_eq10(
      log_OR = treatment_effect,
      prop_categories = pi_vec,
      ICC = WH_DE_ICC,
      cluster_size = cluster_size,
      m = (number_cluster*cluster_size/2),
      alpha = alpha
    )
  }
  
  
  ## Independence + WH
  power_GEE_indep <- NA
  power_WH <- NA
  
  if (add_indep) {
    indep_V <- exch_V_theta(
      pi_T1 = p0_T, pi_T2 = p1_T,
      pi_C1 = p0_C, pi_C2 = p1_C,
      rho_mat = matrix(0, 2, 2),
      cluster_size = cluster_size,
      number_cluster = number_cluster
    )
    
    power_GEE_indep <- power_gee_est(
      V_theta = indep_V,
      theta_star = treatment_effect,
      alpha = alpha
    )
    
    power_WH <- wh_power(
      theta_R = treatment_effect,
      n = number_cluster * cluster_size,
      pr = pi_vec,
      A = 1,
      alpha = alpha
    )
  }
  
  list(
    power_GEE_exch = power_GEE_exch,
    power_GEE_bin1 = power_GEE_bin1,
    power_GEE_bin2 = power_GEE_bin2,
    power_WH_DE = power_WH_DE,
    power_GEE_indep = power_GEE_indep,
    power_WH = power_WH,
    p_C_vec = p_C_vec,
    p_T_vec = p_T_vec,
    pi_vec = pi_vec,
    treatment_effect = treatment_effect,
    bin_rho1 = bin_rho1,
    bin_rho2 = bin_rho2,
    theta_bin1 = theta_bin1,
    theta_bin2 = theta_bin2,
    WH_DE_ICC = WH_DE_ICC,
    specified = resolved$specified
  )
}


# examples: 
# calculate_power_3cat(
#   number_cluster = 20,
#   cluster_size = 50,
#   p_C_vec = c(0.1, 0.8, 0.1),
#   treatment_effect = log(1.5),
#   rho11 = 0.005,
#   rho12 = -0.005,
#   rho22 = 0.005
# )







