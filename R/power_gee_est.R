## ---------- GEE ordinal (3-cat) building block ----------
gee_info_mat_3cat <- function(pi_T1, pi_T2, pi_C1, pi_C2,
                              rho_mat, cluster_size) {
  if (cluster_size <= 1) stop("cluster_size must be > 1.")
  if (!is.matrix(rho_mat) || any(dim(rho_mat) != c(2, 2))) {
    stop("rho_mat must be 2x2 matrix.")
  }
  
  cdf_T2 <- pi_T1 + pi_T2
  cdf_C2 <- pi_C1 + pi_C2
  
  a_T1 <- (pi_T1 * (1 - pi_T1))^(-1/2)
  a_T2 <- (pi_T2 * (1 - pi_T2))^(-1/2)
  a_C1 <- (pi_C1 * (1 - pi_C1))^(-1/2)
  a_C2 <- (pi_C2 * (1 - pi_C2))^(-1/2)
  
  a_mat_T <- matrix(c(a_T1, 0, 0, a_T2), nrow = 2)
  a_mat_C <- matrix(c(a_C1, 0, 0, a_C2), nrow = 2)
  
  r_T <- -pi_T1 * pi_T2 * a_T1 * a_T2
  r_C <- -pi_C1 * pi_C2 * a_C1 * a_C2
  
  Q_mat_T <- matrix(c(1, r_T,
                      r_T, 1), nrow = 2, byrow = TRUE)
  Q_mat_C <- matrix(c(1, r_C,
                      r_C, 1), nrow = 2, byrow = TRUE)
  
  d_mat_T <- matrix(c(
    pi_T1 * (1 - pi_T1), -pi_T1 * (1 - pi_T1),
    0,                    cdf_T2 * (1 - cdf_T2),
    pi_T1 * (1 - pi_T1),  cdf_T2 * (1 - cdf_T2) - pi_T1 * (1 - pi_T1)
  ), ncol = 3)
  
  d_mat_C <- matrix(c(
    pi_C1 * (1 - pi_C1), -pi_C1 * (1 - pi_C1),
    0,                    cdf_C2 * (1 - cdf_C2),
    0,                    0
  ), ncol = 3)
  
  inv_part_T <- solve(Q_mat_T + (cluster_size - 1) * rho_mat)
  inv_part_C <- solve(Q_mat_C + (cluster_size - 1) * rho_mat)
  
  DVD_T <- cluster_size * t(d_mat_T) %*% a_mat_T %*% inv_part_T %*% a_mat_T %*% d_mat_T
  DVD_C <- cluster_size * t(d_mat_C) %*% a_mat_C %*% inv_part_C %*% a_mat_C %*% d_mat_C
  
  DVD_T + DVD_C
}

is_valid_gee_info_mat <- function(M, tol = 1e-10) {
  if (!is.matrix(M) || nrow(M) != ncol(M)) return(FALSE)
  if (any(!is.finite(M))) return(FALSE)
  
  ev <- tryCatch(
    eigen(M, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  
  all(is.finite(ev)) && all(ev > tol)
}

# Returns V_theta (variance of theta-hat) for exchangeable / independence working corr
exch_V_theta = function(pi_T1, pi_T2, pi_C1, pi_C2,
                        rho_mat, # 2 by 2 # if rho_mat is 2 by 2 0, this is independence.
                        # J: cluster size 
                        cluster_size,number_cluster){
  # checks
  if (cluster_size <= 1) stop("cluster_size must be > 1.")
  if (number_cluster <= 1) stop("number_cluster must be > 1.")
  if (!is.matrix(rho_mat) || any(dim(rho_mat) != c(2, 2))) stop("rho_mat must be 2x2 matrix.")
  
  M <- gee_info_mat_3cat(
    pi_T1 = pi_T1, pi_T2 = pi_T2,
    pi_C1 = pi_C1, pi_C2 = pi_C2,
    rho_mat = rho_mat,
    cluster_size = cluster_size
  )
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(ev)) || any(ev <= 1e-10)) {
    stop("DVD_T + DVD_C is not positive definite.")
  }
  
  V_mat <- (2 / number_cluster) * solve(M)
  V_mat[3, 3]
}


power_gee_est = function(alpha=0.05,V_theta,theta_star){
  if (!is.numeric(V_theta) || length(V_theta) != 1 || is.na(V_theta) || V_theta <= 0) {
    stop("V_theta must be a positive scalar.")
  }
  if (!is.numeric(theta_star) || length(theta_star) != 1 || is.na(theta_star)) {
    stop("theta_star must be a scalar numeric.")
  }
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0,1).")
  
  z_alpha_2 <- qnorm(1 - alpha/2)
  pnorm(sqrt(theta_star^2/V_theta) - z_alpha_2)
}
