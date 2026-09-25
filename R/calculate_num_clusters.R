# helpers
make_a_mat <- function(pi_vec) {
  K <- length(pi_vec)
  a_vec <- (pi_vec * (1 - pi_vec))^(-1/2)
  diag(a_vec, nrow = K, ncol = K)
}

make_Q_mat_fast <- function(pi_vec) {
  a <- (pi_vec * (1 - pi_vec))^(-1/2)
  
  r <- -outer(pi_vec, pi_vec) * outer(a, a)
  
  diag(r) <- 1
  r
}

make_d_mat <- function(pi_vec, treat_last_col = TRUE) {
  K <- length(pi_vec)
  
  vartheta <- cumsum(pi_vec)
  v <- vartheta * (1 - vartheta)
  
  d <- matrix(0, nrow = K, ncol = K + 1)
  
  for (k in 1:K) {
    d[k, k] <- v[k]
    
    if (k >= 2)
      d[k, k - 1] <- -v[k - 1]
    
    if (treat_last_col) {
      prev <- if (k >= 2) v[k - 1] else 0
      d[k, K + 1] <- v[k] - prev
    }
  }
  d
}

solve_number_cluster <- function(p_C_vec,
                                 theta,
                                 cluster_size,
                                 rho_matrix,
                                 alpha = 0.05,
                                 target_power = 0.80,
                                 n_cat_index = NULL,
                                 round_even = FALSE) {
  # p_C_vec: full category probabilities, length = #categories
  # for 3-cat: c(p1, p2, p3)
  
  if (!is.numeric(p_C_vec) || length(p_C_vec) < 2) {
    stop("p_C_vec must be a numeric probability vector of length >= 2.")
  }
  if (any(!is.finite(p_C_vec)) || any(p_C_vec <= 0) || abs(sum(p_C_vec) - 1) > 1e-8) {
    stop("p_C_vec must have positive entries and sum to 1.")
  }
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("theta must be a single finite number.")
  }
  if (!is.numeric(cluster_size) || length(cluster_size) != 1 || !is.finite(cluster_size) || cluster_size <= 1) {
    stop("cluster_size must be > 1.")
  }
  if (!is.matrix(rho_matrix)) {
    stop("rho_matrix must be a matrix.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).")
  }
  if (!is.numeric(target_power) || length(target_power) != 1 || target_power <= 0 || target_power >= 1) {
    stop("target_power must be in (0, 1).")
  }
  
  # Number of thresholds = number of categories - 1
  K <- length(p_C_vec) - 1
  
  if (any(dim(rho_matrix) != c(K, K))) {
    stop(sprintf("rho_matrix must be %d x %d.", K, K))
  }
  
  if (is.null(n_cat_index)) n_cat_index <- K + 1
  
  # 1) Treatment category probabilities
  p_T_vec <- pi_T_from_PO(p_C_vec, theta)
  
  # 2) Convert to threshold probabilities for the matrix helpers
  pi_C_hat <- p_C_vec[1:K]
  pi_T_hat <- p_T_vec[1:K]
  
  # 3) Build pieces
  a_mat_T <- make_a_mat(pi_T_hat)
  a_mat_C <- make_a_mat(pi_C_hat)
  
  Q_mat_T <- make_Q_mat_fast(pi_T_hat)
  Q_mat_C <- make_Q_mat_fast(pi_C_hat)
  
  d_mat_T <- make_d_mat(pi_T_hat, treat_last_col = TRUE)
  d_mat_C <- make_d_mat(pi_C_hat, treat_last_col = FALSE)
  
  inv_part_T <- solve(Q_mat_T + (cluster_size - 1) * rho_matrix)
  inv_part_C <- solve(Q_mat_C + (cluster_size - 1) * rho_matrix)
  
  dvd_T <- (t(d_mat_T) %*% a_mat_T %*% inv_part_T %*% a_mat_T) %*% d_mat_T
  dvd_C <- (t(d_mat_C) %*% a_mat_C %*% inv_part_C %*% a_mat_C) %*% d_mat_C
  
  M <- solve(dvd_T + dvd_C)
  v0 <- (2 / cluster_size) * M[n_cat_index, n_cat_index]
  
  if (!is.finite(v0) || v0 <= 0) {
    stop("Computed v0 <= 0. Check rho_matrix, Q/A matrices, and indexing.")
  }
  
  z_alpha_2 <- qnorm(1 - alpha / 2)
  z_pow <- qnorm(target_power)
  z_sum <- z_alpha_2 + z_pow
  
  N_req <- v0 * (z_sum^2) / (theta^2)
  
  N_int <- ceiling(N_req)
  if (round_even && (N_int %% 2 == 1)) N_int <- N_int + 1
  
  var_theta_hat <- v0 / N_int
  achieved_power <- pnorm(abs(theta) / sqrt(var_theta_hat) - z_alpha_2)
  
  list(
    N_required_continuous = N_req,
    N_required_integer = N_int,
    achieved_power = achieved_power,
    v0 = v0,
    p_T_vec = p_T_vec,
    p_C_vec = p_C_vec
  )
}

solve_clusters_binary <- function(P0, theta_R, rho, cluster_size,
                                  alpha = 0.05, target_power = 0.80,
                                  allocate_prop_C = 0.5,
                                  round_even = F) {
  OR <- exp(theta_R)
  P1 <- OR * P0 / (1 - P0 + OR * P0)
  # print(P0)
  # print(P1)
  DE <- 1 + (cluster_size - 1) * rho
  var_term <- 1/(allocate_prop_C * P0 * (1 - P0)) +
    1/((1 - allocate_prop_C) * P1 * (1 - P1))
  
  z_alpha_2 <- qnorm(1 - alpha/2)
  z_pow     <- qnorm(target_power)
  z_sum     <- z_alpha_2 + z_pow
  
  m_cont <- (DE * var_term / (cluster_size * theta_R^2)) * (z_sum^2)
  
  m_int <- ceiling(m_cont)
  if (round_even && (m_int %% 2 == 1)) m_int <- m_int + 1
  
  final_term1 <- sqrt(cluster_size * m_int * theta_R^2 / (DE * var_term))
  achieved_power <- pnorm(final_term1 - z_alpha_2)
  
  list(
    N_required_continuous = m_cont,
    N_required_integer = m_int,
    achieved_power = achieved_power,
    v0 = (DE * var_term) / cluster_size,   # so Var(theta_hat)= v0 / m
    pi_T_hat = P1 ,                         # treatment prob
    pi_C_hat = P0                      # control prob
  )
}




