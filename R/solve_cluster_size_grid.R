solve_cluster_size_grid <- function(p_C_vec,
                                    theta,
                                    n_clusters,
                                    rho_matrix,
                                    alpha = 0.05,
                                    target_power = 0.80,
                                    n_cat_index = NULL,
                                    m_grid = 2:500) {

  res <- lapply(m_grid, function(m) {
    out <- solve_number_cluster(
      p_C_vec = p_C_vec,
      theta = theta,
      cluster_size = m,
      rho_matrix = rho_matrix,
      alpha = alpha,
      target_power = target_power,
      n_cat_index = n_cat_index,
      round_even = FALSE
    )

    data.frame(
      cluster_size = m,
      N_required = out$N_required_continuous,
      achieved_power = out$achieved_power
    )
  })

  res_df <- do.call(rbind, res)

  feasible <- res_df[res_df$N_required <= n_clusters, , drop = FALSE]

  if (nrow(feasible) == 0) {
    return(list(
      message = "No feasible cluster size found in the search range.",
      results = res_df
    ))
  }

  best <- feasible[1, , drop = FALSE]

  list(
    cluster_size_integer = best$cluster_size,
    n_clusters_target = n_clusters,
    achieved_clusters_required = best$N_required,
    achieved_power = best$achieved_power,
    results = res_df
  )
}