compute_clusters_assumed_df <- function(p_C_vec,
                                        or_grid,
                                        cluster_size,
                                        alpha,
                                        target_power,
                                        rho_entry_grid,
                                        add_bin = TRUE,
                                        round_even = FALSE) {
  if (!is.numeric(p_C_vec) || length(p_C_vec) != 3) {
    stop("p_C_vec must be a numeric vector of length 3.")
  }
  if (any(!is.finite(p_C_vec)) || any(p_C_vec <= 0) || abs(sum(p_C_vec) - 1) > 1e-8) {
    stop("p_C_vec must have positive entries and sum to 1.")
  }
  if (!is.numeric(or_grid) || length(or_grid) == 0 || any(!is.finite(or_grid)) || any(or_grid <= 0)) {
    stop("or_grid must be a non-empty numeric vector of positive values.")
  }
  if (!is.numeric(cluster_size) || length(cluster_size) != 1 || !is.finite(cluster_size) || cluster_size <= 1) {
    stop("cluster_size must be > 1.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).")
  }
  if (!is.numeric(target_power) || length(target_power) != 1 || target_power <= 0 || target_power >= 1) {
    stop("target_power must be in (0, 1).")
  }
  if (!is.numeric(rho_entry_grid) || length(rho_entry_grid) == 0 || any(!is.finite(rho_entry_grid))) {
    stop("rho_entry_grid must be a non-empty numeric vector.")
  }
  
  out_list <- list()
  
  for (or_val in or_grid) {
    theta <- log(or_val)
    
    for (rho_entry in rho_entry_grid) {
      res <- list()
      
      # GEE-Ordinal
      ord_out <- #tryCatch(
        solve_number_cluster(
          p_C_vec = p_C_vec,
          theta = theta,
          cluster_size = cluster_size,
          rho_matrix = build_rho_assumed(3, rho_entry),
          alpha = alpha,
          target_power = target_power,
          round_even = round_even
        #),
       # error = function(e) NULL
      )
      
      if (!is.null(ord_out)) {
        res[[length(res) + 1]] <- data.frame(
          method = "Ordinal-3cat",
          OR = or_val,
          OR_label = sprintf("log(%s)", format(round(or_val, 2), nsmall = 2, trim = TRUE)),
          rho_entry = rho_entry,
          rho_label = sprintf("\u03c1 = %s", format(round(rho_entry, 3), nsmall = 3, trim = TRUE)),
          N_continuous = ord_out$N_required_continuous,
          N_integer = ord_out$N_required_integer,
          achieved_power = 100 * ord_out$achieved_power,
          check.names = FALSE
        )
      }
      
      if (isTRUE(add_bin)) {
        # Binary1
        bin1_out <- tryCatch(
          solve_clusters_binary(
            P0 = as.numeric(p_C_vec[1]),
            theta_R = theta,
            rho = rho_entry,
            cluster_size = cluster_size,
            alpha = alpha,
            target_power = target_power,
            round_even = round_even
          ),
          error = function(e) NULL
        )
        
        if (!is.null(bin1_out)) {
          res[[length(res) + 1]] <- data.frame(
            method = "Binary1",
            OR = or_val,
            OR_label = sprintf("log(%s)", format(round(or_val, 2), nsmall = 2, trim = TRUE)),
            rho_entry = rho_entry,
            rho_label = sprintf("\u03c1 = %s", format(round(rho_entry, 3), nsmall = 3, trim = TRUE)),
            N_continuous = bin1_out$N_required_continuous,
            N_integer = bin1_out$N_required_integer,
            achieved_power = 100 * bin1_out$achieved_power,
            check.names = FALSE
          )
        }
        
        # Binary2
        bin2_out <- tryCatch(
          solve_clusters_binary(
            P0 = as.numeric(sum(p_C_vec[1:2])),
            theta_R = theta,
            rho = rho_entry,
            cluster_size = cluster_size,
            alpha = alpha,
            target_power = target_power,
            round_even = round_even
          ),
          error = function(e) NULL
        )
        
        if (!is.null(bin2_out)) {
          res[[length(res) + 1]] <- data.frame(
            method = "Binary2",
            OR = or_val,
            OR_label = sprintf("log(%s)", format(round(or_val, 2), nsmall = 2, trim = TRUE)),
            rho_entry = rho_entry,
            rho_label = sprintf("\u03c1 = %s", format(round(rho_entry, 3), nsmall = 3, trim = TRUE)),
            N_continuous = bin2_out$N_required_continuous,
            N_integer = bin2_out$N_required_integer,
            achieved_power = 100 * bin2_out$achieved_power,
            check.names = FALSE
          )
        }
      }
      
      if (length(res) > 0) {
        out_list[[length(out_list) + 1]] <- do.call(rbind, res)
      }
    }
  }
  
  if (length(out_list) == 0) {
    return(data.frame())
  }
  
  df <- do.call(rbind, out_list)
  rownames(df) <- NULL
  
  method_levels <- c("Binary1", "Binary2", "Ordinal-3cat")
  df$method <- factor(df$method, levels = method_levels[method_levels %in% unique(df$method)])
  
  df$OR_label <- factor(
    df$OR_label,
    levels = sprintf("log(%s)", format(round(or_grid, 2), nsmall = 2, trim = TRUE))
  )
  
  df
}

# # examples
# p_C_vec <- c(0.876, 0.0465, 0.0775)
# or_grid <- c(1.5, 1.75, 2.0)
# cluster_size <- 50
# alpha <- 0.05
# target_power <- 0.80
# rho_entry_grid <- c(0.006, 0.009, 0.012)
# 
# # Run the function
# df_clusters <- compute_clusters_assumed_df(
#   p_C_vec = p_C_vec,
#   or_grid = or_grid,
#   cluster_size = cluster_size,
#   alpha = alpha,
#   target_power = target_power,
#   rho_entry_grid = rho_entry_grid,
#   add_bin = TRUE,
#   round_even = FALSE
# )
# 
# # View results
# df_clusters

