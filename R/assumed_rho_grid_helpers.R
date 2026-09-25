build_rho_assumed <- function(n_cat, rho_entry) {
  if (n_cat == 2) {
    matrix(rho_entry, nrow = 1, ncol = 1)
  } else if (n_cat == 3) {
    matrix(c(rho_entry, -rho_entry,
             -rho_entry, rho_entry),
           nrow = 2, byrow = TRUE)
  } else {
    stop("This grid is only intended for binary (n_cat = 2) and ordinal-3cat (n_cat = 3).")
  }
}

parse_rho_grid <- function(mode = c("sequence", "custom"),
                           rho_min = NULL,
                           rho_max = NULL,
                           rho_step = NULL,
                           rho_custom = NULL,
                           tol = 1e-12) {
  mode <- match.arg(mode)
  
  if (mode == "sequence") {
    if (!is.numeric(rho_min) || length(rho_min) != 1 || !is.finite(rho_min)) {
      stop("rho_min must be a single finite number.")
    }
    if (!is.numeric(rho_max) || length(rho_max) != 1 || !is.finite(rho_max)) {
      stop("rho_max must be a single finite number.")
    }
    if (!is.numeric(rho_step) || length(rho_step) != 1 || !is.finite(rho_step) || rho_step <= 0) {
      stop("rho_step must be a single positive finite number.")
    }
    if (rho_min > rho_max) {
      stop("rho_min must be less than or equal to rho_max.")
    }
    
    grid <- seq(rho_min, rho_max, by = rho_step)
    grid <- unique(round(grid, 10))
    return(grid)
  }
  
  if (is.null(rho_custom) || !nzchar(trimws(rho_custom))) {
    stop("Please provide at least one custom rho entry.")
  }
  
  parts <- unlist(strsplit(rho_custom, ","))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  
  vals <- suppressWarnings(as.numeric(parts))
  
  if (length(vals) == 0 || any(is.na(vals))) {
    stop("Custom rho entries must be comma-separated numeric values.")
  }
  
  vals <- unique(round(vals, 10))
  vals <- sort(vals)
  vals
}

