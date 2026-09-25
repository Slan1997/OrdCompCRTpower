# ============================================================
# Helpers for CRT 4+ Cat. tab
# ============================================================

parse_crtK_dist_text <- function(txt, label_prefix = "Distribution", tol = 1e-8) {
  x <- unlist(strsplit(txt, ","))
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- as.numeric(x)
  
  shiny::validate(
    need(length(x) >= 4, paste0(label_prefix, " must contain at least 4 probabilities.")),
    need(length(x) <= 7, paste0(label_prefix, " must contain at most 7 probabilities.")),
    need(all(is.finite(x)), paste0(label_prefix, " entries must all be numeric.")),
    need(all(x > 0), paste0(label_prefix, " entries must all be > 0.")),
    need(abs(sum(x) - 1) < tol, paste0(label_prefix, " must sum to 1."))
  )
  
  x
}

default_3cat_collapse_text <- function(K) {
  if (K < 4) stop("K must be at least 4.")
  paste0("1; 2; ", paste(3:K, collapse = ","))
}

default_binary_collapse_text <- function(K) {
  if (K < 2) stop("K must be at least 2.")
  paste0("1; ", paste(2:K, collapse = ","))
}

parse_collapse_text <- function(txt, K, G = 3) {
  if (is.null(txt) || !nzchar(trimws(txt))) {
    stop("Collapse specification cannot be empty.")
  }
  
  group_txt <- strsplit(txt, ";")[[1]]
  group_txt <- trimws(group_txt)
  group_txt <- group_txt[nzchar(group_txt)]
  
  if (length(group_txt) != G) {
    stop(sprintf("Collapse specification must contain exactly %d groups separated by semicolons.", G))
  }
  
  groups <- lapply(group_txt, function(g) {
    vals <- strsplit(g, ",")[[1]]
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    
    if (length(vals) == 0) {
      stop("Each collapsed group must contain at least one category.")
    }
    
    suppressWarnings(as.integer(vals))
  })
  
  if (any(unlist(lapply(groups, function(g) any(is.na(g)))))) {
    stop("Collapse specification must contain only integer category labels.")
  }
  
  all_vals <- unlist(groups)
  
  if (any(all_vals < 1 | all_vals > K)) {
    stop(sprintf("Category labels must be integers between 1 and %d.", K))
  }
  
  if (any(duplicated(all_vals))) {
    stop("Each original category can appear only once.")
  }
  
  missing_vals <- setdiff(seq_len(K), all_vals)
  
  if (length(missing_vals) > 0) {
    stop(
      paste0(
        "Collapse specification is missing category/categories: ",
        paste(missing_vals, collapse = ", ")
      )
    )
  }
  
  expected_start <- 1
  
  for (g in groups) {
    g_sorted <- sort(g)
    
    if (!identical(g_sorted, seq(min(g_sorted), max(g_sorted)))) {
      stop("Categories within each group must be adjacent.")
    }
    
    if (min(g_sorted) != expected_start) {
      stop("Groups must be adjacent and ordered, e.g., 1,2; 3,4; 5,6.")
    }
    
    expected_start <- max(g_sorted) + 1
  }
  
  if (expected_start != K + 1) {
    stop("Groups must cover all categories in order.")
  }
  
  groups
}

collapse_probs_by_groups <- function(p, groups) {
  as.numeric(sapply(groups, function(g) sum(p[g])))
}

collapse_label_from_groups <- function(groups) {
  paste(
    sapply(groups, function(g) {
      paste0("[", paste(g, collapse = ","), "]")
    }),
    collapse = " | "
  )
}

safe_logit <- function(p, eps = 1e-10) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}


# -----------------------------
# Binary power wrapper using user's power_binary()
# -----------------------------
power_binary_crt_custom <- function(p_C_event,
                                    p_T_event,
                                    theta_R,
                                    rho,
                                    cluster_size,
                                    number_cluster,
                                    alpha = 0.05) {
  
  shiny::validate(
    need(
      exists("power_binary"),
      "power_binary() not found. Did you source your R/ files?"
    ),
    need(
      is.function(power_binary),
      "power_binary exists but is not a function."
    )
  )
  
  out <- tryCatch(
    power_binary(
      P0 = p_C_event,
      P1 = p_T_event,
      theta_R = theta_R,
      rho = rho,
      alpha = alpha,
      cluster_size = cluster_size,
      number_cluster = number_cluster,
      allocate_prop_C = 0.5
    ),
    error = function(e) NA_real_
  )
  
  as.numeric(out)
}


