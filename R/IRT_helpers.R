# ============================================================
# General ordinal helpers for IRT tab
# ============================================================
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

parse_ordinal_dist_text <- function(txt, label_prefix = "Distribution", tol = 1e-8) {
  x <- unlist(strsplit(txt, ","))
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- as.numeric(x)
  
  shiny::validate(
    need(length(x) >= 2, paste0(label_prefix, " must contain at least 2 probabilities.")),
    need(all(is.finite(x)), paste0(label_prefix, " entries must all be numeric.")),
    need(all(x > 0), paste0(label_prefix, " entries must all be > 0.")),
    need(abs(sum(x) - 1) < tol, paste0(label_prefix, " must sum to 1."))
  )
  
  x
}

.check_prob_vec_K <- function(p, name, K = NULL, tol = 1e-8) {
  if (is.null(p)) return(invisible(NULL))
  
  if (!is.numeric(p) || any(!is.finite(p))) {
    stop(sprintf("%s must be a numeric vector.", name))
  }
  
  if (!is.null(K) && length(p) != K) {
    stop(sprintf("%s must be a numeric vector of length %d.", name, K))
  }
  
  if (length(p) < 2) {
    stop(sprintf("%s must contain at least 2 categories.", name))
  }
  
  if (any(p <= 0) || abs(sum(p) - 1) > tol) {
    stop(sprintf("%s must have strictly positive entries and sum to 1.", name))
  }
  
  invisible(NULL)
}

.prob_to_cum_K <- function(p) {
  cumsum(p)[1:(length(p) - 1)]
}

.cum_to_prob_K <- function(v) {
  c(v[1], diff(v), 1 - v[length(v)])
}

pi_T_from_PO_K <- function(pi_C, theta) {
  .check_prob_vec_K(pi_C, "pi_C")
  
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("theta must be a single finite numeric value.")
  }
  
  vartheta_C <- .prob_to_cum_K(pi_C)
  vartheta_T <- expit(logit(vartheta_C) + theta)
  pi_T <- .cum_to_prob_K(vartheta_T)
  
  if (any(pi_T <= 0)) {
    warning("Computed pi_T has non-positive entries; check inputs.")
  }
  
  pi_T
}

pi_C_from_PO_K <- function(pi_T, theta) {
  .check_prob_vec_K(pi_T, "pi_T")
  
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("theta must be a single finite numeric value.")
  }
  
  vartheta_T <- .prob_to_cum_K(pi_T)
  vartheta_C <- expit(logit(vartheta_T) - theta)
  pi_C <- .cum_to_prob_K(vartheta_C)
  
  if (any(pi_C <= 0)) {
    warning("Computed pi_C has non-positive entries; check inputs.")
  }
  
  pi_C
}

theta_from_probs_K <- function(p_C_vec, p_T_vec, tol = 1e-4) {
  .check_prob_vec_K(p_C_vec, "p_C_vec")
  .check_prob_vec_K(p_T_vec, "p_T_vec", K = length(p_C_vec))
  
  vartheta_C <- .prob_to_cum_K(p_C_vec)
  vartheta_T <- .prob_to_cum_K(p_T_vec)
  
  theta_vec <- logit(vartheta_T) - logit(vartheta_C)
  
  if (length(theta_vec) > 1 && sd(theta_vec) > tol) {
    warning("Treatment effect differs across cutpoints; PO assumption may not hold.")
  }
  
  mean(theta_vec)
}

pT_from_pC_pi_K <- function(p_C_vec, pi_vec) {
  .check_prob_vec_K(p_C_vec, "p_C_vec")
  .check_prob_vec_K(pi_vec, "pi_vec", K = length(p_C_vec))
  
  p_T_vec <- 2 * pi_vec - p_C_vec
  .check_prob_vec_K(p_T_vec, "Derived p_T_vec", K = length(p_C_vec))
  
  p_T_vec
}

pC_from_pT_pi_K <- function(p_T_vec, pi_vec) {
  .check_prob_vec_K(p_T_vec, "p_T_vec")
  .check_prob_vec_K(pi_vec, "pi_vec", K = length(p_T_vec))
  
  p_C_vec <- 2 * pi_vec - p_T_vec
  .check_prob_vec_K(p_C_vec, "Derived p_C_vec", K = length(p_T_vec))
  
  p_C_vec
}

pC_from_pi_theta_K <- function(pi_vec, theta, tol = 1e-10) {
  .check_prob_vec_K(pi_vec, "pi_vec")
  
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("theta must be a single finite numeric value.")
  }
  
  v_overall <- .prob_to_cum_K(pi_vec)
  v_C <- numeric(length(v_overall))
  
  for (k in seq_along(v_overall)) {
    target <- v_overall[k]
    
    f <- function(u) {
      0.5 * (u + expit(logit(u) + theta)) - target
    }
    
    lower <- 1e-10
    upper <- 1 - 1e-10
    
    f_low <- f(lower)
    f_up <- f(upper)
    
    if (f_low * f_up > 0) {
      stop(sprintf(
        "Could not solve for control cumulative probability at threshold %d. Inputs may be incompatible.",
        k
      ))
    }
    
    v_C[k] <- uniroot(f, interval = c(lower, upper), tol = tol)$root
  }
  
  if (any(diff(v_C) <= 0)) {
    warning("Derived control cumulative probabilities are not strictly increasing; inputs may be near boundary.")
  }
  
  p_C_vec <- .cum_to_prob_K(v_C)
  .check_prob_vec_K(p_C_vec, "Derived p_C_vec", K = length(pi_vec))
  
  p_C_vec
}

pT_from_pi_theta_K <- function(pi_vec, theta) {
  p_C_vec <- pC_from_pi_theta_K(pi_vec, theta)
  pi_T_from_PO_K(p_C_vec, theta)
}

resolve_ordinal_inputs <- function(p_C_vec = NULL,
                                   p_T_vec = NULL,
                                   pi_vec = NULL,
                                   treatment_effect = NULL,
                                   tol = 1e-4) {
  
  supplied <- c(
    p_C_vec = !is.null(p_C_vec),
    p_T_vec = !is.null(p_T_vec),
    pi_vec = !is.null(pi_vec),
    treatment_effect = !is.null(treatment_effect)
  )
  
  n_supplied <- sum(supplied)
  
  if (n_supplied != 2) {
    stop("Exactly two of p_C_vec, p_T_vec, pi_vec, treatment_effect must be provided.")
  }
  
  prob_lengths <- c(
    if (!is.null(p_C_vec)) length(p_C_vec) else NA_integer_,
    if (!is.null(p_T_vec)) length(p_T_vec) else NA_integer_,
    if (!is.null(pi_vec)) length(pi_vec) else NA_integer_
  )
  
  prob_lengths <- prob_lengths[!is.na(prob_lengths)]
  
  if (length(prob_lengths) > 0 && length(unique(prob_lengths)) > 1) {
    stop("All supplied probability vectors must have the same length.")
  }
  
  K <- if (length(prob_lengths) > 0) prob_lengths[1] else NULL
  
  .check_prob_vec_K(p_C_vec, "p_C_vec", K = K)
  .check_prob_vec_K(p_T_vec, "p_T_vec", K = K)
  .check_prob_vec_K(pi_vec, "pi_vec", K = K)
  
  if (!is.null(treatment_effect)) {
    if (!is.numeric(treatment_effect) || length(treatment_effect) != 1 || !is.finite(treatment_effect)) {
      stop("treatment_effect must be a single finite numeric value.")
    }
  }
  
  if (!is.null(p_C_vec) && !is.null(p_T_vec)) {
    treatment_effect <- theta_from_probs_K(p_C_vec, p_T_vec, tol = tol)
    pi_vec <- (p_C_vec + p_T_vec) / 2
  }
  
  if (!is.null(p_C_vec) && !is.null(treatment_effect)) {
    p_T_vec <- pi_T_from_PO_K(p_C_vec, treatment_effect)
    pi_vec <- (p_C_vec + p_T_vec) / 2
  }
  
  if (!is.null(p_T_vec) && !is.null(treatment_effect)) {
    p_C_vec <- pi_C_from_PO_K(p_T_vec, treatment_effect)
    pi_vec <- (p_C_vec + p_T_vec) / 2
  }
  
  if (!is.null(p_C_vec) && !is.null(pi_vec)) {
    p_T_vec <- pT_from_pC_pi_K(p_C_vec, pi_vec)
    treatment_effect <- theta_from_probs_K(p_C_vec, p_T_vec, tol = tol)
  }
  
  if (!is.null(p_T_vec) && !is.null(pi_vec)) {
    p_C_vec <- pC_from_pT_pi_K(p_T_vec, pi_vec)
    treatment_effect <- theta_from_probs_K(p_C_vec, p_T_vec, tol = tol)
  }
  
  if (!is.null(pi_vec) && !is.null(treatment_effect)) {
    p_C_vec <- pC_from_pi_theta_K(pi_vec, treatment_effect)
    p_T_vec <- pi_T_from_PO_K(p_C_vec, treatment_effect)
  }
  
  .check_prob_vec_K(p_C_vec, "Resolved p_C_vec")
  .check_prob_vec_K(p_T_vec, "Resolved p_T_vec", K = length(p_C_vec))
  .check_prob_vec_K(pi_vec, "Resolved pi_vec", K = length(p_C_vec))
  
  pi_check <- (p_C_vec + p_T_vec) / 2
  theta_check <- theta_from_probs_K(p_C_vec, p_T_vec, tol = tol)
  
  if (max(abs(pi_vec - pi_check)) > 1e-6) {
    warning("Resolved pi_vec is slightly inconsistent with (p_C_vec + p_T_vec)/2; using implied average.")
    pi_vec <- pi_check
  }
  
  if (abs(treatment_effect - theta_check) > 1e-3) {
    warning(sprintf(
      "Resolved treatment_effect (%.6f) differs slightly from value implied by p_C_vec and p_T_vec (%.6f). Using implied value.",
      treatment_effect, theta_check
    ))
    treatment_effect <- theta_check
  }
  
  list(
    p_C_vec = p_C_vec,
    p_T_vec = p_T_vec,
    pi_vec = pi_vec,
    treatment_effect = treatment_effect,
    K = length(pi_vec),
    specified = names(supplied)[supplied]
  )
}