parse_or_grid <- function(mode = c("sequence", "custom"),
                          or_min = NULL,
                          or_max = NULL,
                          or_step = NULL,
                          or_custom = NULL) {
  mode <- match.arg(mode)
  
  if (mode == "sequence") {
    if (!is.numeric(or_min) || length(or_min) != 1 || !is.finite(or_min) || or_min <= 0) {
      stop("or_min must be a single positive finite number.")
    }
    if (!is.numeric(or_max) || length(or_max) != 1 || !is.finite(or_max) || or_max <= 0) {
      stop("or_max must be a single positive finite number.")
    }
    if (!is.numeric(or_step) || length(or_step) != 1 || !is.finite(or_step) || or_step <= 0) {
      stop("or_step must be a single positive finite number.")
    }
    if (or_min > or_max) {
      stop("or_min must be less than or equal to or_max.")
    }
    
    grid <- seq(or_min, or_max, by = or_step)
    grid <- unique(round(grid, 10))
    return(grid)
  }
  
  if (is.null(or_custom) || !nzchar(trimws(or_custom))) {
    stop("Please provide at least one custom OR entry.")
  }
  
  parts <- unlist(strsplit(or_custom, ","))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  
  vals <- suppressWarnings(as.numeric(parts))
  
  if (length(vals) == 0 || any(is.na(vals)) || any(vals <= 0)) {
    stop("Custom OR entries must be comma-separated positive numeric values.")
  }
  
  vals <- unique(round(vals, 10))
  vals <- sort(vals)
  vals
}

