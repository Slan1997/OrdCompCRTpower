# -----------------------------
# Helper: blank placeholder plot
# -----------------------------
plot_blank_message <- function(msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 6, color = "#666666") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

# -----------------------------
# COVID structure plot
# -----------------------------
plot_covid_vine_structure <- function(
    fam1, tau1,
    fam2, tau2,
    fam3a, tau3a,
    fam3b, tau3b,
    fam3c, tau3c
) {
  edge_label_from_tau <- function(fam, tau, conditional_text = NULL) {
    fam_id <- name_to_family(fam)
    
    par_val <- if (is.na(fam_id) || fam_id == 0 || is.na(tau) || tau <= 0) {
      0
    } else {
      BiCopTau2Par(family = fam_id, tau = tau)
    }
    
    par_txt <- format(round(par_val, 3), nsmall = 3, trim = TRUE)
    
    if (is.null(conditional_text)) {
      paste0(fam, "\npar = ", par_txt)
    } else {
      paste0(fam, "\npar = ", par_txt, "\n", conditional_text)
    }
  }
  
  # -------------------------
  # Subgroup 1 and 2
  # -------------------------
  nodes12 <- rbind(
    data.frame(
      subgroup = "Subgroup 1: Hospitalization — Vent/ECMO",
      node = c("Hospital", "Vent/ECMO"),
      x = c(1, 3),
      y = c(1, 1),
      stringsAsFactors = FALSE
    ),
    data.frame(
      subgroup = "Subgroup 2: Hospitalization — Supp. Oxygen",
      node = c("Hospital", "Oxygen"),
      x = c(1, 3),
      y = c(1, 1),
      stringsAsFactors = FALSE
    )
  )
  
  edges12 <- rbind(
    data.frame(
      subgroup = "Subgroup 1: Hospitalization — Vent/ECMO",
      x = 1, y = 1, xend = 3, yend = 1,
      label = edge_label_from_tau(fam1, tau1),
      lx = 2, ly = 1.12,
      stringsAsFactors = FALSE
    ),
    data.frame(
      subgroup = "Subgroup 2: Hospitalization — Supp. Oxygen",
      x = 1, y = 1, xend = 3, yend = 1,
      label = edge_label_from_tau(fam2, tau2),
      lx = 2, ly = 1.12,
      stringsAsFactors = FALSE
    )
  )
  
  # -------------------------
  # Subgroup 3: Tree 1 on top, Tree 2 below
  # -------------------------
  # Tree 1 top
  sg3_nodes_t1 <- data.frame(
    subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
    tree = "Tree 1",
    node = c("Hospital", "Symptoms", "Activity"),
    x = c(1, 2, 3),
    y = c(2, 2, 2),
    stringsAsFactors = FALSE
  )
  
  # Tree 2 below
  sg3_nodes_t2 <- data.frame(
    subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
    tree = "Tree 2",
    node = c("Hospital|Symptoms", "Activity|Symptoms"),
    x = c(1.5, 2.5),
    y = c(1, 1),
    stringsAsFactors = FALSE
  )
  
  sg3_nodes <- rbind(sg3_nodes_t1, sg3_nodes_t2)
  
  sg3_edges_t1 <- rbind(
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = 1, y = 2, xend = 2, yend = 2,
      label = edge_label_from_tau(fam3a, tau3a),
      lx = 1.5, ly = 2.14,
      stringsAsFactors = FALSE
    ),
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = 2, y = 2, xend = 3, yend = 2,
      label = edge_label_from_tau(fam3b, tau3b),
      lx = 2.5, ly = 2.14,
      stringsAsFactors = FALSE
    )
  )
  
  sg3_edges_t2 <- data.frame(
    subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
    x = 1.5, y = 1, xend = 2.5, yend = 1,
    label = edge_label_from_tau(fam3c, tau3c, conditional_text = "cond. on Symptoms"),
    lx = 2.0, ly = 1.22,
    stringsAsFactors = FALSE
  )
  
  # dotted guides from Tree 1 down to Tree 2
  sg3_guides <- rbind(
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = 1, y = 1.92, xend = 1.5, yend = 1.08
    ),
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = 2, y = 1.92, xend = 1.5, yend = 1.08
    ),
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = 2, y = 1.92, xend = 2.5, yend = 1.08
    ),
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = 3, y = 1.92, xend = 2.5, yend = 1.08
    )
  )
  
  sg3_tree_labels <- data.frame(
    subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
    x = c(0.72, 0.72),
    y = c(2, 1),
    label = c("Tree 1", "Tree 2"),
    stringsAsFactors = FALSE
  )
  
  # blank rows to control facet heights
  facet_space <- rbind(
    data.frame(
      subgroup = "Subgroup 1: Hospitalization — Vent/ECMO",
      x = c(0.5, 3.5),
      y = c(0.75, 1.35)
    ),
    data.frame(
      subgroup = "Subgroup 2: Hospitalization — Supp. Oxygen",
      x = c(0.5, 3.5),
      y = c(0.75, 1.35)
    ),
    data.frame(
      subgroup = "Subgroup 3: Hospital — Symptoms — Activity",
      x = c(0.5, 3.5),
      y = c(0.65, 2.45)
    )
  )
  
  ggplot() +
    geom_blank(
      data = facet_space,
      aes(x = x, y = y)
    ) +
    
    # subgroup 1 and 2
    geom_segment(
      data = edges12,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.8,
      color = "#4F6D7A"
    ) +
    geom_label(
      data = edges12,
      aes(x = lx, y = ly, label = label),
      size = 3.5,
      fill = "white",
      label.size = 0.15
    ) +
    geom_point(
      data = nodes12,
      aes(x = x, y = y),
      size = 5,
      color = "#2C3E50"
    ) +
    geom_label(
      data = nodes12,
      aes(x = x, y = y - 0.12, label = node),
      size = 3.6,
      fill = "#F8FBFA",
      label.size = 0.15
    ) +
    
    # subgroup 3
    geom_segment(
      data = sg3_edges_t1,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.8,
      color = "#4F6D7A"
    ) +
    geom_segment(
      data = sg3_edges_t2,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.8,
      linetype = "dashed",
      color = "#B85C6B"
    ) +
    geom_segment(
      data = sg3_guides,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.35,
      linetype = "dotted",
      color = "#999999"
    ) +
    geom_label(
      data = sg3_edges_t1,
      aes(x = lx, y = ly, label = label),
      size = 3.4,
      fill = "white",
      label.size = 0.15
    ) +
    geom_label(
      data = sg3_edges_t2,
      aes(x = lx, y = ly, label = label),
      size = 3.3,
      fill = "white",
      label.size = 0.15
    ) +
    geom_point(
      data = sg3_nodes,
      aes(x = x, y = y),
      size = 5,
      color = "#2C3E50"
    ) +
    geom_label(
      data = sg3_nodes,
      aes(x = x, y = y - 0.14, label = node),
      size = 3.5,
      fill = "#F8FBFA",
      label.size = 0.15
    ) +
    geom_text(
      data = sg3_tree_labels,
      aes(x = x, y = y, label = label),
      hjust = 0,
      fontface = "bold",
      size = 4.0,
      color = "#555555"
    ) +
    
    facet_wrap(~ subgroup, ncol = 1, scales = "free_y") +
    coord_cartesian(xlim = c(0.5, 3.5), clip = "off") +
    theme_void(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(0.9, "lines"),
      plot.margin = margin(10, 20, 10, 30)
    )
}

# -----------------------------
# Helper: one bi-copula density grid
# -----------------------------
build_bicop_density_df <- function(panel, family_name, tau, n_grid = 60) {
  fam_id <- name_to_family(family_name)
  
  u <- seq(0.01, 0.99, length.out = n_grid)
  grid <- expand.grid(u1 = u, u2 = u)
  
  if (is.na(fam_id) || fam_id == 0 || is.na(tau) || tau <= 0) {
    grid$z <- 1
  } else {
    par_val <- BiCopTau2Par(family = fam_id, tau = tau)
    grid$z <- VineCopula::BiCopPDF(
      u1 = grid$u1,
      u2 = grid$u2,
      family = fam_id,
      par = par_val
    )
  }
  
  grid$panel <- panel
  grid$family <- family_name
  grid$tau <- tau
  grid
}

# -----------------------------
# COVID bi-copula shape plot
# -----------------------------
plot_covid_bicop_shapes <- function(
    fam1, tau1,
    fam2, tau2,
    fam3a, tau3a,
    fam3b, tau3b,
    fam3c, tau3c
) {
  df <- rbind(
    build_bicop_density_df("Hospital — Vent/ECMO", fam1, tau1),
    build_bicop_density_df("Hospital — Supp. Oxygen", fam2, tau2),
    build_bicop_density_df("Hospital — Symptoms", fam3a, tau3a),
    build_bicop_density_df("Symptoms — Activity", fam3b, tau3b),
    build_bicop_density_df("Hospital — Activity | Symptoms", fam3c, tau3c)
  )
  
  ggplot(df, aes(x = u1, y = u2, fill = z)) +
    geom_raster() +
    geom_contour(aes(z = z), color = "white", bins = 6, linewidth = 0.25) +
    facet_wrap(~ panel, ncol = 2) +
    coord_equal() +
    scale_fill_gradient(low = "#EAF4F1", high = "#4F9F8A") +
    labs(
      title = "Bi-copula Shapes for COVID Subgroups",
      x = "u1",
      y = "u2",
      fill = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}