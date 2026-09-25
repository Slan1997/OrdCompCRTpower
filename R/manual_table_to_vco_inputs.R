manual_table_to_vco_inputs <- function(tbl) {
  tbl <- tbl[order(tbl$Tree), ]
  
  method_list <- split(tbl$Copula, tbl$Tree)
  
  param_list <- lapply(split(tbl, tbl$Tree), function(df_tree) {
    mapply(
      function(fam, tau) {
        fam_id <- name_to_family(fam)
        if (is.na(fam_id) || fam_id == 0 || is.na(tau) || tau <= 0) {
          return(0)
        } else {
          return(BiCopTau2Par(family = fam_id, tau = tau))
        }
      },
      fam = df_tree$Copula,
      tau = df_tree$Tau,
      SIMPLIFY = TRUE
    )
  })
  
  list(
    method_list = method_list,
    param_list  = param_list
  )
}