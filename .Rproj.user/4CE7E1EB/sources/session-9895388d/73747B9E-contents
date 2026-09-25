# rho_binary2 <- function(p1, p2, rho11, rho12, rho22) {
#   p3 <- 1 - p1 - p2
#   
#   bad_prob <- p1 < 0 | p2 < 0 | p3 < 0
#   bad_var  <- p3 == 0 | (p1 + p2) == 0
#   
#   out <- (
#     rho11 * p1 * (1 - p1) +
#       2 * rho12 * sqrt(p1 * (1 - p1) * p2 * (1 - p2)) +
#       rho22 * p2 * (1 - p2)
#   ) / ((p1 + p2) * p3)
#   
#   out[bad_prob | bad_var] <- NA_real_
#   out
# }
# 
# dt = read_csv("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/rhomat_from_superpop.csv")
# dt <- dt %>%
#   mutate(
#     rhobin2 = rho_binary2(
#       p1 = mean.pi_0,
#       p2 = mean.pi_1,
#       rho11 = rho11,
#       rho12 = rho12,
#       rho22 = rho22
#     )
#   )
# dt %>% select(bin_rho01, rhobin2) %>% as.data.frame
# rho_binary2(p1=0.33, p2=0.33, rho11=0.005057269, rho12=-0.0003543616, rho22=0.0008008906)

rho_binary2_from_p <- function(p, rho11, rho12, rho22) {
  if (length(p) != 3) stop("p must be a length-3 vector: c(p1, p2, p3)")
  if (abs(sum(p) - 1) > 1e-8) stop("Probabilities must sum to 1.")
  
  p1 <- p[1]
  p2 <- p[2]
  p3 <- p[3]
  
  if (any(p < 0)) stop("Probabilities must be nonnegative.")
  if ((p1 + p2) == 0 || p3 == 0) {
    stop("Binary2 variance is zero, so correlation is undefined.")
  }
  
  (
    rho11 * p1 * (1 - p1) +
      2 * rho12 * sqrt(p1 * (1 - p1) * p2 * (1 - p2)) +
      rho22 * p2 * (1 - p2)
  ) / ((p1 + p2) * p3)
}

# rho_binary2_from_p(p=c(0.33,0.33,0.34),rho11=0.005057269, rho12=-0.0003543616, rho22=0.0008008906)

# apply(as.data.frame(dt),1,\(x) rho_binary2(p1=x$pi_0, p2=x$pi_1, rho11=x$rho11 , rho12=x$rho12, rho22=x$rho22))
# dt %>% mutate(rhobin2 = rho_binary2(p1=pi_0, p2=pi_1, rho11=rho11 , rho12=rho12, rho22=rho22))
