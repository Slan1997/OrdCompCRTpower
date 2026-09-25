wh_power = function(theta_R,n, pr,A=1,alpha=0.05){
  sum_cube = sum(pr^3) # pr is the average proportions for all categories
  # A is the ratio between control and treatment (randomization: B is treatment, C is placebo)
  mu_beta1 = theta_R*sqrt( A*n^3*(1-sum_cube) /3/(n+1)^2/(A+1)^2 ) - qnorm(1-alpha/2)
  pnorm(mu_beta1)
}
