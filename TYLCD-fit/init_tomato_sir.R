init_sir <- function(){
  list(
    theta = c(
      beta_p = 0.01,
      r_1 = 0.01, 
      r_2 = 0.01,
      b = 0.75,
      beta_v = 0.003,
      gamma = 0.15,
      gamma_f = 0.15,
      theta_1 = 0.4,
      mu = 1.0,
      C = 50
    ),
    Lp0 = 1,
    Ip0 = 0.0001,
    Sv0 = 25, 
    Iv0 = 800,
    Dv0 = 985
  )
}