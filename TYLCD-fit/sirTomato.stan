functions {
  real[] SIR(real t,  // time
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {

    real dy_dt[5];
    real Sp = y[1];
    real Lp = y[2];
    real Ip = y[3];
    real Sv = y[4];
    real Iv = y[5];
    // real Dv = y[6];
    // real Np = y[7];
    // real Nv = y[8];
    // real Ypc = y[9];
    real N_v;
    real N_p;
//
    real beta_p = theta[1];
    real r_1 = theta[2];
    real r_2 = theta[3];
    real b = theta[4];
    real beta_v = theta[5];
    real gamma = theta[6];
    real gamma_f = theta[7];
    real theta_1 = theta[8];
    real mu = theta[9];
//
    N_v = Sv + Iv;
    N_p = Sp + Lp + Ip;
    
    dy_dt[1] = -beta_p * Sp * Iv/N_v + r_1 * Lp + r_2 * Ip;
    dy_dt[2] = beta_p * Sp * Iv/N_v - (b + r_1) * Lp;
    dy_dt[3] = b * Lp - r_2 * Ip;
    dy_dt[4] = -beta_v * Sv * Ip/N_p - (gamma + gamma_f) * Sv + (1 - theta_1) * mu;
    dy_dt[5] = beta_v * Sv * Ip/N_p - (gamma + gamma_f) * Iv + theta_1 * mu;
    // dy_dt[6] = (gamma + gamma_f) * (Sv + Iv) - mu;
    // dy_dt[7] = -beta_p * Sp * Iv/N_v + r_1 * Lp + r_2 * Ip +
    // beta_p * Sp * Iv/N_v - (b + r_1) * Lp + b * Lp - r_2 * Ip;
    // dy_dt[8] = -beta_v * Sv * Ip/N_p - (gamma + gamma_f) * Sv + (1 - theta_1) * mu +
    // beta_v * Sv * Ip/N_p - (gamma + gamma_f) * Iv + theta_1 * mu +
    // (gamma + gamma_f) * (Sv + Iv) - mu;
    // dy_dt[9] = b * Lp;
    return dy_dt;
  }

}
data {
    int<lower = 1> n_obs;       // number of days observed
    int<lower = 1> n_theta;     // number of model parameters
    int<lower = 1> n_difeq;     // number of differential equations
    int<lower = 1> n_pop;       // population
    int y[n_obs];               // data, total number of infected
    real t0;                    // initial time point (zero)
    real ts[n_obs];             // time points observed
}

transformed data {
  real x_r[0];
  int x_i[0];
  //real<lower = 0.0001, upper = 0.1> beta_p;
  // real<lower = 0.0083, upper = 0.013> r_1;
  // real<lower = 0.0083, upper = 0.013> r_2;
  //real<lower = 0.05, upper = 0.1> b;
  //real<lower = 0.0001, upper = 0.2> beta_v;
  // real<lower = 0.03, upper = 1> gamma;
  // real<lower = 0.0001, upper = 0.5> gamma_f;
  // real<lower = 0.0001, upper = 1> theta_1;
  // real<lower = 0.0001, upper = 1> mu;
  //
  //beta_p = 0.1;
  // r_1 = 0.01;
  // r_2 = 0.01;
  //b = 0.075;
  //beta_v = 0.01;
  // gamma = 0.06;
  // gamma_f = 0.06;
  // theta_1 = 0.4;
  // mu =1.0;

}
parameters {
    //real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real<lower = 0, upper = 10> beta_p;
    real<lower = 0.0083, upper = 10> r_1;
    real<lower = 0.0083, upper = 10> r_2;
    real<lower = 0.05, upper = 10> b;
    real<lower = 0, upper = 10> beta_v;
    real<lower = 0.0, upper = 10> gamma;
    real<lower = 0, upper = 10> gamma_f;
    real<lower = 0, upper = 10> theta_1;
    real<lower = 0, upper = 10> mu;
    //real <lower = 0, upper =1000> Lp0;
    real <lower = 1, upper = 20> Lp0;
    real <lower = 1, upper= 10> Ip0;
    // real <lower = 0, upper= 50> Sv0;
    //real <lower = 0, upper= 2> Iv0;
    //real <lower = 0> N_v;
    real<lower=0> phi_inv;
}

transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];
    real incidence[n_obs - 1];
    //initial conditions
    // real Lp0 = 1;
    // real Ip0 = 1;
    real Sv0 = 100;
    real Iv0 = 10;
    // real N_v0 = 5.0;
    // real Ypc0 = 5.0;
    
    real theta[n_theta];
    real<lower=0> phi = 1. / phi_inv;

    theta[1] = beta_p;
    theta[2] = r_1;
    theta[3] = r_2;
    theta[4] = b;
    theta[5] = beta_v;
    theta[6] = gamma;
    theta[7] = gamma_f;
    theta[8] = theta_1;
    theta[9] = mu;

    y_init[1] = n_pop - (Lp0 + Ip0);
    y_init[2] = Lp0;
    y_init[3] = Ip0;
    y_init[4] = Sv0;
    y_init[5] = Iv0;
    // y_init[6] = N_v0 - (Sv0 + Iv0);
    // y_init[7] = n_pop;
    // y_init[8] = N_v0;
    // y_init[9] = Ypc0;

    y_hat = integrate_ode_adams(SIR, y_init, t0, ts, theta, x_r, x_i);
    for(i in 1:n_obs-1){
    incidence[i] =  -(y_hat[i+1, 2] - y_hat[i, 2] + y_hat[i+1, 1] - y_hat[i, 1]);
    //print(incidence)
    }
}

model {
    real lambda[n_obs];      //poisson parameter
    //priors
    beta_p ~ uniform(0,10);//normal(0.01,0.001);//uniform(0,1);
    r_1 ~ uniform(0,10);//gamma(0.1,0.1);//
    r_2 ~ uniform(0,10);//gamma(0.1,0.1);//
    b ~ uniform(0,10);//gamma(0.225,0.333);//
    beta_v ~ uniform(0,10);//normal(0.003,0.001);
    gamma ~ uniform(0,10);//gamma(0.6,0.1);//
    gamma_f ~ uniform(0,10);//gamma(0.36, 0.167);
    theta_1 ~ uniform(0, 10);
    mu ~ uniform(0, 10);
    // 
    Lp0 ~ uniform(1,30);//gamma(12.5, 0.4);
    Ip0 ~ uniform(1,30);//gamma(12.5, 0.4);
    //Sv0 ~ uniform(0, 50);
    //Iv0 ~ gamma(40, 0.05);
    phi_inv ~ exponential(5);
       //likelihood
     for (i in 1:n_obs-1){
      lambda[i] = incidence[i];
    }
    //y ~ poisson(lambda);
    y[1:n_obs] ~ neg_binomial_2(y_hat[,3],phi);
}

generated quantities {
    real R_0;      // Basic reproduction number
    //real pred_cases[n_obs];
    R_0 = beta_p * beta_v * b /((b + r_1) * (gamma + gamma_f)* r_2);
    
    //pred_cases = neg_binomial_2_rng(col(to_matrix(y_hat),9), phi);
}