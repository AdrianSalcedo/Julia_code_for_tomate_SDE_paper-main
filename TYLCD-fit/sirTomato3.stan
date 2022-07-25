functions {
  real[] SIR(real t,  // time
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {

    real dy_dt[7];
    real Sp = y[1];
    real Lp = y[2];
    real Ip = y[3];
    real Sv = y[4];
    real Iv = y[5];
    real Dv = y[6];
    real Yp = y[7];
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
    real Np = Sp + Lp + Ip;
    real Nv = Sv + Iv;
    
    dy_dt[1] = -beta_p * Sp * Iv/ Nv + r_1 * Lp + r_2 * Ip;
    dy_dt[2] = beta_p * Sp * Iv /Nv - (b + r_1) * Lp;
    dy_dt[3] = b * Lp - r_2 * Ip;
    dy_dt[4] = -beta_v * Sv * Ip / Np - (gamma + gamma_f) * Sv + (1 - theta_1) * mu;
    dy_dt[5] = beta_v * Sv * Ip / Np - (gamma + gamma_f) * Iv + theta_1 * mu;
    dy_dt[6] = (gamma + gamma_f) * (Sv + Iv) - mu;
    dy_dt[7] = b * Lp;
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
  
}
parameters {
    //real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real<lower = 0.0001, upper = 0.1> beta_p;
    real<lower = 0.0001, upper = 0.1> r_1;
    real<lower = 0.0001, upper = 0.1> r_2;
    real<lower = 0.05, upper = 0.1> b;
    real<lower = 0.0001, upper = 0.2> beta_v;
    // real<lower = 0, upper = 1> gamma;
    // real<lower = 0, upper = 1> gamma_f;
    // real<lower = 0.0, upper = 1.0> theta_1;
    // real<lower = 0, upper = 0.5> mu;
    //real<lower = 0.000001, upper = 10> Lp0;
    //real<lower = 0.000001, upper = 10> Ip0;
    //real<lower = 0.000001, upper = 30000> Iv0;
    //real<lower = 0.0001, upper = 1200> Dv0;
    real<lower=0> phi_inv;
    //real<lower = 40, upper = 110> C;
    
}

transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];
    real incidence[n_obs-1];
    //initial conditions
     real Lp0 = 0.1;
     real Ip0 = 0.1;
    real Sv0 = 15;
     real Iv0 = 0; //0-5 per plant
     real Dv0 = 985;
    real Nv = 1 * n_pop;
    //real Yp0 = 1.0;
    
    //beta_p = 0.1;
    //r_1 = 0.01;
    // r_2 = 0.01;
    //b = 0.075;
    //beta_v = 0.01;
    real gamma = 0.15;
    real gamma_f = 0.15;
    real theta_1 = 0.4;
    real mu = 1.0;
    
    
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
    //theta[10] = C;
    
    y_init[1] = n_pop - (Lp0 + Ip0);
    y_init[2] = Lp0;
    y_init[3] = Ip0;
    y_init[4] = Nv - (Dv0 + Iv0);
    y_init[5] = Iv0;
    y_init[6] = Nv - (y_init[4] + Iv0);
    y_init[7] = 1;

    y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
    
}

model {
    real lambda[n_obs-1];      //poisson parameter
    //priors 
    //Gamma ~(alpha,beta) mean alpha/beta variance alpha / beta²
    //Invense Gamma ~(alpha,beta) mean beta/(alpha-1), alpha >1 
    //variance beta² / (alpha- 1)²(alpha -2)
    
    beta_p ~ normal(0.01,0.005);
    r_1 ~ inv_gamma(2.005,0.01005);
    r_2 ~ inv_gamma(2.005,0.01005);
    b ~ inv_gamma(7.625,0.496875);
    beta_v ~ normal(0.001,0.001);
    //C ~ uniform(50,100);
    //gamma ~ inv_gamma(2.23,0.073);
    //gamma_f ~ inv_gamma(2.23,0.073);
    //theta_1 ~ normal(0.4, 0.1);
    //mu ~ normal(0.39, 0.005);

    // 
    //Lp0 ~ uniform(0,10);
    //Ip0 ~ uniform(0,10);
    //Iv0 ~ normal(25, 25);
    //Dv0 ~ normal(1183, 1);
    phi_inv ~ exponential(1.9);
       //likelihood
     for (i in 1:n_obs-1){
      lambda[i] = y_hat[i,7];
    }
    y[1:(n_obs-1)] ~ neg_binomial_2(lambda,phi);

}

generated quantities {
    real R_0;      // Basic reproduction number
    
    R_0 = (beta_p * beta_v * b /((b + r_1) * (gamma + gamma_f)* r_2))^(0.5);
    
}