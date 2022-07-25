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
  //real<lower = 0.0001, upper = 0.1> beta_p;
  //real<lower = 0, upper = 1> r_1;
  // real<lower = 0.0083, upper = 0.013> r_2;
  //real<lower = 0.05, upper = 0.1> b;
  //real<lower = 0.0001, upper = 0.2> beta_v;
  // real<lower = 0.03, upper = 1> gamma;
  // real<lower = 0.0001, upper = 0.5> gamma_f;
  // real<lower = 0.0001, upper = 1> theta_1;
  // real<lower = 0.0001, upper = 1> mu;

}
parameters {
    //real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real<lower = 0, upper = 0.1> beta_p;
    real<lower = 0.0083, upper = 0.013> r_1;
    real<lower = 0.0083, upper = 0.013> r_2;
     real<lower = 0.05, upper = 0.1> b;
     real<lower = 0, upper = 0.2> beta_v;
    //real<lower = 0.03, upper = 1> gamma;
    // real<lower = 0.03, upper = 1> gamma_f;
    // real<lower = 0.0, upper = 1.0> theta_1;
    //real<lower = 0, upper = 1> mu;
    real<lower = 0, upper = 150> Lp0;
    // real<lower = 0.0001, upper = 10> Ip0;
    // real<lower = 0, upper = 3000> Iv0;
    real<lower = 0, upper = 300> Yp0;
    //real <lower = 0> N_v;
    real<lower=0> phi_inv;
    //real<lower = 40, upper = 110> C;
    real<lower=0, upper=1> p_reported; // proportion of infected (symptomatic) people reported
}

transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];
    real incidence[n_obs-1];
    //initial conditions
    //real Lp0 = 20;
    real Ip0 = 5;
    real Sv0 = 1000;
    real Iv0 = 40000; //0-5 per plant
    real Dv0 = 9000;
    real Nv = 50 * n_pop;
    //real Yp0 = 1.0;
    
    //real beta_p = 0.1;
    //real r_1 = 0.01;
    //real r_2 = 0.01;
    //real b = 0.065;
    //real beta_v = 0.01;
    real gamma = 0.03;
    real gamma_f = 0.03;
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
    y_init[7] = Yp0;

    y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);

}

model {
    real lambda[n_obs-1];      //poisson parameter
    //priors 
    //Gamma ~(alpha,beta) mean alpha/beta variance alpha / beta²
    //Invense Gamma ~(alpha,beta) mean beta/(alpha-1), alpha >1 
    //variance beta² / (alpha- 1)²(alpha -2)
    
    beta_p ~ normal(0.1,0.1);
    r_1 ~ normal(0.01,0.01);//inv_gamma(2.02,0.0102);//inv_gamma(3.125,0.159375);//inv_gamma(2.01,0.0101);
    r_2 ~ normal(0.01,0.01);//inv_gamma(2.01,0.0101);//inv_gamma(2.01,0.01002);//inv_gamma(2.005,0.01005);
     b ~ normal(0.075,0.01);//inv_gamma(27,1.3);//inv_gamma(58.25,4.29375);//inv_gamma(2.05625,0.0792188);
    beta_v ~ normal(0.01,0.001);
    //C ~ uniform(50,100);
    //gamma ~ inv_gamma(2.018,0.03054);
    //gamma_f ~ inv_gamma(2.018,0.03054);//inv_gamma(12,3);
    //theta_1 ~ normal(0.4, 0.0001);
    //mu ~ normal(1.0, 0.0001);

    // 
    Lp0 ~ normal(120,20);
    // Ip0 ~ normal(50,10);
    // Iv0 ~ normal(5000, 500);
    Yp0 ~ normal(70,20);//normal(100,20);
    phi_inv ~ exponential(1.9);
  
       //likelihood
     for (i in 1:n_obs-1){
      lambda[i] = y_hat[i,7];
    }
    //y[1:(n_obs-1)] ~ poisson(lambda);
    y[1:(n_obs-1)] ~ neg_binomial_2(lambda,phi);

}

generated quantities {
    real R_0;      // Basic reproduction number

    R_0 = (beta_p * beta_v * b /((b + r_1) * (gamma + gamma_f)* r_2))^(0.5);

}