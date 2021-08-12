function Verify_rs_persistence_condition(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   stochastic R0: verify Rs0 persistence condition

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    r_2 = par.r_2[1]
    beta_v = par.beta_v[1]
    gamma = par.gamma[1]
    gamma_f = par.gamma_f[1]
    theta = par.theta[1]
    sigma_L = par.sigma_L[1]
    sigma_I = par.sigma_I[1]
    sigma_v = par.sigma_v[1]

    Rs0 = Compute_stochastic_R0(par)
    test = Rs0 > 1
    return test
end
