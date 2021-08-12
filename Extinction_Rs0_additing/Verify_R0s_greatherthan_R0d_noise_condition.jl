function Verify_R0s_greatherthan_R0d_noise_condition(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   No noise condition: verify that satisfy no noise extinction conditions

    beta_v = par.beta_v[1]
    theta = par.theta[1]
    gamma = par.gamma[1]
    gamma_f = par.gamma_f[1]
    sigma_L = par.sigma_L[1]
    sigma_I = par.sigma_I[1]
    sigma_v = par.sigma_v[1]

    cond = (sigma_L + sigma_I) ^ 2 -
        sigma_v ^ 2 / (2 * (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f)))
    test = sign(cond) == -1.0

    return test
end
