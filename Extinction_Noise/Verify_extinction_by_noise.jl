function Verify_extinction_by_noise(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   No noise condition: verify that satisfy no noise extinction conditions

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    r_2 = par.r_2[1]
    b = par.b[1]
    beta_v = par.beta_v[1]
    theta = par.theta[1]
    mu = par.mu[1]
    gamma = par.gamma[1]
    gamma_f = par.gamma_f[1]
    sigma_L = par.sigma_L[1]
    sigma_I = par.sigma_I[1]
    sigma_v = par.sigma_v[1]

    cond_1 = beta_p ^ 2 / (2 * sigma_L ^ 2) +
        r_2 ^ 2 / (2 * sigma_I ^ 2) + 2 * beta_p - r_1
    #sigma_I ^ 2 - (sigma_L ^ 2 * r_2 ^ 2) /(2 * sigma_L ^ 2* r_1 - 4 * beta_p * sigma_L ^ 2 - beta_p ^ 2)
    cond_2 = sigma_v ^ 2 - beta_v ^ 2 / (2 * ((gamma + gamma_f) - theta * mu - beta_v))
    test1 = sign(cond_1) == -1.0
    test2 = sign(cond_2) == -1.0
    ans = test1 && test2

    return ans
end
