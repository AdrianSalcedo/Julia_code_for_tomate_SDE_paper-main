function Compute_stochastic_R0(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   stochastic R0: stochastic reproductive number

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    r_2 = par.r_2[1]
    b = par.b[1]
    beta_v = par.beta_v[1]
    gamma = par.gamma[1]
    gamma_f = par.gamma_f[1]
    theta = par.theta[1]
    sigma_L = par.sigma_L[1]
    sigma_I = par.sigma_I[1]
    sigma_v = par.sigma_v[1]

    Rs0 = (beta_p * beta_v * b) / ((gamma + gamma_f) * (b + r_1) * r_2) -
        (1 / 2) * ((sigma_L + sigma_I) ^ 2 - sigma_v ^ 2 /
        (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f)))

    return Rs0
end
