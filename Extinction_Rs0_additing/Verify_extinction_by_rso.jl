function Verify_extinction_by_rso(par)

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    r_2 = par.r_2[1]
    b = par.b[1]
    beta_v = par.beta_v[1]
    mu = par.mu[1]

    cond1 = r_2 - r_1
    cond2 = beta_p * beta_v - 2 * mu * r_2

    test1 = sign(cond1) == 1.0
    test2 = sign(cond2) == 1.0
    ans = test1 && test2
    return ans
end
