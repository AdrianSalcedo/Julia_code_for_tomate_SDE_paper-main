function Compute_auxiliar_constants_ci(par,endemic_fixed_point)
    # Input:
    #   par: vector with parameters
    #   endemic_fixed_point: endemic_fixed_point vector
    # Output:
    #   auxiliar_constants_c_i: vector with auxiliar constants C, c_i

    beta_p = par.beta_p[1]

    N_p = par.N_p[1]
    N_v = par.N_v[1]
    sigma_L = par.sigma_L[1]
    sigma_I = par.sigma_I[1]
    sigma_v = par.sigma_v[1]
    L_p_aster = endemic_fixed_point.L_p_aster[1]

    c_1 = (sigma_L) / (N_p ^ 2)
    c_2 = 1 / (N_p ^ 2)
    c_3 = 1 / (N_v ^ 2)
    C = 2 * beta_p * N_p * c_1 +
        (1 / 2) * c_1 * sigma_L ^ 2 * L_p_aster +
            (1 / 2) *
                (c_2 * sigma_I ^ 2 * N_p ^ 2 + c_3 * sigma_v ^ 2 * N_v ^ 2)
    auxiliar_constants_c_i = DataFrame(C = C, c_1 = c_1, c_2 = c_2, c_3 = c_3)

return auxiliar_constants_c_i
end
