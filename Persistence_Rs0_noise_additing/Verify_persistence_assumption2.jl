function Verify_persistence_assumption2(
    auxiliar_constants_c_i,auxiliar_constants_rho_i,par,endemic_fixed_point
    )
    # Input:
    #   auxiliar_constants_c_i:
    #   auxiliar_constants_rho_i:
    #   par: vector with parameters
    #   endemic_fixed_point: endemic_fixed_point vector
    # Output:
    #   cond: verify minimum condition

    r_2 = par.r_2[1]
    b = par.b[1]
    gamma = par.gamma[1]
    gamma_f = par.gamma_f[1]
    S_p_aster = endemic_fixed_point.S_p_aster[1]
    L_p_aster = endemic_fixed_point.L_p_aster[1]
    L_p_aster = endemic_fixed_point.I_p_aster[1]
    S_v_aster = endemic_fixed_point.S_v_aster[1]
    I_v_aster = endemic_fixed_point.I_v_aster[1]
    C = auxiliar_constants_c_i.C[1]
    c_1 = auxiliar_constants_c_i.c_1[1]
    c_2 = auxiliar_constants_c_i.c_2[1]
    c_3 = auxiliar_constants_c_i.c_3[1]
    rho_1 = auxiliar_constants_rho_i.rho_1[1]
    rho_2 = auxiliar_constants_rho_i.rho_2[1]

    a_1 = r_2 - (b + 2 * r_2) / (4 * rho_1)
    a_2 = b + r_2 - rho_1 * ((b + 2 * r_2) / 2)
    a_3 = (1 - 1 / (2 * rho_2))
    a_4 = (1 - rho_2)
    auxiliar_constants_a_i = DataFrame(a_1 = a_1, a_2 = a_2, a_3 = a_3,
     a_4 = a_4)
    alpha_1 = c_2 * a_1
    alpha_2 = c_2 * a_2
    alpha_3 = c_3 * a_3
    alpha_4 = c_3 * a_4

    LC = min(c_2 * a_1 * (S_p_aster) ^ 2, c_2 * a_2 * (L_p_aster) ^ 2,
        c_3 * (gamma + gamma_f) * a_3 * (S_v_aster) ^ 2,
        c_3 * (gamma + gamma_f) * a_4 * (I_v_aster) ^ 2)

    cond = LC >= C
    return cond, auxiliar_constants_a_i
end
