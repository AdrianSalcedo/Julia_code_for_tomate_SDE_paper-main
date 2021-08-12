function Compute_auxiliar_constants_rho_i(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   auxiliar_constants_rho_i: vector with auxiliar constants rho_i

    r_2 = par.r_2[1]
    b = par.b[1]

    rho_1 = rand(Uniform((1/2)+b/(4r_2),1+b/(b+2*r_2)))
    rho_2 = rand(Uniform((1/2),1))
    
    auxiliar_constants_rho_i = DataFrame(rho_1 = rho_1, rho_2 = rho_2)
    return auxiliar_constants_rho_i
end
