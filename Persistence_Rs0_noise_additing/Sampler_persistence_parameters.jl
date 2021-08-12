using CSV
using IterableTables, DataFrames, DataTables
using Distributions

path1 = "D:/Julia_code_for_tomate_SDE_paper-main/"
path2 = "Persistence_Rs0_noise_additing_term/"
path = path1 * path2

include(path * "Compute_fixed_points.jl")
include(path * "Verify_no_extinction_by_noise.jl")
include(path * "Verify_R0s_additin_noise_condition.jl")
include(path * "Compute_stochastic_R0.jl")
include(path * "Verify_rs_persistence_condition.jl")
include(path * "Compute_auxiliar_constants_ci.jl")
include(path * "Compute_auxiliar_constants_rho_i.jl")
include(path * "Verify_persistence_assumption2.jl")

function Sampler_persistence_parameters(N_p)
    # Input:
    #   N_p: plant size
    # Output:
    #   parameters: vector with parameters that satisfies persistence conditions

    Test1 = false
    while Test1 == false
        beta_p = 0.1#rand(Uniform(0,1))
        r_1 = 0.01#rand(Uniform(0,1))
        b = 0.075#rand(Uniform(0,1))
        r_2 = rand(Uniform((1/2)*(-1+sqrt(2))*b,1))
        beta_v = 0.01#rand(Uniform(0,1))
        theta = 0.4#rand(Uniform(0,1))
        mu = 1.0#rand(Uniform(0,1))
        gamma = 0.06#rand(Uniform(0,1))
        gamma_f = 0.01#rand(Uniform(0,1))
        sigma_v = rand(Uniform(0,1))
        epsilon = rand(Uniform(0,1))
        sigma_L = rand(Uniform(0,1))
        sigma_I = rand(Uniform(0,1))#epsilon*sigma_L
        N_p = N_p
        N_v = mu / (gamma + gamma_f)

        par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
         beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
         gamma_f = gamma_f, sigma_L = sigma_L, sigma_I = sigma_I,
         sigma_v = sigma_v, N_v = N_v, N_p = N_p, epsilon = epsilon)

        endemic_fixed_point = Compute_fixed_points(par)
        cond0 = Verify_no_extinction_by_noise(par)
        cond1 = Verify_R0s_additin_noise_condition(par)
        cond2 = Verify_rs_persistence_condition(par)
        auxiliar_constants_c_i =
            Compute_auxiliar_constants_ci(par,endemic_fixed_point)
        auxiliar_constants_rho_i = Compute_auxiliar_constants_rho_i(par)
        cond3, auxiliar_constants_a_i = Verify_persistence_assumption2(
                auxiliar_constants_c_i,auxiliar_constants_rho_i,par,
                endemic_fixed_point
            )
        Test1 = cond0 && cond1 && cond2 && cond3
        if Test1 == true
            return par, auxiliar_constants_rho_i, auxiliar_constants_a_i,
            auxiliar_constants_c_i, endemic_fixed_point
        end
    end
end

par, auxiliar_constants_rho_i, auxiliar_constants_a_i,
    auxiliar_constants_c_i, endemic_fixed_point =
        Sampler_persistence_parameters(1)

CSV.write(path * "Parameter_Persistence_2.csv", par)
CSV.write(path * "Constant_rho_i_2.csv", auxiliar_constants_rho_i)
CSV.write(path * "Constants_a_i_2.csv", auxiliar_constants_a_i)
CSV.write(path * "Constants_c_i_2.csv",auxiliar_constants_c_i)
CSV.write(path * "endemic_fixed_point_2.csv",endemic_fixed_point)
