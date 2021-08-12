using CSV
using IterableTables, DataFrames, DataTables
using Distributions

path = "D:/Julia_code_for_tomate_SDE_paper-main/Extinction_Noise/"

include(path * "Verify_extinction_by_noise.jl")
include(path * "Compute_deterministic_R0.jl")
include(path * "Verify_R0d_greatherthan_one_condition.jl")


function Sampler_extinction_by_noise_parameters(N_p)
    # Input:
    #   N_p: plant size
    # Output:
    #   parameters: vector with parameters that satisfies persistence conditions

    Test1 = false
    while Test1 == false
        beta_p = rand(Uniform(0,1))
        r_1 = rand(Uniform(0,1))
        b = rand(Uniform(0,1))
        r_2 = rand(Uniform(0,1))
        beta_v = rand(Uniform(0,1))
        theta = rand(Uniform(0,1))
        mu = rand(Uniform(0,1))
        gamma = rand(Uniform(0,1))
        gamma_f = rand(Uniform(0,1))
        sigma_v = rand(Uniform(0,1))
        epsilon = rand(Uniform(0,1))
        sigma_L = rand(Uniform(0,1))
        sigma_I = rand(Uniform(0,1))
        N_p = N_p
        N_v = mu / (gamma + gamma_f)

        par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
         beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
         gamma_f = gamma_f, sigma_L = sigma_L, sigma_I = sigma_I,
         sigma_v = sigma_v, N_v = N_v, N_p = N_p, epsilon = epsilon)

        cond0 = Verify_R0d_greatherthan_one_condition(par)
        cond1 = Verify_extinction_by_noise(par)

        Test1 = cond0 && cond1
        if Test1 == true
            return par
        end
    end
end

par = Sampler_extinction_by_noise_parameters(1)

CSV.write(path * "Parameter_extinction_noise.csv", par)
