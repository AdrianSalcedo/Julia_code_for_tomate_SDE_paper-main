path1 = "/home/gabrielsalcedo/Github/Julia_code_for_tomate_SDE_paper/"
path2 = "Persistence_Rs0_noise_additing_term/"
path = path1 * path2
include(path * "Compute_fixed_points.jl")
include(path * "Verify_no_extinction_by_noise.jl")
include(path * "Compute_deterministic_R0.jl")
include(path * "Compute_stochastic_R0.jl")
include(path * "Verify_rs_persistence_condition.jl")
include(path * "Compute_auxiliar_constants_ci.jl")
include(path * "Compute_auxiliar_constants_rho_i.jl")
include(path * "Verify_persistence_assumption2.jl")

Persistence_parameters =
    CSV.read(path * "Parameter_Persistence.csv", DataFrame)
Constanst_rho_i = CSV.read(path * "Constant_rho_i.csv", DataFrame)
Constanst_a_i = CSV.read(path * "Constants_a_i.csv", DataFrame)
Constanst_c_i = CSV.read(path * "Constants_c_i.csv", DataFrame)
Endemic_fixed_point = CSV.read(path * "endemic_fixed_point.csv", DataFrame)

R0 = Compute_deterministic_R0(Persistence_parameters)
Rs0 = Compute_stochastic_R0(Persistence_parameters)
Condition_no_extinction_By_noise =
    Verify_no_extinction_by_noise(Persistence_parameters)
    Verify_R0s_additin_noise_condition(Persistence_parameters)
    Verify_rs_persistence_condition(Persistence_parameters)
    Verify_persistence_assumption2(
            Constanst_c_i,Constanst_rho_i,Persistence_parameters,
            Endemic_fixed_point)
