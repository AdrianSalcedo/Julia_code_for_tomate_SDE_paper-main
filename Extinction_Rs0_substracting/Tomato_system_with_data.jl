using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using DataFrames
using StochasticDiffEq

path = "/media/gabrielsalcedo/ADATA HV620S/"
path2 ="Julia_code_for_tomate_SDE_paper-main/Extinction_Rs0_substracting"
path = path1 * path2

include(path * "Dynamics.jl")
circles = CSV.read(path * "datos_circle.csv", DataFrame)
squarts = CSV.read(path * "datos_squares.csv", DataFrame)

beta_p = 0.1
r_1 =  0.01
r_2 =  0.01
b =  0.075
beta_v = 0.01
theta = 0.4
mu = 1.0
gamma = 0.06
gamma_f = 0.0

sigma_L = 0.621265
sigma_I = 0.990762
sigma_v = 0.5

u_0 = [1.0,0.0,0.0,0.3,0.2]
T = 70.0
time = (0.0,T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt=0.001

R0_without_fumigation =
    (beta_p * beta_v * b) / ( (gamma + gamma_f) * ( b + r_1) * r_2)

######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_Holt = ODEProblem(F_Holt,u_0,time)
det_Holt = solve(prob_Holt, Tsit5(), dt = dt)
prob_det = ODEProblem(F_Drift,u_0,time)
det_sol = solve(prob_det, Tsit5(), dt = dt)

gamma_f = 0.553183
mu = 0.0687
beta_v = 0.15444
R0_with_fumigation =
    (beta_p * beta_v * b) / ( (gamma + gamma_f) * ( b + r_1) * r_2)

R0 = (beta_p * beta_v * b) / ((gamma + gamma_f) * (b + r_1) * r_2)
Rs0 = R0 -
    (1 / 2) * (
        (sigma_L + sigma_I) ^ 2 -
            sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f))
        )

########################## Stochastis Solution #################################

prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,
noise_rate_prototype=zeros(5,2))
stc_sol = solve(prob_sde_tomato_sys,SROCKC2(),dt = dt)
################################################################################
############################ PLot variables ####################################
p1 = plot(
    det_Holt, vars=(3), color="blue", title = "Holt System", legend = false
        )
p1 = plot!(
    circles.t, circles.y, seriestype = :scatter, color="red", legend = false
        )
p1 = plot!(
    squarts.t, squarts.y, seriestype = :scatter, color="black", legend = false
        )
p2 = plot(
    det_sol, vars=(3), color="purple", title = "Adimensional System"
        )
p2 = plot!(
    circles.t, circles.y, seriestype = :scatter, color="red", legend = false
        )
p2 = plot!(
    squarts.t, squarts.y, seriestype = :scatter, color="black", legend = false
        )
p2 = plot!( stc_sol, vars=(3), color="darkred", legend = false)
ylims!(-0.01, 1.05)
plot(p1, p2, layout = @layout([A B]), label = "")

################################################################################
######################    data  Det Solution    ################################
################################################################################

det_Time = det_Holt.t
det_xu = det_Holt.u
det_xu_glued = hcat(det_xu...)
Xu1 = det_xu_glued[1:5:end]
Xu2 = det_xu_glued[2:5:end]
Xu3 = det_xu_glued[3:5:end]
Xu4 = det_xu_glued[4:5:end]
Xu5 = det_xu_glued[5:5:end]

det_DF1 = DataFrame(t = det_Time,S_p = Xu1, L_p =Xu2, I_p = Xu3, S_v = Xu4, I_v = Xu5)
det_DF1_red = det_DF1#[1:10:end,1:end]
CSV.write(path * "Det_Rs0_ext_Holt.csv",det_DF1_red)

det_Time = det_sol.t
det_xu = det_sol.u
det_xu_glued = hcat(det_xu...)
Xu1 = det_xu_glued[1:5:end]
Xu2 = det_xu_glued[2:5:end]
Xu3 = det_xu_glued[3:5:end]
Xu4 = det_xu_glued[4:5:end]
Xu5 = det_xu_glued[5:5:end]

det_DF1 = DataFrame(t = det_Time,S_p = Xu1, L_p =Xu2, I_p = Xu3, S_v = Xu4, I_v = Xu5)
det_DF1_red = det_DF1#[1:10:end,1:end]
CSV.write(path * "Det_Rs0_ext_good.csv",det_DF1_red)

################################################################################
######################    data  Stc Solution    ################################
################################################################################
stc_Time = stc_sol.t
det_yu = stc_sol.u
det_yu_glued = hcat(det_yu...)
Yu1 = det_yu_glued[1:5:end]
Yu2 = det_yu_glued[2:5:end]
Yu3 = det_yu_glued[3:5:end]
Yu4 = det_yu_glued[4:5:end]
Yu5 = det_yu_glued[5:5:end]


stc_DF1 = DataFrame(t = stc_Time, S_p = Yu1, L_p = Yu2, I_p = Yu3, S_v = Yu4, I_v = Yu5)
stc_DF1_red = stc_DF1[1:20:end,1:end]
CSV.write(path *"Stc_Rs0_ext_good.csv",stc_DF1_red)
