using DifferentialEquations
using Plots, PlotlyBase; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using IterableTables, DataFrames, DataTables
using StochasticDiffEq
using Distributions

path = "D:/Julia_code_for_tomate_SDE_paper-main/Persistence/"
include(path * "Dynamics.jl")

#Parameters = CSV.read(path * "Parameter_Persistence_1.csv", DataFrame)
#=
beta_p = Parameters.beta_p[1]
r_1 =  Parameters.r_1[1]
r_2 =  Parameters.r_2[1]
b =  Parameters.b[1]
beta_v =  Parameters.beta_v[1]
theta =  Parameters.theta[1]
mu =  Parameters.mu[1]
gamma =  Parameters.gamma[1]
gamma_f =  Parameters.gamma_f[1]
N_v =  Parameters.N_v[1]
sigma_L =  Parameters.sigma_L[1]
sigma_I =  Parameters.sigma_I[1]
sigma_v = Parameters.sigma_v[1]
=#

beta_p = 0.1
r_1 =  0.01
r_2 =  0.01
b =  0.075
beta_v = 0.01
theta = 0.4
mu = 1.0
gamma = 0.06
gamma_f = 0.01
N_v =  mu / (gamma + gamma_f)
sigma_L = 0.118417
sigma_I = 0.990762
sigma_v = 0.652238


u_0 = [ 1.0, 0.0, 0.0, 0.0 , 0.0]
T = 150.0
time = (0.0, T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt = 0.01
t_s = range( 0.0, T, step = 1.0)

R0 = (beta_p * beta_v * b) / ((gamma + gamma_f) * (b + r_1) * r_2)
Rs0 = R0 -
    (1 / 2) * (
        (sigma_L + sigma_I) ^ 2 -
            sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f))
        )

################################################################################
######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Det,u_0,time)
det_sol = solve(prob_det, Tsit5(), dt = dt)
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,
noise_rate_prototype=zeros(5,2))
stc_sol = solve(prob_sde_tomato_sys,SROCKC2(), dt = dt)
################################################################################
############################ PLot variables ####################################
################################################################################

title = plot(title = "R_s =$Rs0", grid = false, showaxis = false,
    bottom_margin = -50Plots.px)
p1=plot(det_sol,vars=(1),color="blue")
p1=plot!(stc_sol,vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(det_sol,vars=(2),color="blue")
p2=plot!(stc_sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(det_sol,vars=(3),color="blue")
p3=plot!(stc_sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(det_sol,vars=(4),color="blue")
p4=plot!(stc_sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(det_sol,vars=(5),color="blue")
p5=plot!(stc_sol,vars=(5),color="red",title ="Infec. v.")

plot(p1,p2,p3,p4,p5,title,layout = @layout([[A B C]; [D E F]]),label = "")

################################################################################
######################    data  Det Solution    ################################
################################################################################
det_Time = det_sol.t
det_xu = det_sol.u
det_xu_glued = hcat(det_xu...)
Xu1 = det_xu_glued[1:5:end]
Xu2 = det_xu_glued[2:5:end]
Xu3 = det_xu_glued[3:5:end]
Xu4 = det_xu_glued[4:5:end]
Xu5 = det_xu_glued[5:5:end]


det_DF1 = DataFrame(t = det_Time,S_p = Xu1, L_p =Xu2, I_p = Xu3, S_v = Xu4,
    I_v = Xu5)
det_DF1_red = det_DF1#[1:10:end,1:end]
CSV.write(path * "Det_solution_Persistence_1.csv",det_DF1_red)

################################################################################
######################    data  Sto Solution    ################################
################################################################################
stc_Time = stc_sol.t
det_yu = stc_sol.u
det_yu_glued = hcat(det_yu...)
Yu1 = det_yu_glued[1:5:end]
Yu2 = det_yu_glued[2:5:end]
Yu3 = det_yu_glued[3:5:end]
Yu4 = det_yu_glued[4:5:end]
Yu5 = det_yu_glued[5:5:end]


stc_DF1 = DataFrame(t = stc_Time, S_p = Yu1, L_p = Yu2, I_p = Yu3, S_v = Yu4,
    I_v = Yu5)
stc_DF1_red = stc_DF1[1:20:end,1:end]
CSV.write(path * "Stc_Solution_Persistence_1.csv",stc_DF1_red)

################################################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################
Datos=DataFrame()
j = 0
trajectories = 1
while j <= 10000
    monte_prob = EnsembleProblem(prob_sde_tomato_sys)
    sim = solve(monte_prob, SROCKC2(),dt = dt,EnsembleThreads(),
        trajectories = trajectories)
    component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
    component = transpose(component) #transpose to obtain any*5 data matrix
    component = vcat(component...) #to obtain shape for dataframe
    component = vcat(component...) # again do a reshape
    variables = DataFrame(component, :auto) # define first data frame
    Datos_aux = DataFrame(t = t_s, S_p = variables[:,1], I_p = variables[:,3], I_v = variables[:,5]) #only some variables
    Datos = append!(Datos, Datos_aux) #append the data in the loop
    j+=1
    println("acepted =",j)
end
CSV.write(path * "Data_persistence_additing.csv",Datos)
