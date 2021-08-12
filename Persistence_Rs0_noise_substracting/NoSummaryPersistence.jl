using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using DataFrames
using StochasticDiffEq


path =
    "/home/gabrielsalcedo/Github/Julia_code_for_tomate_SDE_paper/Persistence/"
include(path * "Dynamics.jl")

Parameters = CSV.read(path * "Parameter_Persistence_2.csv", DataFrame)

beta_p = Parameters.beta_p[1]
r_1 =  Parameters.r_1[1]
r_2 =  Parameters.r_2[1]
b =  Parameters.b[1]
beta_v =  Parameters.beta_v[1]
theta =  Parameters.theta[1]
mu =  Parameters.mu[1]
gamma =  Parameters.gamma[1]
N_v =  Parameters.N_v[1]
sigma_L =  Parameters.sigma_L[1]
sigma_I =  Parameters.sigma_I[1]
sigma_v = Parameters.sigma_v[1]



u_0 = [ 1.0, 0.0, 0.0, 0.03 , 0.04]
T = 100.0
time = (0.0, T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt = 0.01
t_s = range( 0.0, T, step = 1.0)


R0 = sqrt((beta_p*beta_v*b)/(gamma*(b+r_1)*r_2))
Rs0 = (R0^2) -(1/2)*((sigma_L+sigma_I)^2-sigma_v^2/(2*(beta_v+sigma_v^2+theta*gamma)))

print("Rs0=",Rs0)
#=
function F_Det(du,u,p,t)
 @inbounds begin
        du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
        du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
        du[3] = b*u[2]-r_2*u[3]
        du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
        du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
    end
    nothing
end

function F_Drift(du,u,p,t)
     @inbounds begin
        du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
        du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
        du[3] = b*u[2]-r_2*u[3]
        du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
        du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
    end
    nothing

end

function G_Diffusion(du,u,p,t)
    @inbounds begin
        du[1,1] = sigma_L*u[2]*u[1]/N_p+sigma_I*u[1]*u[3]/N_p
        du[1,2] = 0
        du[2,1] = -sigma_L*u[1]*u[2]/N_p
        du[2,2] = 0
        du[3,1] = -sigma_I*u[1]*u[3]/N_p
        du[3,2] = 0
        du[4,1] = 0
        du[4,2] = -sigma_v*u[4]
        du[5,1] = 0
        du[5,2] = -sigma_v*u[5]
    end
    nothing
end
=#
################################################################################
######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Det,u_0,time)
det_sol = solve(prob_det,Tsit5())
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,noise_rate_prototype=zeros(5,2))
stc_sol = solve(prob_sde_tomato_sys,SROCK2(),dt=dt)
################################################################################
############################ PLot variables ####################################
################################################################################

title = plot(title = "Rs0 = $Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
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


det_DF1 = DataFrame(t = det_Time,S_p = Xu1, L_p =Xu2, I_p = Xu3, S_v = Xu4, I_v = Xu5)
det_DF1_red = det_DF1#[1:10:end,1:end]
CSV.write(path * "Det_solution_Persistence.csv",det_DF1_red)

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


stc_DF1 = DataFrame(t = stc_Time, S_p = Yu1, L_p = Yu2, I_p = Yu3, S_v = Yu4, I_v = Yu5)
stc_DF1_red = stc_DF1[1:20:end,1:end]
CSV.write(path * "Stc_Solution_Persistence.csv",stc_DF1_red)
