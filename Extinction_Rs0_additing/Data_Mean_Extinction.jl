using Plots; plotly()
using CSV
using IterableTables, DataFrames, DataTables
using Statistics

path1 = "/home/gabrielsalcedo/Github/"
path2 = "Julia_code_for_tomate_SDE_paper/Extinction_Rs0/"
path = path1 * path2
Data = CSV.read(path * "Data_Rs0_3.csv", DataFrame)


Data_Mean = DataFrame()

for i in 0:500
    Sub = Data[Data.t .== i, :]
    Mean_t = mean(Sub.t)
    Mean_S_p = mean(Sub.S_p)
    Mean_I_p = mean(Sub.I_p)
    Mean_I_v = mean(Sub.I_v)
    Quartile_S_p = quantile!(Sub.S_p, [0.05, 0.25, 0.50, 0.75, 0.95])
    Quartile_I_p = quantile!(Sub.I_p, [0.05, 0.25, 0.50, 0.75, 0.95])
    Quartile_I_v = quantile!(Sub.I_v, [0.05, 0.25, 0.50, 0.75, 0.95])
    Data_aux = DataFrame(
    t = Mean_t,
    Q05_S_p = Quartile_S_p[1],
    Q25_S_p = Quartile_S_p[2],
    Q50_S_p = Quartile_S_p[3],
    Mean_S_p = Mean_S_p,
    Q75_S_p = Quartile_S_p[4],
    Q95_S_p = Quartile_S_p[5],
    Q05_I_p = Quartile_I_p[1],
    Q25_I_p = Quartile_I_p[2],
    Q50_I_p = Quartile_I_p[3],
    Mean_I_p = Mean_I_p,
    Q75_I_p = Quartile_I_p[4],
    Q95_I_p = Quartile_I_p[5],
    Q05_I_v = Quartile_I_v[1],
    Q25_I_v = Quartile_I_v[2],
    Q50_I_v = Quartile_I_v[3],
    Mean_I_v = Mean_I_v,
    Q75_I_v = Quartile_I_v[4],
    Q95_I_v = Quartile_I_v[5],
      )
    Data_Mean = append!(Data_Mean,Data_aux)
end

CSV.write(path * "Data_Mean_Extinction_rs_3.csv",Data_Mean)
