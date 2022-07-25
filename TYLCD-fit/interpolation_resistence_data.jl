using CSV, DataFrames
using DataInterpolations, Plots
path = "/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/"

N_p = 1000
df = CSV.read(path * "resistence_data.csv", DataFrame)
#interpolation:
t1_end = round(df.Time_LA_1582[end])
t1_init = round(df.Time_LA_1582[1])
timeLine1 = t1_init:t1_end
t2_end = round(df.Time_F1_Tyking[end-1])
t2_init = round(df.Time_F1_Tyking[1])
timeLine2 = t2_init:t2_end

scaled_i_p_La1582 = N_p * df.LA_1582
scaled_i_p_f1tyking = N_p * df.F1_Tyking

obs_t1 = df.Time_LA_1582
obs_t2 = df.Time_F1_Tyking

itp_La1582 = BSplineApprox(scaled_i_p_La1582, obs_t1,3,8,:ArcLen,:Average)
#(scaled_i_p_La1582, obs_t1)
itp_f1tyking = BSplineApprox(scaled_i_p_f1tyking, obs_t2,3,7,:ArcLen,:Average)

scatter(obs_t1, scaled_i_p_La1582)
plot!(itp_La1582, legend = false)
scatter!(obs_t2, scaled_i_p_f1tyking)
plot!(itp_f1tyking, legend = false)

interpolated_infected_La1582 = [itp_La1582(t) for t in timeLine1]
interpolated_infected_f1tyking = [itp_f1tyking(t) for t in timeLine2]

#difference_La1582 = diff(interpolated_infected_La1582)
#difference_f1tyking = diff(interpolated_infected_f1tyking)

#plot(difference_La1582)
#plot(difference_f1tyking)

Resistence_interpolate_df1 = DataFrame(
    Time_LA_1582 = timeLine1,
    interpolated_La1582 = interpolated_infected_La1582
    )
Resistence_interpolate_df2 = DataFrame(
    Time_F1_Tyking = timeLine2,
    interpolated_f1tyking = interpolated_infected_f1tyking
    )
plot(
    Resistence_interpolate_df1.Time_LA_1582,
    Resistence_interpolate_df1.interpolated_La1582
    )
scatter!(
    Resistence_interpolate_df1.Time_LA_1582,
    Resistence_interpolate_df1.interpolated_La1582
    )
plot!(
    Resistence_interpolate_df2.Time_F1_Tyking,
    Resistence_interpolate_df2.interpolated_f1tyking
    )
scatter!(
    Resistence_interpolate_df2.Time_F1_Tyking,
    Resistence_interpolate_df2.interpolated_f1tyking
    )

CSV.write(
    path * "Resistence_interpolated_data1.csv",
    Resistence_interpolate_df1
)
CSV.write(
    path * "Resistence_interpolated_data2.csv",
    Resistence_interpolate_df2
)
plot(
    Resistence_interpolate_df1.Time,
    Resistence_interpolate_df.interpolated_La1582
    )
scatter!(
    Resistence_interpolate_df.Time,
    Resistence_interpolate_df.interpolated_La1582
)
