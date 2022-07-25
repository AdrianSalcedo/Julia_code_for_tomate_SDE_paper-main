library(lubridate)
mcmcm_stain_summary <- function(fit = nuts_fit){
    # Set initial values:
    parameters = c(
        "theta", "R_0", "y_init"
                   )
    nuts_fit_summary <- summary(nuts_fit, pars = parameters)$summary
    df_pars <- data.frame(nuts_fit_summary)
    names(df_pars)[4] <-"2.5%"
    names(df_pars)[5] <-"25%"
    names(df_pars)[6] <-"50%"
    names(df_pars)[7] <-"75%"
    names(df_pars)[8] <-"97.5%"
    colnames(df_pars) <- c('mean', 'se_mean','sd',
                            '2.5%', '25%','50%',
                            '75%' ,'97.5%', 'n_eff', 'Rhat')
    row.names(df_pars) <- c("beta_p", "r_1", "r_2", "b", "beta_v", 
                            "gamma", "gamma_f", "theta_1", "mu",
                            "R_zero",
                            "Sp0", "Lp0", "Ip0", 
                            "Sv0", "Iv0"
                            )
    print(df_pars)
    prefix_time <- Sys.time()
    prefix_time <- paste(as.Date(prefix_time),
                    hour(prefix_time),
                    minute(prefix_time), sep="_")
    sub_path_1 <- "/home/gabrielsalcedo"
    sub_path_2 <- "Documentos/TYLCVD-fit"
    sub_path_3 <- "Model5/estimation"
    file_name <- paste("parameters", prefix_time,".csv", sep="")
    path <- paste(sub_path_1,
                    sub_path_2,
                    sub_path_3,
                    file_name,
                    sep = "/")
    path_last <- paste(sub_path_1,
                        sub_path_2,
                        sub_path_3, 
                        "last_est.csv",
                            sep = "/")
    write.csv(df_pars, path)
    write.csv(df_pars, path_last)
}
