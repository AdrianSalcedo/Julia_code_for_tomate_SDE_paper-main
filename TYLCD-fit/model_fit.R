library(lubridate)
model_fit <- function(fit = nuts_fit, dates = onset, dates_real = onset_real){
    posts <-  rstan::extract(nuts_fit)
    mod_diagnostics  <- rstan::get_sampler_params(nuts_fit)
    color_scheme_set("viridisE")
    fit_CIS <- posts$y_hat[, ,7]
    fit_SIR <- fit_CIS
    median_I = apply(fit_SIR, 2, median)
    low_I = apply(fit_SIR, 2, quantile, probs = c(0.025))
    high_I = apply(fit_SIR, 2, quantile, probs = c(0.975))
    df_sample_N = data.frame(cum_cases, dates$Time_F1_Tyking)
     df_sample_N_real = data.frame(cum_cases_real, dates_real)
    df_fit_CIS = data.frame(median_I, low_I, high_I, dates)
    #
    sub_path_1 <- "/home/gabrielsalcedo"
    sub_path_2 <- "Documentos/TYLCVD-fit"
    sub_path_3 <- "Model5/data"
    file_name <- "data.Rda"
    file_name_real <- "data_real.Rda"
    data_path <- 
        paste(sub_path_1, 
              sub_path_2, 
              sub_path_3, 
              file_name, 
              sep = "/")
    data_path_real <- 
      paste(sub_path_1, 
            sub_path_2, 
            sub_path_3, 
            file_name_real, 
            sep = "/")
    save(df_sample_N, file = data_path)
    save(df_sample_N_real, file = data_path_real)
        file_name <- "df_I_det_Poiss.Rda"
    data_path <- 
        paste(sub_path_1, 
              sub_path_2, 
              sub_path_3, 
              file_name, 
              sep = "/")
    save(df_fit_CIS, file = data_path)
    write.csv(df_sample_N, "/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/Data_LA_1582.csv")
    write.csv(df_sample_N_real, "/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/Data_real_LA_1582.csv")
    write.csv(df_fit_CIS, "/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/BinNeg_LA_1582.csv")
    #
    plt <- ggplot(df_sample_N,
                    aes(x = dates.Time_F1_Tyking, y = cum_cases)) +
            geom_ribbon(aes(x = dates.Time_F1_Tyking,
                            ymin = low_I,
                            ymax = high_I),
                        fill = "orange",
                        alpha = 0.6)  +
            geom_line(data = df_fit_CIS,
                        aes(x = Time_F1_Tyking,
                            y = median_I,
                            color = "Median"),
                            size = 1.3) +
      geom_point(data = df_sample_N_real, 
                 aes(x = dates_real,y = cum_cases_real), 
                 color = "red") +
            geom_point(shape = 1,
                        size = 2,
                        (aes(color = "Data"))) +
            scale_colour_manual(name = '',
                                values = c('Data' = 'black',
                                           'Median' = 'darkorange3')) +
            guides(colour = 
                            guide_legend(override.aes = 
                                            list(shape = c(16, NA),
                                                    linetype = c(0, 1)))) +
            labs(x = "Time (days)",
                    y = "Cumulative Infected Cases")# + 
      #xlim(0, 50) #+ ylim(0,30)
            #+
            #  scale_x_continuous(limits = c(0, 32)) +
            # scale_y_continuous(limits = c(0, 1),
            #          breaks = seq(from = 0, to = 32, by = 1)) #+
            #theme_bw() + theme(text = element_text(size = 20))
    #
    show(plt)
    sub_path_1 <- "/home/gabrielsalcedo"
    sub_path_2 <- "Documentos/TYLCVD-fit"
    sub_path_3 <- "Model5/plots"
    file_name <- "Tomato_data_begining_fit.pdf"
    plot_path <- 
        paste(sub_path_1, 
              sub_path_2, 
              sub_path_3, 
              file_name, sep = "/")
    ggsave(plot_path)
}
