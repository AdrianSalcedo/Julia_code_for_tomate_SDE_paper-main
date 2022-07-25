library(bayesplot)
library(ggplot2)
library(GGally)
library(dplyr)

mcmcm_post_analysis <- function(nuts_sample = nuts_fit){
    #' Create and save post analysis plots
    #' @param nuts_sample a list with the output of the sampling method
    #' @return trace, pairs, and density plots for evaluate the quality
    #' of the estimation
    #'
    posts <-  rstan::extract(nuts_fit)
    mod_diagnostics  <- rstan::get_sampler_params(nuts_fit)
    #
    # Check for divergent transitions
    #
    rstan::check_divergences(nuts_fit)
    posterior <- as.array(nuts_fit)
    color_scheme_set("viridisE")
    #
    # Markov chain traceplots
    #
    pars1 <- c("theta[1]",
            "theta[2]",
            "theta[3]",
            "theta[4]"
            )
    pars2 <- c("theta[5]",
            "theta[6]",
            "theta[7]",
            "theta[8]",
            "theta[9]")
    pars3 <-c(
            "y_init[1]",
            "y_init[2]",
            "y_init[3]",
            "y_init[4]"
    )
    pars4<-c("y_init[5]",
            "y_init[6]",
            "y_init[7]"
            )
    p1 <- mcmc_trace_highlight(posterior, pars = "lp__")
    p2 <- mcmc_trace_highlight(posterior, pars = "R_0")
    p31 <- mcmc_trace(posterior, pars = pars1)
    p32 <- mcmc_trace(posterior, pars = pars2)
    p33 <- mcmc_trace(posterior, pars = pars3)
    p34 <- mcmc_trace(posterior, pars = pars4)
    #
#
    p4 <- mcmc_pairs(posterior,
                        pars = pars1,
                        off_diag_args = list(size = 0.75),
                        np_style = pairs_style_np(
                                            div_color = "red",
                                            div_shape = 4,
                                            div_size = 1,
                                            div_alpha = 1,
                                            td_color = "yellow2",
                                            td_shape = 3,
                                            td_size = 1,
                                            td_alpha = 0.5))
    p41 <- mcmc_pairs(posterior,
                     pars = pars1,
                     off_diag_args = list(size = 0.75),
                     np_style = pairs_style_np(
                       div_color = "red",
                       div_shape = 4,
                       div_size = 1,
                       div_alpha = 1,
                       td_color = "yellow2",
                       td_shape = 3,
                       td_size = 1,
                       td_alpha = 0.5))
    p42 <- mcmc_pairs(posterior,
                     pars = pars2,
                     off_diag_args = list(size = 0.75),
                     np_style = pairs_style_np(
                       div_color = "red",
                       div_shape = 4,
                       div_size = 1,
                       div_alpha = 1,
                       td_color = "yellow2",
                       td_shape = 3,
                       td_size = 1,
                       td_alpha = 0.5))
    p43 <- mcmc_pairs(posterior,
                     pars = pars3,
                     off_diag_args = list(size = 0.75),
                     np_style = pairs_style_np(
                       div_color = "red",
                       div_shape = 4,
                       div_size = 1,
                       div_alpha = 1,
                       td_color = "yellow2",
                       td_shape = 3,
                       td_size = 1,
                       td_alpha = 0.5))
    p44 <- mcmc_pairs(posterior,
                      pars = pars4,
                      off_diag_args = list(size = 0.75),
                      np_style = pairs_style_np(
                        div_color = "red",
                        div_shape = 4,
                        div_size = 1,
                        div_alpha = 1,
                        td_color = "yellow2",
                        td_shape = 3,
                        td_size = 1,
                        td_alpha = 0.5))
#
    p51 <- stan_dens(nuts_fit, pars = pars1, separate_chains = TRUE)
    p52 <- stan_dens(nuts_fit, pars = pars2, separate_chains = TRUE)
    p53 <- stan_dens(nuts_fit, pars = pars3, separate_chains = TRUE)
    p54 <- stan_dens(nuts_fit, pars = pars4, separate_chains = TRUE)
#
    p61 <- mcmc_intervals(posterior, pars = pars1)
    p62 <- mcmc_intervals(posterior, pars = pars2)
    p63 <- mcmc_intervals(posterior, pars = pars3)
    p64 <- mcmc_intervals(posterior, pars = pars4)
#
    sub_path_1 <- "/home/gabrielsalcedo"
    sub_path_2 <- "Documentos/TYLCVD-fit"
    sub_path_3 <- "Model5/plots"
    file_name_1 <- "trace_posterior.pdf"
    file_name_2 <- "trace_r_zero.pdf"
    
    file_name_31 <- "trace_parameters1.pdf"
    file_name_32 <- "trace_parameters2.pdf"
    file_name_33 <- "trace_parameters3.pdf"
    file_name_34 <- "trace_parameters4.pdf"
    
    file_name_41 <- "pairs1.pdf"
    file_name_42 <- "pairs2.pdf"
    file_name_43 <- "pairs3.pdf"
    file_name_44 <- "pairs4.pdf"
    
    file_name_51 <- "posterior_densities1.pdf"
    file_name_52 <- "posterior_densities2.pdf"
    file_name_53 <- "posterior_densities3.pdf"
    file_name_54 <- "posterior_densities4.pdf"
    
    file_name_61 <- "intervals1.pdf"
    file_name_62 <- "intervals2.pdf"
    file_name_63 <- "intervals3.pdf"
    file_name_64 <- "intervals4.pdf"
#    
    plot_path_1 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_1, sep = "/")
    plot_path_2 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_2, sep = "/")
    
    plot_path_31 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_31, sep = "/")
    plot_path_32 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_32, sep = "/")
    plot_path_33 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_33, sep = "/")
    plot_path_34 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_34, sep = "/")
    
    plot_path_41 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_41, sep = "/")
    plot_path_42 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_42, sep = "/")
    plot_path_43 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_43, sep = "/")
    plot_path_44 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_44, sep = "/")
    
    plot_path_51 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_51, sep = "/")
    plot_path_52 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_52, sep = "/")
    plot_path_53 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_53, sep = "/")
    plot_path_54 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_54, sep = "/")
    
    plot_path_61 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_61, sep = "/")
    plot_path_62 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_62, sep = "/")
    plot_path_63 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_63, sep = "/")
    plot_path_64 <- 
      paste(sub_path_1, sub_path_2, sub_path_3, file_name_64, sep = "/")
    
    ggsave(plot_path_1, plot = p1)
    ggsave(plot_path_2, plot = p2)
    
    ggsave(plot_path_31 , plot = p31)
    ggsave(plot_path_32 , plot = p32)
    ggsave(plot_path_33 , plot = p33)
    ggsave(plot_path_34 , plot = p34)
    
    ggsave(plot_path_41,
            plot = p41,
            width = 11,
            height = 8,
            dpi = 100,
            units = "in")
    ggsave(plot_path_42,
           plot = p42,
           width = 11,
           height = 8,
           dpi = 100,
           units = "in")
    ggsave(plot_path_43,
           plot = p43,
           width = 11,
           height = 8,
           dpi = 100,
           units = "in")
    ggsave(plot_path_44,
           plot = p44,
           width = 11,
           height = 8,
           dpi = 100,
           units = "in")
#    
    ggsave(plot_path_51, plot = p51)
    ggsave(plot_path_52, plot = p52)
    ggsave(plot_path_53, plot = p53)
    ggsave(plot_path_54, plot = p54)
    
    ggsave(plot_path_61, plot = p61)
    ggsave(plot_path_62, plot = p62)
    ggsave(plot_path_63, plot = p63)
    ggsave(plot_path_64, plot = p64)
    }
