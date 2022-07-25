library(deSolve)
library(dplyr)
library(rstan)
library(gridExtra)
library(outbreaks)
library(bayesplot)
library(data.table)
library(knitr)
library(kableExtra)
library(tidyverse)
library(plotly)
library(ggplot2)

source("load_data.R")
source("init_tomato_sir.R")
source("mcmc_stain_summary.R")
source("mcmc_post_analysis.R")
source("divergence_plots.R")
source("model_fit.R")
#
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#
data_tomato <- load_data()
onset <- data_tomato[[1]]
cum_cases <- round(data_tomato[[2]], digits = 0)
cum_cases <- unlist(cum_cases, use.names = FALSE)

####################

data_df <- read_csv("/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/resistence_data.csv")
onset_real <- round(data_df$Time_F1_Tyking)
cum_cases_real <- round(1000 * data_df$F1_Tyking)



######################

N <- nrow(onset)
pop <- 1000 # Plant population
sample_time <- 1:N

#### SIR rstan model ####
mdl <- stan_model("sirTomato3.stan")
#
#### mcmc parameters ####
n_chains <- 5
n_warmups <- 500
n_iter <- 2000
n_thin <- 50
n_total_iter = n_warmups + n_iter * n_thin
set.seed(972198)
tomato_data <- list(n_obs = N,
                     n_theta = 9,
                     n_difeq = 7,
                     n_pop = pop,
                     y = cum_cases,
                     t0 = 0,
                     ts = sample_time)
parameters <- c("y_hat", "y_init", "theta", "R_0")
time.start_nuts <- Sys.time()
nuts_fit <-
  sampling(mdl,
           data = tomato_data,
           pars = parameters,
           init = init_sir,
           chains = n_chains,
           warmup = n_warmups,
           iter = n_iter,
           thin = n_thin#,
           #control = list(adapt_delta = 0.99)
           )

time.end_nuts <- Sys.time()
duration_nuts <- time.end_nuts - time.start_nuts
parameters = c("theta", "R_0", "y_init")
nuts_fit_summary <- summary(nuts_fit, pars = parameters)$summary
print(nuts_fit_summary,
      scientific = FALSE,
      digits = 4)

sub_path_1 <- "/home/gabrielsalcedo"
sub_path_2 <- "Documentos/TYLCVD-fit"
sub_path_3 <- "Model5/runs"

prefix_time <- Sys.time()
prefix_time <- paste(date(prefix_time),
                     hour(prefix_time),
                     minute(prefix_time), sep="_")
file_name <- paste("run-", prefix_time,".RData", sep="")
runs_path <- 
  paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/") 
save.image(file=runs_path)
#
#### Post analysis ####
#
mcmcm_post_analysis(nuts_sample = nuts_fit)
#
### Model Fit ####
# Model fitted values across the observed time period
#
model_fit()
#### Divergence analysis ####
divergence_analysis()




chain1_3 <- traceplot(nuts_fit, pars = c("theta[1]", "theta[2]", "theta[3]"))
chain4_6 <- traceplot(nuts_fit, pars = c("theta[4]", "theta[5]", "theta[6]"))
chain7_9 <- traceplot(nuts_fit, pars = c("theta[7]", "theta[8]", "theta[9]"))
chain_lp <- traceplot(nuts_fit, pars = c("theta[10]", "lp__"))

file_name1_3 <- "chain1_3.pdf"
file_name4_6 <- "chain4_6.pdf"
file_name7_9 <- "chain7_9.pdf"
file_name_lp <- "chain_lp.pdf"
plot_path1_3 <- 
  paste(sub_path_1, 
        sub_path_2, 
        "Model5/plots", 
        file_name1_3, sep = "/")
ggsave(plot_path1_3)

plot_path4_6 <- 
  paste(sub_path_1, 
        sub_path_2, 
        "Model5/plots", 
        file_name4_6, sep = "/")
ggsave(plot_path4_6)

plot_path7_9 <- 
  paste(sub_path_1, 
        sub_path_2, 
        "Model5/plots", 
        file_name7_9, sep = "/")
ggsave(plot_path7_9)

plot_path_lp <- 
  paste(sub_path_1, 
        sub_path_2, 
        "Model5/plots", 
        file_name_lp, sep = "/")
ggsave(plot_path_lp)




