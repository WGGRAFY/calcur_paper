# setwd("../../")

require(dplyr)
require(stringr)
require(ggplot2)
theme_set(theme_light())
require(tidyr)
require(nmfspalette)
require(cmdstanr)
require(bridgesampling)
require(posterior)
remotes::install_github("WGGRAFY/sarla", force = TRUE, ref = "addcovars")
require(sarla)
require(loo)

load("./data/WareHouse_2019.RData")
# load the temperature data
ex <- new.env()
load("./data/AR1_temperature_regions_2020.RData", env = ex)
load("./data/temp_objects_for_ss_lvb_temp.RData", env = ex)
match_area <- read.csv("./results/ss_lvb_temp/first_age_area_depth.csv") %>%
  mutate(region = paste(area,depth, sep="_"))

ls(ex)
temperature_by_region <- ex$tempest
colnames(temperature_by_region) <- ex$regions
str(WareHouse.All.Ages.Env)

# Peek at which spp has the most data
WareHouse.All.Ages.Env %>%
  group_by(common_name) %>%
  count() 

#8527 total for lingcod
#Look at lingcod
WareHouse.All.Ages.Env %>%
  filter(common_name=="lingcod", latitude_dd > 40+1/6) 

source("./code/sarla_model/functions/preprocess_cal_curr.R")

source("./code/sarla_model/functions/make_plots.R")
source("./code/sarla_model/functions/run_model.R")

spp <- read.csv("code/sarla_model/process_config.csv")
spp <- left_join(x = spp, match_area, by = c("spp" = "comname"))
model_data <- vector("list")

for (i in seq_len(length(spp$spp))[-1]) {
  processed_data <- preprocess_cal_curr(
    data__ = WareHouse.All.Ages.Env,
    common_ = spp[i, "spp"],
    sex_ = "F",
    survey_string = spp[i, "surv_string"],
    years_ = spp[i, "years"]
  )


  processed_data <- process_length_data(
    spp_data = processed_data,
    common_ = spp[i, "spp"],
    minimum_n = 10
  )

  # widen data so it is by row = age and columns = year
  model_data[[i]] <- processed_data %>%
    select(age_years, year, standardl) %>%
    unique() %>%
    arrange(age_years) %>%
    pivot_wider(names_from = age_years, values_from = standardl) %>%
    arrange(year) %>%
    mutate_all(~ replace(., is.na(.), 999)) %>%
    t()
  
  cohort_temp <- temperature_by_region[,spp[i,"region"]]
  

fit_all <- run_model(i = i, 1L, 1L, 1L)
fit_year <- run_model(i = i, 0L, 0L, 1L)
fit_cohort <- run_model(i = i, 1L, 0L, 0L, cohort_cov = cohort_temp)
fit_init <- run_model(i = i, 0L, 1L, 0L)
 fit_null <- run_model(i = i, 0L, 0L, 0L)
  
 
}

strings <- paste0("./code/sarla_model/output/",rep(spp$spp, each = 5),rep(c("y1i1c1"
,"y0i1c0", "y1i0c0", "y0i0c1", "y0i0c0"),7),"model.RData")

loo_list <- vector("list")
loo_table <- matrix(nrow = 7, ncol=5)
j <- k <- 1
L <- 10
for(i in seq_len(length(strings))){
  tosave <- get(load(strings[i]))
  fit <- tosave$fit
  loo_list[[i]]<- fit$loo(variables = "log_lik")
  loo_table[j,k] <- loo_list[[i]]$looic
  if(k==5){
    j <- j+1
    k <- 1
  } else{
    k <- k+1
  }
}
save(loo_table, file="lootable.RData")
load("lootable.RData")

k_thres = 0.7
approx_elpds_1sap <- rep(NA, N)
N <- length(model_data[[i]])
stanfit <- rstan::read_stan_csv(fit$output_files())
stanfit$call <- function(x, ...){
  x$formula
}
L = 10
# initialize the process for i = L
past <- 1:L
oos <- L + 1
df <- model_data[[i]]
df_past <- df[, past, drop = FALSE]
df_oos <- df[,c(past, oos), drop = FALSE]
model_data[[i]] <- df_past
fit_past <- run_model(i = i, 0L, 0L, 0L)
###
#Stopping point 12/6 - there's no oos argument to fit$loo
# when passed to the brms::log_lik function via ...
# this is passed to brms::prepare_predictions to tell it to compute out-of-sample
#rather than in-sample predictions. Not sure how to do this when using fit$loo rather than 
#log_lik
loglikpast <- fit_past$loo(variables = "log_lik", newdata = df_oos)
, oos = oos)
approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])

# iterate over i > L
i_refit <- L
refits <- L
ks <- NULL
for (i in (L + 1):(N - 1)) {
  past <- 1:i
  oos <- i + 1
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]

  loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
  
  logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
  psis_obj <- suppressWarnings(psis(logratio))
  k <- pareto_k_values(psis_obj)
  ks <- c(ks, k)
  if (k > k_thres) {
    # refit the model based on the first i observations
    i_refit <- i
    refits <- c(refits, i)
    fit_past <- update(fit_past, newdata = df_past, recompile = FALSE)
    loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
    approx_elpds_1sap[i + 1] <- log_mean_exp(loglik[, oos])
  } else {
    lw <- weights(psis_obj, normalize = TRUE)[, 1]
    approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
  }
} 

plot_ks <- function(ks, ids, thres = 0.6) {
    dat_ks <- data.frame(ks = ks, ids = ids)
    ggplot(dat_ks, aes(x = ids, y = ks)) + 
      geom_point(aes(color = ks > thres), shape = 3, show.legend = FALSE) + 
      geom_hline(yintercept = thres, linetype = 2, color = "red2") + 
      scale_color_manual(values = c("cornflowerblue", "darkblue")) + 
      labs(x = "Data point", y = "Pareto k") + 
      ylim(-0.5, 1.5)
  }
