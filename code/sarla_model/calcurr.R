setwd("../../")

require(dplyr)
require(stringr)
require(ggplot2)
theme_set(theme_light())
require(tidyr)
require(nmfspalette)
require(cmdstanr)
require(bridgesampling)
require(posterior)
detach("package:sarla", unload = TRUE)

## This is the script version that has temperature covariates

devtools::load_all("pers-git/sarla")
remotes::install_github("WGGRAFY/sarla", force = TRUE, ref = "addcovars",
                        dependencies = FALSE)

require(sarla)
require(loo)


load("pers-git/calcur_paper/data/WareHouse_2019.RData")
# load the temperature data
ex <- new.env()
load("pers-git/calcur_paper/results/AR1_temperature_regions_2023.RData", env = ex)

load("pers-git/calcur_paper/data/temp_objects_for_ss_lvb_temp.RData", env = ex)
match_area <- read.csv("pers-git/calcur_paper/results/ss_lvb_temp/first_age_area_depth.csv") %>%
  mutate(region = paste(area,depth, sep="_"))

ls(ex)
temperature_by_region <- ex$tempest
temperature_by_region <- data.frame(temperature_by_region)

colnames(temperature_by_region) <- dimnames(ex$tempest)[[2]]
temperature_by_region$year <- 1977:2018

#source helper functions that are in the functions/ folder
source("pers-git/calcur_paper/code/sarla_model/functions/preprocess_cal_curr.R")
source("pers-git/calcur_paper/code/sarla_model/functions/make_plots.R")
source("pers-git/calcur_paper/code/sarla_model/functions/run_model.R")

# process config sets up the model inputs
spp <- read.csv("pers-git/calcur_paper/code/sarla_model/process_config.csv")
spp <- left_join(x = spp, match_area, by = c("spp" = "comname"))
model_data <- vector("list")

#Loop through the data, process it, and fit the model
for (i in 1:length(spp$spp)) {
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
  model_data[[i]] <- processed_data |>
    select(age_years, year, standardl) |>
    unique() |>
    arrange(age_years) |>
    pivot_wider(names_from = age_years, values_from = standardl) |>
    arrange(year) |>
    mutate_all(~ replace(., is.na(.), 999)) |>
    t()

  ind <- which(temperature_by_region$year %in% model_data[[i]]["year",])
  cohort_temp <- temperature_by_region[ind,spp[i,"region"]]


fit_all <- run_model(i = i, 1L, 1L, 1L, cov_effects = 0, cohort_cov = cohort_temp)
make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_all,
          species = paste0(spp$spp[i], "y", 1L,"t",0, "i", 1L,
                           "c", 1L, "a", 1L))

fit_year <- run_model(i = i, 0L, 0L, year_effects = 1L,
                      cov_effects = 0, cohort_cov = rep(0,length(cohort_temp)))


make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_year,
           species = paste0(spp$spp[i], "y", 1L, "t", 0, "i", 0L,
                            "c", 0L, "a", 1L))

fit_year_temp <- run_model(i = i, 0L, 0L, year_effects = 1L,
                    cov_effects = 1, cohort_cov = cohort_temp)


make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_year_temp,
           species = paste0(spp$spp[i], "y", 1L,"t", 1L, "i", 0L,
                            "c", 0L, "a", "1L", "dec"))

fit_cohort <- run_model(i = i, cohort_effects = 1L,
                        year_effects = 0L, init_effects = 0L,
                        cov_effects = 0L, cohort_cov = cohort_temp)
make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_cohort,
           species = paste0(spp$spp[i], "y", 0L,0, "i", 0L,
                            "c", 1L, "a", 1L))


fit_init <- run_model(i = i, 0L, 1L, 0L, cov_effects = 0, cohort_cov = rep(0,length(cohort_temp)))

make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_init,
          species = paste0(spp$spp[i], "y", 0L, "i", 1L,"t",0,
                           "c", 0L, "a", 1))

fit_init_temp <- run_model(i = i, 0L, 1L, 0L, cov_effects = 1, cohort_cov = cohort_temp)
make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_init_temp,
           species = paste0(spp$spp[i], "y", 0L, "i", 1L,"t",1L,
                            "c", 0L, "a", 1))
fit_null <- run_model(i = i, 0L, 0L, 0L, cov_effects = 0, cohort_cov = cohort_temp)

make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = fit_null,
           species = paste0(spp$spp[i], "y", 0L, "i", 0L,"t",0L,
                            "c", 0L, "a", 1))

}



strings <- paste0("./pers-git/calcur_paper/code/sarla_model/output/",rep(spp$spp, each = 7),rep(c("y1i1c1t0"
,"y0i1c0t0","y0i1c0t1", "y1i0c0t0", "y1i0c0t1", "y0i0c1t0","y0i0c0t0"),7),"model.RData")

##Plot only
j<- k <- 1
for(i in 1:length(strings)){
  tosave <- get(load(strings[i]))
  make_plots(dir = "pers-git\\calcur_paper\\code\\sarla_model\\plots\\", fit_obj = tosave,
             species = paste0(spp$spp[j], "model", k))
  if(k==7){
    j <- j+1
    k <- 1
  } else{
    k <- k+1
  }
}


loo_list <- vector("list")
loo_table <- matrix(nrow = 7, ncol=7)
diags <- vector("list")
j <- k <- 1
for(i in 1:length(strings)){
  tosave <- get(load(strings[i]))
  # if('try-error' %in% class(tosave)){
  #   loo_list[[i]] <- NULL
  #   loo_table[j,k] <- NA
  #   print(i)
  # } else{
  fit <- tosave$fit
  loo_list[[i]]<- fit$loo(variables = "log_lik")
  #diags[[i]] <- fit$sampler_diagnostics()
  loo_table[j,k] <- loo_list[[i]]$looic
#  }
  if(k==7){
    j <- j+1
    k <- 1
  } else{
    k <- k+1
  }
}
save(loo_table, file="lootable.RData")
load("lootable.RData")


apply(loo_table, 1,which.min)


### Everything below this is for LFO, which I am trying to avoid having to figure out.

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
