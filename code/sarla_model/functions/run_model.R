run_model <- function(i, cohort_effects, init_effects, year_effects){
  # Put lingcod data into stan format
  realdat <- vector("list")
  
  SPECIES <- spp$spp[i]
  
  realdat$xaa_observed <- realdat$laa_observed <- model_data[[i]][-1,]
  realdat$Nages <- nrow(realdat$xaa_observed)
  realdat$Nyears <- ncol(realdat$xaa_observed)
  realdat$Ncohorts <- realdat$Nages + realdat$Nyears - 1
  stan_dat <- plot_and_fill_data(realdat, init_effects = init_effects,
                                 year_effects = year_effects, cohort_effects = cohort_effects,
                                 plot=T)
  names(stan_dat)[1] <- "laa_obs"
  
  stan_dat$sigma_o_prior <- c(log(0.5), 0.5) # arbitrary!?
  stan_dat$sigma_p_prior <- c(log(0.2), 0.05) # arbitrary!?
  
  # quick testing, adjust path as needed:
  #library(cmdstanr)
  # file.remove("../sarla/inst/stan/sarla")
  #mod <- cmdstan_model("../sarla/inst/stan/sarla.stan")
  
  
  # cohort effects:
  stan_dat$est_cohort_effects <- cohort_effects
  stan_dat$est_init_effects <- init_effects
  stan_dat$est_year_effects <- year_effects
  stan_dat$N_delta_c <- stan_dat$Ncohorts*cohort_effects
  stan_dat$N_gamma_y <- stan_dat$Ncohorts*year_effects
  stan_dat$N_eta_c <- stan_dat$Ncohorts*init_effects
  
  # fit <- mod$sample(
  #   data = stan_dat,
  #   chains = 4,
  #   parallel_chains = parallel::detectCores(),
  #   iter_warmup = 1000,
  #   iter_sampling = 1000,
  #   adapt_delta = 0.95,
  #   max_treedepth = 10
  # )
  
  
  
  fit <- sarla::fit_sarla(stan_dat,
                          chains = 4,
                          parallel_chains = parallel::detectCores(),
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          adapt_delta = 0.95,
                          max_treedepth = 10)
  
  post <- posterior::as_draws_df(fit$draws())
  pars <- names(post)
  pars <- pars[!grepl("raw", pars)]
  pars_main <- pars[unique(c(
    grep("_sd", pars),
    grep("sigma_", pars), grep("beta", pars), grep("sigma_", pars)
  ))]
  
  fit$summary(variables = c("sigma_p", "sigma_o", "beta", "xaa[1,10]", "xaa[2,10]", "gamma_y[1]", "gamma_y[2]", "gamma_y_sd"))
  fit$cmdstan_diagnose()
  
  tosave<- list(stan_dat, post, fit)
  save(tosave, file = paste0("SPECIES","model.RData"))
  
  
  make_plots(dir = ".\\code\\sarla_model\\plots\\", fit_obj = fit, species = SPECIES, pars_main = pars_main)
  
}