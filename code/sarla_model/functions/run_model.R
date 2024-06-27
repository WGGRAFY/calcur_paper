run_model <- function(i, cohort_effects, init_effects, year_effects, 
                      cohort_cov = NULL, cov_effects = NULL){
  # Put lingcod data into stan format
  SPECIES <- spp$spp[i]
  filename <- paste0(".\\pers-git\\calcur_paper\\code\\sarla_model\\output\\",SPECIES,
                     "y", year_effects,
                     "i", init_effects,
                     "c", cohort_effects, "t", cov_effects,
                     "model.RData")

  # if(file.exists(filename)){
  #   load(filename)
  #   stan_dat <- tosave$stan_dat
  #   fit <- tosave$fit
  #   post <- tosave$post
  # } else{


  realdat <- vector("list")

  realdat$xaa_observed <- realdat$laa_observed <- model_data[[i]][-1,]
  realdat$Nages <- nrow(realdat$xaa_observed)
  realdat$Nyears <- ncol(realdat$xaa_observed)
  realdat$Ncohorts <- realdat$Nages + realdat$Nyears - 1
  #We need a temp for every cohort which means going back longer in time
  #now calculating inside stan
  #extra_temps <- rep(mean(cohort_temp), realdat$Ncohorts - length(cohort_temp))
  realdat$cov_effect <- cohort_cov
  stan_dat <- plot_and_fill_data(realdat, init_effects = init_effects,
                                 year_effects = year_effects, 
                                 cohort_effects = cohort_effects,
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
  stan_dat$est_cov_effects <- cov_effects
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
  post_mat <- posterior::as_draws_matrix(fit$draws())


  pars <- names(post)
  pars <- pars[!grepl("raw", pars)]
  pars_main <- pars[unique(c(
    grep("_sd", pars),
    grep("sigma_", pars), grep("beta", pars), grep("sigma_", pars)
  ))]
  
  tosave<- list("stan_dat" = stan_dat, "post" = post, "fit" = fit, "pars_main" = pars_main)
  
  save(tosave, file = filename)
  
  summary_vars <- c("sigma_p", "sigma_o", "beta", "xaa[1,10]", "xaa[2,10]")
  if(year_effects==1L){
    summary_vars <- c(summary_vars, "gamma_y[1]", "gamma_y[2]", "gamma_y_sd")
  }
  if(cohort_effects==1L){
    summary_vars <- c(summary_vars, "delta_c[1]", "delta_c[2]", "delta_c_sd")
  }
  if(init_effects==1L){
    summary_vars <- c(summary_vars, "eta_c[1]", "eta_c[2]", "eta_c_sd")
  }

  fit$cmdstan_diagnose()
  
  return(tosave)
  
}
