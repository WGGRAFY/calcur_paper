make_plots <- function(dir, species, fit_obj){
  
  stan_dat <- fit_obj$stan_dat

  f <- paste0(dir, species, "\\")
  dir.create(f, showWarnings = FALSE)
  bayesplot::bayesplot_theme_set(theme_light())
  browser()
  bayesplot::mcmc_areas_ridges(fit_obj$fit$draws(fit_obj$pars_main))
  ggsave(paste0(f, "ridges.pdf"), width = 4, height = 5)
  
  bayesplot::mcmc_trace(fit_obj$fit$draws(fit_obj$pars_main))
  
  g <- bayesplot::mcmc_pairs(fit_obj$fit$draws(fit_obj$pars_main), off_diag_fun = "hex")
  g
  ggsave(paste0(f, "pairs.pdf"), plot = g, width = 8, height = 8)

  post_xaa <- tidybayes::gather_draws(fit_obj$fit, xaa[i, y]) %>%
    rename(age = i, year = y)
  
  #add actual age and birth year to the tibble
  post_xaa$ageyear <- as.numeric(dimnames(model_data[[i]])[[1]][-1])[post_xaa$age] 
  post_xaa$birthyear <- post_xaa$age+post_xaa$year-stan_dat$Nages + max(model_data[[i]][1,]) - stan_dat$Ncohorts
  
  quantile_summary <- post_xaa %>%
    # filter(y >= stan_dat$Nages) %>%
    # mutate(y = y - stan_dat$Nages + 1) %>%
    group_by(year, age) %>%
    summarise(
      lwr = quantile(.value, 0.1),
      upr = quantile(.value, 0.9),
      med = median(.value), .groups = "drop"
    )
  
  # legend_def <- c("median" = "black", "95 quantile" = "gray")
  quantile_summary %>%
    filter(!(lwr == 0 & upr == 0)) %>%
    ggplot(aes(age, med, ymin = lwr, ymax = upr)) +
    facet_wrap(~year) +
    geom_line(alpha = 1) +
    geom_ribbon(alpha = 0.5) +
    theme_minimal() +
    # scale_colour_manual(values = legend_def) +
    theme(
      legend.position = c(.90, .05),
      legend.key.size = unit(0.1, "cm"),
      legend.title = element_text(size = "6"),
      legend.text = element_text(size = "6")
    )
  ggsave(paste0(f, "dev-by-age.pdf"), width = 10, height = 10)
  yearlabels <- as.character(seq(1980,2020,10))
  
  quantile_summary %>%
    filter(!(lwr == 0 & upr == 0)) %>%
    ggplot(aes(year, med, ymin = lwr, ymax = upr)) +
    facet_wrap(~age) +
    geom_line(alpha = 1) +
    geom_ribbon(alpha = 0.5) +
    theme_minimal() +
    # scale_colour_manual(values = legend_def) +
    theme(
      legend.position = c(.90, .05),
      legend.key.size = unit(0.1, "cm"),
      legend.title = element_text(size = "6"),
      legend.text = element_text(size = "6")
    )
  ggsave(paste0(f, "dev-by-year.pdf"), width = 10, height = 10)
  
  if (stan_dat$est_cohort_effects) {
    post_delta_c <- tidybayes::gather_draws(fit_obj$fit, delta_c[y])
    delta_hat <- post_delta_c %>%
      group_by(y) %>%
      summarise(med = median(.value),
                lwr = quantile(.value, probs = 0.1),
                upr = quantile(.value, probs = 0.9)
      )
    delta_hat$y <- 1940:2018
    g1 <- ggplot(delta_hat, aes(y, med, ymin = lwr, ymax = upr)) +
      geom_pointrange() + ggtitle("Cohort effects") 
    g1
    ggsave(paste0(f, "cohort-effects.pdf"), width = 5, height = 3)
  }
  
  if (stan_dat$est_init_effects) {
    post_eta_c <- tidybayes::gather_draws(fit_obj$fit, eta_c[y])
    eta_hat <- post_eta_c %>%
      group_by(y) %>%
      summarise(med = median(.value),
                lwr = quantile(.value, probs = 0.1),
                upr = quantile(.value, probs = 0.9)
      )
    
    eta_hat$y <- 1940:2018
    g2 <- ggplot(eta_hat, aes(y, med, ymin = lwr, ymax = upr)) +
      geom_pointrange() + ggtitle("Init effects")
    g2
    ggsave(paste0(f, "init-effects.pdf"), width = 5, height = 3)
  }
  
  if (stan_dat$est_year_effects) {
    post_gamma_y <- tidybayes::gather_draws(fit_obj$fit, gamma_y[y])
    #post_gamma_y$birthyear <- post_gamma_y$age+post_gamma_y$year-stan_dat$Nages + max(model_data[[i]][1,]) - stan_dat$Ncohorts
    gamma_hat <- post_gamma_y %>%
      group_by(y) %>%
      summarise(
        med = median(.value),
        lwr = quantile(.value, probs = 0.1),
        upr = quantile(.value, probs = 0.9)
      )
    gamma_hat$y <- seq(1940,2018)
    g3 <- ggplot(gamma_hat, aes(y, med, ymin = lwr, ymax = upr)) +
      geom_pointrange() + ggtitle("Year effects")  + xlim(c(1977, 2020))
    g3
    ggsave(paste0(f, "year-effects.pdf"), width = 5, height = 3)
  }
  if(stan_dat$est_year_effects&stan_dat$est_cohort_effects&stan_dat$est_init_effects){
  cowplot::plot_grid(g1, g2, g3, ncol = 1L)
  ggsave(paste0(f, "all-effects.pdf"), width = 5, height = 7)
  }
  postpred <- tidybayes::gather_draws(fit_obj$fit,laa_postpred[i,y]) %>%
  #      ungroup() %>%
  #     pivot_wider(names_from = c(.iteration, .draw), values_from = .value) %>%
  #  filter(.chain==1) %>%
  #  select(-c(i,y,.chain,.variable))
  #%>%
    group_by(i,y) %>%
    summarise(
     med = median(.value),
      lwr = quantile(.value, probs = 0.1),
      upr = quantile(.value, probs = 0.9)
    )
  #set missing values to NA
  stan_dat$laa_obs[which(stan_dat$laa_obs==999)] <- 0
  df <- stan_dat$laa_obs
  df <- cbind(i = rownames(df), data.frame(df, row.names=NULL))

  laa_dat <- df %>% 
    pivot_longer(cols = starts_with("X"), names_to = "y",names_prefix="X") %>%
    mutate(i = as.numeric(i), y = as.numeric(y))

    #  ppc_stat_grouped(y = laa_dat$value,
    #                  as.matrix(postpred),
    #                 group = laa_dat$y,
    #                 stat = "median")
  
  
  g4 <- ggplot(data = postpred, aes(i, med, ymin = lwr, ymax = upr)) +
    facet_wrap(~y) + 
    geom_smooth(stat = "identity", size = 0.5) +
    geom_point(data= laa_dat, aes(x=i, y=value), inherit.aes=FALSE,
              alpha = 1/2) +
    ggtitle("Postpred") + theme_minimal() 
  g4
  ggsave(paste0(f, "postpred.pdf"), width = 20, height = 12)
  
  #plot(density(extract(fit, "laa_postpred")), xlim=c(1,2000),
  #     xlab="Loss", col=grey(0, 0.8),
  #     main="Predicitive distribution")
  
  
}