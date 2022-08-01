make_plots <- function(dir, species, fit_obj, pars_main){
  
  f <- paste0(dir, species, "\\")
  dir.create(f, showWarnings = FALSE)
  bayesplot::bayesplot_theme_set(theme_light())
  
  bayesplot::mcmc_areas_ridges(fit_obj$draws(pars_main))
  ggsave(paste0(f, "ridges.pdf"), width = 4, height = 5)
  
  bayesplot::mcmc_trace(fit_obj$draws(pars_main))
  
  g <- bayesplot::mcmc_pairs(fit_obj$draws(pars_main), off_diag_fun = "hex")
  g
  ggsave(paste0(f, "pairs.pdf"), plot = g, width = 8, height = 8)
  
  post_xaa <- tidybayes::gather_draws(fit_obj, xaa[i, y]) %>%
    rename(age = i, year = y)
  
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
    post_delta_c <- tidybayes::gather_draws(fit_obj, delta_c[y])
    delta_hat <- post_delta_c %>%
      group_by(y) %>%
      summarise(med = median(.value),
                lwr = quantile(.value, probs = 0.1),
                upr = quantile(.value, probs = 0.9)
      )
    g1 <- ggplot(delta_hat, aes(y, med, ymin = lwr, ymax = upr)) +
      geom_pointrange() + ggtitle("Cohort effects")
    g1
    ggsave(paste0(f, "cohort-effects.pdf"), width = 5, height = 3)
  }
  
  if (stan_dat$est_init_effects) {
    post_eta_c <- tidybayes::gather_draws(fit_obj, eta_c[y])
    eta_hat <- post_eta_c %>%
      group_by(y) %>%
      summarise(med = median(.value),
                lwr = quantile(.value, probs = 0.1),
                upr = quantile(.value, probs = 0.9)
      )
    g2 <- ggplot(eta_hat, aes(y, med, ymin = lwr, ymax = upr)) +
      geom_pointrange() + ggtitle("Init effects")
    g2
    ggsave(paste0(f, "init-effects.pdf"), width = 5, height = 3)
  }
  
  if (stan_dat$est_year_effects) {
    post_gamma_y <- tidybayes::gather_draws(fit_obj, gamma_y[y])
    gamma_hat <- post_gamma_y %>%
      group_by(y) %>%
      summarise(
        med = median(.value),
        lwr = quantile(.value, probs = 0.1),
        upr = quantile(.value, probs = 0.9)
      )
    g3 <- ggplot(gamma_hat, aes(y, med, ymin = lwr, ymax = upr)) +
      geom_pointrange() + ggtitle("Year effects")
    g3
    ggsave(paste0(f, "year-effects.pdf"), width = 5, height = 3)
  }
  
  cowplot::plot_grid(g1, g2, g3, ncol = 1L)
  ggsave(paste0(f, "all-effects.pdf"), width = 5, height = 7)
  
  
}