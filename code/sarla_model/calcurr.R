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
# remotes::install_github("WGGRAFY/sarla", ref="missingdata")
require(sarla)

load("./data/WareHouse_2019.RData")
# load the temperature data
ex <- new.env()
load("./data/AR1_temperature_regions_2020.RData", env = ex)
ls(ex)
temperature_by_region <- ex$tempest
names(temperature_by_region) <- ex$regions
str(WareHouse.All.Ages.Env)

# Peek at which spp has the most data
WareHouse.All.Ages.Env %>%
  group_by(common_name) %>%
  count() %>%
  filter(n > 15000)

source("./code/sarla_model/functions/preprocess_cal_curr.R")

spp <- read.csv("code/sarla_model/process_config.csv")
model_data <- vector("list")

for (i in seq_len(length(spp$spp))) {
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
    mutate_all(~ replace(., is.na(.), 999))
}

petrale_data <- t(model_data[[7]][, -c(1, 2, 22)])
lingcod_data <- t(model_data[[6]][, -c(1,13,14)])

# Put lingcod data into stan format
realdat <- vector("list")

SPECIES <- "lingcod"
SPECIES <- "petrale"

if (SPECIES == "petrale") {
  realdat$xaa_observed <- realdat$laa_observed <- petrale_data
} else if (SPECIES == "lingcod") {
  realdat$xaa_observed <- realdat$laa_observed <- lingcod_data
} else {
  stop("Species not found")
}
realdat$Nages <- nrow(realdat$xaa_observed)
realdat$Nyears <- ncol(realdat$xaa_observed)
realdat$Ncohorts <- realdat$Nages + realdat$Nyears - 1
stan_dat <- plot_and_fill_data(realdat, init_effects = 0, plot=T)
names(stan_dat)[1] <- "laa_obs"

stan_dat$sigma_o_prior <- c(log(0.5), 0.5) # arbitrary!?
stan_dat$sigma_p_prior <- c(log(0.2), 0.2) # arbitrary!?

# quick testing, adjust path as needed:
library(cmdstanr)
# file.remove("../sarla/inst/stan/sarla")
mod <- cmdstan_model("../sarla/inst/stan/sarla.stan")

# year effects:
stan_dat$N_eta_c <- 0L

# cohort effects:
stan_dat$est_cohort_effects <- 1L
stan_dat$est_init_effects <- 0L
stan_dat$est_year_effects <- 0L
stan_dat$N_delta_c <- stan_dat$Ncohorts
stan_dat$N_gamma_y <- 0L
stan_dat$N_eta_c <- 0L

# init effects:
stan_dat$est_init_effects <- 1L
stan_dat$est_year_effects <- 0L
stan_dat$est_cohort_effects <- 0L
stan_dat$N_delta_c <- 0L
stan_dat$N_gamma_y <- 0L
stan_dat$N_eta_c <- stan_dat$Ncohorts

# year and cohort effects??
stan_dat$est_cohort_effects <- 1L
stan_dat$est_init_effects <- 0L
stan_dat$est_year_effects <- 1L
stan_dat$N_delta_c <- stan_dat$Ncohorts
stan_dat$N_gamma_y <- stan_dat$Ncohorts
stan_dat$N_eta_c <- 0L

# all 3!?
stan_dat$est_cohort_effects <- 1L
stan_dat$est_init_effects <- 1L
stan_dat$est_year_effects <- 1L
stan_dat$N_delta_c <- stan_dat$Ncohorts
stan_dat$N_gamma_y <- stan_dat$Ncohorts
stan_dat$N_eta_c <- stan_dat$Ncohorts

fit <- mod$sample(
  data = stan_dat,
  chains = 4,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 10
)
fit

post <- posterior::as_draws_df(fit$draws())
pars <- names(post)
pars <- pars[!grepl("raw", pars)]
pars_main <- pars[unique(c(
  grep("_sd", pars),
  grep("sigma_", pars), grep("beta", pars), grep("sigma_", pars)
))]

fit$summary(variables = c("sigma_p", "sigma_o", "beta", "xaa[1,10]", "xaa[2,10]", "gamma_y[1]", "gamma_y[2]", "gamma_y_sd"))
fit$cmdstan_diagnose()

f <- paste0("code/sarla_model/plots/", SPECIES, "/")
dir.create(f, showWarnings = FALSE)
bayesplot::bayesplot_theme_set(theme_light())

bayesplot::mcmc_areas_ridges(fit$draws(pars_main))
ggsave(paste0(f, "ridges.pdf"), width = 4, height = 5)

bayesplot::mcmc_trace(fit$draws(pars_main))

g <- bayesplot::mcmc_pairs(fit$draws(pars_main), off_diag_fun = "hex")
g
ggsave(paste0(f, "pairs.pdf"), plot = g, width = 8, height = 8)

post_xaa <- tidybayes::gather_draws(fit, xaa[i, y]) %>%
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
  post_delta_c <- tidybayes::gather_draws(fit, delta_c[y])
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
  post_eta_c <- tidybayes::gather_draws(fit, eta_c[y])
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
  post_gamma_y <- tidybayes::gather_draws(fit, gamma_y[y])
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


# --------------

# fit <- sarla::fit_sarla(
#   data = stan_dat,
#   chains = 4,
#   iter = 2000,
#   parallel_chains = parallel::detectCores(),
#   adapt_delta = 0.95
# )
#
# summ <- fit$summary()
# require(shinystan)
# stanfit <- rstan::read_stan_csv(fit$output_files())
# launch_shinystan(stanfit)
#
#
# summ[1:81,]
# unique(summ$variable)
# save(fit, file = "code/sarla_model/output/LingcodAnnual.Rds")
# posterior <- fit$draws()
# gamma_draws <- subset_draws(posterior, "gamma_y")
# #Make some plots
# require(bayesplot)
#
# plot_title <- ggtitle("Posterior distributions",
#                       "with medians and 80% intervals")
# mcmc_intervals(gamma_draws,
#            # pars = c(paste("gamma_y[",1:24,"]",sep="")),
#            prob = 0.8) + plot_title
#
# readin <- load("code/sarla_model/output/LingcodAnnual.Rds")
#
# # create a stanfit object - not working
# init_stanfit <- rstan::read_stan_csv(fit$output_files())
# bridge_sampler(samples = init_stanfit)
