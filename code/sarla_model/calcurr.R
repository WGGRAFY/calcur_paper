# setwd("../../")

require(dplyr)
require(stringr)
require(ggplot2)
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
realdat$xaa_observed <- realdat$laa_observed <- lingcod_data
realdat$Nages <- nrow(realdat$xaa_observed)
realdat$Nyears <- ncol(realdat$xaa_observed)
realdat$Ncohorts <- realdat$Nages + realdat$Nyears - 1
stan_dat <- plot_and_fill_data(realdat, init_effects = 0, plot=T)
names(stan_dat)[1] <- "laa_obs"

# quick testing, adjust path as needed:
library(cmdstanr)
file.remove("../sarla/inst/stan/sarla")
mod <- cmdstan_model("../sarla/inst/stan/sarla.stan")

stan_dat$N_eta_c <- 0L
stan_dat$sigma_o_prior <- c(log(0.5), 0.5)
stan_dat$sigma_p_prior <- c(log(0.2), 0.2)

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

# Look at fitted model
fit$summary(variables = c("sigma_p", "sigma_o", "beta", "xaa[1,10]", "xaa[2,10]", "gamma_y[1]", "gamma_y[2]", "X0[1]"))
fit$cmdstan_diagnose()

post <- posterior::as_draws_df(fit$draws())
pars <- names(post)
pars <- pars[!grepl("raw", pars)]
pars_main <- pars[unique(c(
  grep("_sd", pars),
  grep("sigma_", pars), grep("beta", pars), grep("sigma_", pars)
))]

bayesplot::mcmc_areas_ridges(fit$draws(pars_main))
bayesplot::mcmc_trace(fit$draws(pars_main))
bayesplot::mcmc_pairs(fit$draws(pars_main), off_diag_fun = "hex")

post_xaa <- tidybayes::gather_draws(fit, xaa[i, y]) %>%
  rename(age = i, year = y)

quantile_summary <- post_xaa %>%
  # filter(y >= stan_dat$Nages) %>%
  # mutate(y = y - stan_dat$Nages + 1) %>%
  group_by(year, age) %>%
  summarise(
    lwr = quantile(.value, 0.1),
    upr = quantile(.value, 0.9),
    med = median(.value)
  )

legend_def <- c("median" = "black", "95 quantile" = "gray")
quantile_summary %>%
  filter(!(lwr == 0 & upr == 0)) %>%
  ggplot(aes(age, med, ymin = lwr, ymax = upr)) +
  facet_wrap(~year) +
  geom_line(alpha = 1, aes(colour = "red")) +
  geom_ribbon(alpha = 0.5) +
  theme_minimal() +
  scale_colour_manual(values = legend_def) +
  theme(
    legend.position = c(.90, .05),
    legend.key.size = unit(0.1, "cm"),
    legend.title = element_text(size = "6"),
    legend.text = element_text(size = "6")
  )

# if (stan_dat$est_cohort_effects) {
#   post_delta_c <- tidybayes::gather_draws(fit, delta_c[y])
#   delta_hat <- post_delta_c %>%
#     group_by(y) %>%
#     summarise(med = median(.value))
# }

# if (stan_dat$est_init_effects) {
#   post_eta_c <- tidybayes::gather_draws(fit, eta_c[y])
#   eta_hat <- post_eta_c %>%
#     group_by(y) %>%
#     summarise(med = median(.value))
#   eta_lwr <- post_eta_c %>%
#     group_by(y) %>%
#     summarise(lwr = quantile(.value, 0.025))
#   eta_upr <- post_eta_c %>%
#     group_by(y) %>%
#     summarise(upr = quantile(.value, 0.975))
# }

if (stan_dat$est_year_effects) {
  post_gamma_y <- tidybayes::gather_draws(fit, gamma_y[y])
  gamma_hat <- post_gamma_y %>%
    group_by(y) %>%
    summarise(
      med = median(.value),
      lwr = quantile(.value, probs = 0.1),
      upr = quantile(.value, probs = 0.9)
    )
  ggplot(gamma_hat, aes(y, med, ymin = lwr, ymax = upr)) +
    geom_pointrange()
}


# --------------

# fit <- sarla::fit_sarla(
#   data = stan_dat,
#   chains = 4,
#   iter = 2000,
#   parallel_chains = parallel::detectCores(),
#   adapt_delta = 0.95
# )

summ <- fit$summary()
require(shinystan)
stanfit <- rstan::read_stan_csv(fit$output_files())
launch_shinystan(stanfit)


summ[1:81,]
unique(summ$variable)
save(fit, file = "code/sarla_model/output/LingcodAnnual.Rds")
posterior <- fit$draws()
gamma_draws <- subset_draws(posterior, "gamma_y")
#Make some plots
require(bayesplot)

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_intervals(gamma_draws,
           # pars = c(paste("gamma_y[",1:24,"]",sep="")),
           prob = 0.8) + plot_title

readin <- load("code/sarla_model/output/LingcodAnnual.Rds")

# create a stanfit object - not working
init_stanfit <- rstan::read_stan_csv(fit$output_files())
bridge_sampler(samples = init_stanfit)
