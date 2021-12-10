require(dplyr)
require(stringr)
require(ggplot2)
require(tidyr)
require(nmfspalette)
require(cmdstanr)
require(bridgesampling)
require(posterior)
remotes::install_github("WGGRAFY/sarla")
require(sarla)
options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)
load("./data/WareHouse_2019.RData")
#load the temperature data
ex <- new.env()
load("./data/AR1_temperature_regions_2020.RData", env=ex)
ls(ex)
temperature_by_region <- ex$tempest
names(temperature_by_region) <- ex$regions
str(WareHouse.All.Ages.Env)

#Peek at which spp has the most data
WareHouse.All.Ages.Env %>%
  group_by(common_name) %>%
  count() %>% filter(n>15000)

source("./code/sarla_model/functions/preprocess_cal_curr.R")

spp <- read.csv("code/sarla_model/process_config.csv")
model_data <- vector("list")

for(i in 1:nrow(spp)){
  processed_data <- preprocess_cal_curr(data__ = WareHouse.All.Ages.Env,
                                        common_ = spp[i,"spp"],
                                        sex_ = "F",
                                        survey_string = spp[i,"surv_string"],
                                        years_ = spp[i,"years"])


   processed_data <- process_length_data(spp_data=processed_data, common_ = spp[i,"spp"],
                      minimum_n = 10)

  #widen data so it is by row = age and columns = year
  model_data[[i]] <- processed_data %>%
    select(age_years, year, standardl) %>%
    unique() %>%
    arrange(age_years) %>%
    pivot_wider(names_from=age_years, values_from=standardl) %>%
    arrange(year) %>% mutate_all(~replace(.,is.na(.),0))

}





petrale_data <- t(model_data[[7]][,-c(1,2,14:22)])
lingcod_data <- t(model_data[[6]][,-c(1)])

#Put lingcod data into stan format
realdat <- vector("list")
realdat$xaa_observed <- realdat$laa_observed <- lingcod_data
realdat$Nages <- nrow(realdat$xaa_observed)
realdat$Nyears <- ncol(realdat$xaa_observed)
realdat$Ncohorts <- realdat$Nages + realdat$Nyears - 1
stan_dat <- plot_and_fill_data(realdat)
sarla::fit_sarla(data = stan_dat)

fit_init <- fit_stan_data(lingcod_data, init_effects = 1, year_effects = 0)
fit_annual <- fit_stan_data()
summ <- fit$summary()
unique(summ$variable)
save(fit,file="output/LingcodInitial.Rds")
readin <- load("output/LingcodAnnual.Rds")

#create a stanfit object - not working
init_stanfit <- rstan::read_stan_csv(fit$output_files())
bridge_sampler(samples=init_stanfit)
