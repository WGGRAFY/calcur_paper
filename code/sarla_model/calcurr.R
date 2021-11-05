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

devtools::load_all(".")


#sablefish so lets look at that first
spp <- c("sablefish", "darkblotched rockfish", "shortbelly rockfish", "Pacific hake",
         "Pacific sanddab", "lingcod", "petrale sole")
years <- c(2003, 2003, rep(NULL, 5))
surv_string <- list()
surv_string[[1]] <- surv_string[[2]] <- c("Triennial","Combination")
surv_string[4:7] <- "Triennial|Combination"
surv_string[[3]] <- "Combination"
model_data <- vector("list")
for(i in 1:length(spp)){
  processed_data <- process_length_data(data__ = WareHouse.All.Ages.Env,
                                        common_ = spp[i],
                                        sex_ = "F",
                                        survey_string = surv_string[[i]],
                                        years_ = years[i],
                                        minimum_n = 10,
                                        plot_bool = F)

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



fit_init <- fit_stan_data(lingcod_data, init_effects = 1, year_effects = 0)
fit_annual <- fit_stan_data()
summ <- fit$summary()
unique(summ$variable)
save(fit,file="output/LingcodInitial.Rds")
readin <- load("output/LingcodAnnual.Rds")

#create a stanfit object - not working
init_stanfit <- rstan::read_stan_csv(fit$output_files())
bridge_sampler(samples=init_stanfit)
