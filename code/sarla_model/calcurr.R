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
remotes::install_github("WGGRAFY/sarla", force = TRUE)
require(sarla)

load("./data/WareHouse_2019.RData")
# load the temperature data
ex <- new.env()
load("./data/AR1_temperature_regions_2020.RData", env = ex)
ls(ex)
temperature_by_region <- ex$tempest
colnames(temperature_by_region) <- ex$regions
str(WareHouse.All.Ages.Env)

# Peek at which spp has the most data
WareHouse.All.Ages.Env %>%
  group_by(common_name) %>%
  count() %>%
  filter(n > 15000)

source("./code/sarla_model/functions/preprocess_cal_curr.R")

source("./code/sarla_model/functions/make_plots.R")
source("./code/sarla_model/functions/run_model.R")

spp <- read.csv("code/sarla_model/process_config.csv")
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

 # fit_all <- run_model(i = i, 1L, 1L, 1L)
#  fit_year <- run_model(i = i, 0L, 0L, 1L)
#  fit_cohort <- run_model(i = i, 1L, 0L, 0L)
#  fit_init <- run_model(i = i, 0L, 1L, 0L)
 fit_null <- run_model(i = i, 0L, 0L, 0L)
  
 
}

strings <- paste0("./code/sarla_model/output/",rep(spp$spp, each = 5),rep(c("y1i1c1"
,"y0i1c0", "y1i0c0", "y0i0c1", "y0i0c0"),7),"model.RData")

loo_list <- vector("list")
loo_table <- matrix(nrow = 7, ncol=5)
j <- k <- 1
for(i in seq_len(length(strings))){
  tosave <- get(load(strings[i]))
  fit <- tosave$fit
  loo_list[[i]]<- fit$loo(cores = 2)
  loo_table[j,k] <- loo_list[[i]]$looic
  if(k==5){
    j <- j+1
    k <- 1
  } else{
    k <- k+1
  }
}
save(loo_table, file="lootable.RData")
