library(TMB)
source("code/ss_lvb_temp/helper_MV.r")
setwd("code/ss_lvb_temp")
dll_name = "ss_lvb_temp"
compile(paste0(dll_name,".cpp"))
dyn.load(dynlib(dll_name))
setwd("../..")
load("data/WareHouse.All.Ages.Env_filteredFeb2021.Rdata")
alldat = WareHouse.All.Ages.Env_filtered_Feb2021
alldat$month = as.numeric(substr(alldat[["datetime_utc_iso"]], 6,7))
modyears = min(alldat$year):max(alldat$year) #1977-2018, but there are gaps in data until 2003.
alldat$yc = alldat$year-alldat$age_years  #cohort year
nyears = length(modyears)
saveRDS(alldat, "data/all_length_data.RDS")

#just needed to fill in the input data. No growth will be estimated.
comname = "sablefish"
dat = subset(alldat, common_name == comname)
dat = subset(dat, yc %in% modyears) #only model obs where we can follow cohort from age 0

#get annual observations
#tempdat = read.csv("data/paul_temp_data_2020.csv", as.is = TRUE)
tempdat = read.csv("data/all_temp_data_2023.csv", as.is = TRUE)
tempdat$std_err[which(tempdat$temperature< -8.9)] = NA
tempdat$temperature <- tempdat$avg_temp
tempdat$temperature[which(tempdat$temperature< -8.9)] = NA
tempdat$region = paste0(tempdat$area, "_", tempdat$depth_bin)
regions = sort(unique(tempdat$region))
tdat.obs = tdat.sd = matrix(NA, length(modyears), length(regions))
for(i in modyears) for(j in regions) 
{
  print(c(i,j))
  if(sum(tempdat$region == j & tempdat$year == i))
  {
    tdat.obs[which(modyears == i),which(regions == j)] = tempdat$temperature[which(tempdat$region == j & tempdat$year == i)][1]
    tdat.sd[which(modyears == i),which(regions == j)] = tempdat$std_err[which(tempdat$region == j & tempdat$year == i)][1]
  }
}
colnames(tdat.obs) = colnames(tdat.sd) = regions
rownames(tdat.obs) = rownames(tdat.sd) = modyears
save(tempdat, regions, tdat.obs, tdat.sd, file = "data/temp_objects_for_ss_lvb_temp_2023.RData")

#tdat.obs = t(t(tdat.obs) - apply(tdat.obs,2,mean, na.rm = TRUE))
#for the first 5 ages, the only unique time series are for 1-3, 4-5
input = make_input(dat, edat=list(tdat.obs, tdat.sd))


#no temperature effects, no growth, just AR1 of the temperature time series...
temp = input
temp$dat$fit_growth = 0
temp$map$beta_sig_W_obs = factor(rep(NA,length(temp$par$beta_sig_W_obs)))
temp$map$beta_sig_Linf = factor(rep(NA,length(temp$par$beta_sig_Linf)))
temp$map$beta_b = factor(rep(NA,length(temp$par$beta_b)))
temp$map$beta_a = factor(rep(NA,length(temp$par$beta_a)))
temp$map$t0 = factor(rep(NA,length(temp$par$t0)))
temp$map$beta_k = factor(rep(NA,length(temp$par$beta_k)))
temp$map$beta_Linf = factor(rep(NA,length(temp$par$beta_Linf)))
Ecov_AR1_mod = MakeADFun(temp$dat,temp$par,DLL=dll_name, map = temp$map, random = "Ecov_re")

Ecov_AR1_mod = fit.tmb.fn(Ecov_AR1_mod, 3)
x = summary(Ecov_AR1_mod$sdrep, "report")
ind = matrix(which(rownames(x) == "Ecov_y"), NROW(tdat.obs), NCOL(tdat.obs))
tempcov = array(NA,dim = c(NCOL(tdat.obs), NROW(tdat.obs),NROW(tdat.obs)))
tempest = matrix(NA,NROW(tdat.obs), NCOL(tdat.obs))
for(p in 1:NCOL(tdat.obs)) {
  tempcov[p,,] = Ecov_AR1_mod$sdrep$cov[ind[,p],ind[,p]]
  tempest[,p] = x[ind[,p],1]
}
colnames(tempest) <- dimnames(tempcov)[[1]] <-  regions
rownames(tempest) <- dimnames(tempcov)[[2]] <- dimnames(tempcov)[[3]] <- modyears
tempcov[1,1:5,1:5]
temp_res <- cbind.data.frame(years = modyears, obs = c(tdat.obs), sd.obs = c(tdat.sd))
temp_res$est <- c(tempest)
temp_res$sd <- c(sapply(1:length(regions), function(x) sqrt(diag(tempcov[x,,]))))
temp_res$area <- c(sapply(regions, function(x) rep(strsplit(x,"_")[[1]][1], length(modyears))))
temp_res$depth <- c(sapply(regions, function(x) rep(strsplit(x,"_")[[1]][2], length(modyears))))
temp_res$interval_est <- qnorm(0.975) * temp_res$sd
temp_res$interval_obs <- qnorm(0.975) * temp_res$sd.obs
x <- ggplot(temp_res, aes(x = years, y = est)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = est - qnorm(0.975)*sd, ymax = est + qnorm(0.975)*sd), alpha = 0.2) +
  geom_point(aes(x = years, y = obs)) +
  geom_errorbar(aes(ymin = obs - qnorm(0.975)*sd.obs, ymax = obs + qnorm(0.975)*sd.obs)) +
  facet_grid(depth ~ area) + 
  theme_bw()
x
temp_res$area = factor(temp_res$area, levels = levels(factor(temp_res$area))[c(2:3,1)])
temp_res$depth = factor(temp_res$depth, levels = levels(factor(temp_res$depth))[c(1,3,2)])
library(ggplot2)
y <- ggplot(temp_res, aes(x = years, y = est)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = est - qnorm(0.975)*sd, ymax = est + qnorm(0.975)*sd), alpha = 0.2) +
  geom_point(aes(x = years, y = obs)) +
  geom_errorbar(aes(ymin = obs - qnorm(0.975)*sd.obs, ymax = obs + qnorm(0.975)*sd.obs)) +
  facet_grid(depth ~ area) + 
  theme_bw() + ylab("Bottom temperature") + xlab("Year")
cairo_pdf('results/AR1_temperature_regions_2023.pdf', family = "Times", height = 10, width = 10)
y
dev.off()


x = summary(Ecov_AR1_mod$sdrep, "fixed")
AR1pars = AR1pars_sd = matrix(NA,3,NCOL(tdat.obs))
AR1pars[1,] = x[rownames(x) == "Ecov_mu",1]
AR1pars_sd[1,] = x[rownames(x) == "Ecov_mu",2]
x = summary(Ecov_AR1_mod$sdrep, "report")
AR1pars[2,] = x[rownames(x) == "Ecov_phi",1]
AR1pars_sd[2,] = x[rownames(x) == "Ecov_phi",2]
AR1pars[3,] = x[rownames(x) == "Ecov_sig",1]
AR1pars_sd[3,] = x[rownames(x) == "Ecov_sig",2]
save(regions, tempest,tempcov, AR1pars, AR1pars_sd, file = 'results/AR1_temperature_regions_2023.RData')

