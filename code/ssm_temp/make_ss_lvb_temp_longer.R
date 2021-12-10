#This helps get parameters and data formated for the growth model component to hindcast LVB k and Ecov back in time far enough for the assessment model

library(TMB)
dll_name = "ss_lvb_temp"
source("code/ss_lvb_temp/helper_MV.r")
dyn.load(dynlib(paste0("code/ss_lvb_temp/",dll_name)))
comname = "petrale sole"
sp = "petsole"
#first years of assessment model: 1969
#oldest age of fish in population: 40
ecov_year1 = 1969-45 #just to make sure it goes back far enough
dat = readRDS("data/all_length_data.RDS")
#get tempdat, regions, tdat.obs, tdat.sd
load(file = "data/temp_objects_for_ss_lvb_temp.RData")

ldata_years = min(tempdat$year):max(tempdat$year) #1977-2018, but there are gaps in data until 2003.
nyears = length(ldata_years)
dat = subset(dat, common_name == comname)# & project == "Groundfish Slope and Shelf Combination Survey")
dat = subset(dat, yc %in% ldata_years) #only model obs where we can follow cohort from age 0
dat = subset(dat, age_years >0) #only allow estimation of of temperature effects as early as age 0! Not before!
modyears = ecov_year1:max(ldata_years)


sp.tdat = subset(tempdat, species == comname & age < 6)
region.at.age = sapply(sort(unique(sp.tdat$age)), function(x) sort(unique(sp.tdat$region[sp.tdat$age == x])))[1] #just looking at the region at earliest age now
sp.regions = sort(unique(region.at.age))

tdat.obs = cbind(tdat.obs[,colnames(tdat.obs) %in% sp.regions])
tdat.obs = cbind(c(rep(NA, length(ecov_year1:(min(tempdat$year)-1))), tdat.obs))
tdat.sd = cbind(tdat.sd[,colnames(tdat.sd) %in% sp.regions])
tdat.sd = cbind(c(rep(NA, length(ecov_year1:(min(tempdat$year)-1))), tdat.sd))
tdat.obs = t(t(tdat.obs) - apply(tdat.obs,2,mean, na.rm = TRUE))
#for the first 5 ages, the only unique time series are for 1-3, 4-5
input = make_input(dat, edat=list(tdat.obs, tdat.sd))
saveRDS(input, file = "code/ssm_temp/growth_input.RDS")
#no temperature effects, constant a, b, k, Linf
temp = input
m0_long = MakeADFun(temp$dat,temp$par,DLL=dll_name, map = temp$map, random = "Ecov_re")
m0_long = fit.tmb.fn(m0_long, 3)
m0 = readRDS(file = paste0('results/ss_lvb_temp/', sp, "_m0.RDS"))

#include AR(1) for log VB k
temp = input
temp$par = m0_long$parList
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
m1_long = MakeADFun(temp$dat,temp$par,random=c("Ecov_re", "k_re"),DLL=dll_name, map = temp$map)
m1_long = fit.tmb.fn(m1_long, 3)
m1 = readRDS(file = paste0('results/ss_lvb_temp/', sp, "_m1.RDS"))
saveRDS(m1_long, file = "results/ssm_temp/m1_long.RDS")

#temperature effect at ages 0 and AR(1) annual residuals on VB k
temp = input
x = c(NROW(temp$par$beta_Ecov_k), NCOL(temp$par$beta_Ecov_k))
temp$map$beta_Ecov_k = matrix(NA, x[1],x[2])
temp$map$beta_Ecov_k[1,1] = 1
temp$map$beta_Ecov_k = factor(temp$map$beta_Ecov_k)
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
temp$par = m1_long$parList
m2_long = MakeADFun(temp$dat,temp$par,random=c("Ecov_re", "k_re"),DLL=dll_name, map = temp$map)
m2_long = fit.tmb.fn(m2_long, 3)
m2 = readRDS(file = paste0('results/ss_lvb_temp/', sp, "_m2.RDS"))
saveRDS(m2_long, file = "results/ssm_temp/m2_long.RDS")
#m2_long = readRDS("results/ssm_temp/m2_long.RDS")
best = "m2_long"

#cairo_pdf(paste0("results/ssm_temp/example_ar1.pdf"), family = "Times", height = 5, width = 5)
png(filename = "results/ssm_temp/example_ar1.png", width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
plot.Ecov.res.fn(get(best), years = modyears, ylab2 = "", xrange = c(1970,2020))
dev.off()

cairo_pdf(paste0("results/ssm_temp/petrale_k_age_0.pdf"), family = "Times", height = 5, width = 5)
nyears = length(modyears)
plt.yrs = modyears
par(mfrow = c(1,1), mar = c(1,1,0,1), oma = c(4,4,3,0))
ymax = round(get.ymax.kaa.fn(1, get(best), years = plt.yrs),2)
plot.kaa.fn(1, get(best), ymax = ymax, do.xlab = TRUE, do.ylab = TRUE, do.nobs = FALSE, years = plt.yrs)
text(max(plt.yrs), ymax, adj = c(1,1), cex = 2, "Age 0")
mtext(side = 2, outer = TRUE, line = 2, "LVB k", cex = 1.5)
mtext(side = 1, outer = TRUE, line = 2, "Year", cex = 1.5)
dev.off()

plt.yrs = min(dat$yc[dat$age<2]):max(dat$year)
#cairo_pdf(paste0("results/ss_lvb_temp/", sp, '_', best, '_k_estimates.pdf'), family = "Times", height = 7, width = 10)
nyears = length(modyears)
par(mfrow = c(1,2), mar = c(1,1,0,1), oma = c(4,4,3,0))
for(i in 1:2) 
{
  plot.kaa.fn(i, get(best), ymax = round(get.ymax.kaa.fn(1:2, get(best), years = plt.yrs),2), do.xlab = TRUE, do.ylab = i == 1, do.nobs = i == 1, years = plt.yrs)
  mtext(side = 3, outer = FALSE, line = 1, cex = 1.5, ifelse(i==1, "Age 0", "Age > 0"))
}
mtext(side = 2, outer = TRUE, line = 2, "LVB k", cex = 1.5)
mtext(side = 1, outer = TRUE, line = 2, "Year", cex = 1.5)
#dev.off()

plt.yrs = min(dat$yc):max(dat$year)
#cairo_pdf(paste0("results/ss_lvb_temp/",sp,'_', best, '_laa_estimates.pdf'), family = "Times", height = 10, width = 10)
par(mfrow = c(3,3), mar = c(1,1,1,1), oma = c(4,4,0,0))
for(i in 1:9) plot.laa.fn(i, ymax = round(get.ymax.laa.fn(1:9,get(best), years = plt.yrs)), mod = get(best), do.xlabs = i %in% 7:9, do.ylabs = i %in% c(1,4,7), years = plt.yrs)
mtext(side = 2, outer = TRUE, line = 2, "Length (cm)", cex = 1.5)
mtext(side = 1, outer = TRUE, line = 2, "Year", cex = 1.5)
#dev.off()
