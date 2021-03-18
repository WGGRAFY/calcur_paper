library(TMB)
dll_name = "ss_lvb_temp"
source("code/ss_lvb_temp/helper_MV.r")
dyn.load(dynlib(paste0("code/ss_lvb_temp/",dll_name)))
comname = "lingcod"
sp = "lingcod"
dat = readRDS("data/all_length_data.RDS")
modyears = min(dat$year):max(dat$year) #1977-2018, but there are gaps in data until 2003.
nyears = length(modyears)
dat = subset(dat, common_name == comname)# & project == "Groundfish Slope and Shelf Combination Survey")
dat = subset(dat, yc %in% modyears) #only model obs where we can follow cohort from age 0
dat = subset(dat, age_years >0) #only allow estimation of of temperature effects as early as age 0! Not before!

#get tempdat, regions, tdat.obs, tdat.sd
load(file = "data/temp_objects_for_ss_lvb_temp.RData")

sp.tdat = subset(tempdat, species == comname & age < 6)
region.at.age = sapply(sort(unique(sp.tdat$age)), function(x) sort(unique(sp.tdat$region[sp.tdat$age == x])))[1] #just looking at the region at earliest age now
sp.regions = sort(unique(region.at.age))

tdat.obs = cbind(tdat.obs[,colnames(tdat.obs) %in% sp.regions])
tdat.sd = cbind(tdat.sd[,colnames(tdat.sd) %in% sp.regions])
tdat.obs = t(t(tdat.obs) - apply(tdat.obs,2,mean, na.rm = TRUE))
#for the first 5 ages, the only unique time series are for 1-3, 4-5
input = make_input(dat, edat=list(tdat.obs, tdat.sd))


#no temperature effects, constant a, b, k, Linf
temp = input
m0 = MakeADFun(temp$dat,temp$par,DLL=dll_name, map = temp$map, random = "Ecov_re")
m0 = fit.tmb.fn(m0, 3)
saveRDS(m0, file = paste0('results/ss_lvb_temp/', sp, "_m0.RDS"))


#include AR(1) for log VB k
temp = input
temp$par = m0$parList
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
m1 = MakeADFun(temp$dat,temp$par,random=c("Ecov_re", "k_re"),DLL=dll_name, map = temp$map)
m1 = fit.tmb.fn(m1, 3)
saveRDS(m1, file = paste0('results/ss_lvb_temp/', sp, "_m1.RDS"))

#temperature effect at ages 0 and AR(1) annual residuals on VB k
temp = input
x = c(NROW(temp$par$beta_Ecov_k), NCOL(temp$par$beta_Ecov_k))
temp$map$beta_Ecov_k = matrix(NA, x[1],x[2])
temp$map$beta_Ecov_k[1,1] = 1
temp$map$beta_Ecov_k = factor(temp$map$beta_Ecov_k)
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
temp$par = m1$parList
m2 = MakeADFun(temp$dat,temp$par,random=c("Ecov_re", "k_re"),DLL=dll_name, map = temp$map)
m2 = fit.tmb.fn(m2, 3)
saveRDS(m2, file = paste0('results/ss_lvb_temp/', sp, "_m2.RDS"))

#temperature effect at ages 0 and AR(1) annual residuals on VB k, and different Linf for cohorts before 2000
temp = input
temp$par = m2$parList
x = c(NROW(temp$par$beta_Ecov_k), NCOL(temp$par$beta_Ecov_k))
temp$map$beta_Ecov_k = matrix(NA, x[1],x[2])
temp$map$beta_Ecov_k[1,1] = 1
temp$map$beta_Ecov_k = factor(temp$map$beta_Ecov_k)
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
temp$dat$X_Linf = cbind(as.integer(dat$yc<2000), as.integer(dat$yc>1999))
temp$par$beta_Linf = rep(temp$par$beta_Linf,2)
m3 = MakeADFun(temp$dat,temp$par,random=c("Ecov_re", "k_re"),DLL=dll_name, map = temp$map)
m3 = fit.tmb.fn(m3, 3)
saveRDS(m3, file = paste0('results/ss_lvb_temp/', sp, "_m3.RDS"))


aic = cbind(sapply(0:3, function(x) 2*(get(paste0("m",x))$opt$obj + length(get(paste0("m",x))$opt$par))))
rownames(aic) = paste0("m",0:3)
write.csv(aic, file = paste0("results/ss_lvb_temp/", sp, "_aic.csv"))

best = rownames(aic)[which(aic == min(aic))]

cairo_pdf(paste0("results/ss_lvb_temp/", sp, '_', best, '_ecov_res.pdf'), family = "Times", height = 10, width = 10)
plot.Ecov.res.fn(get(best), years = modyears, ylab2 = colnames(tdat.obs))
dev.off()

write.csv(summary(get(best)$sdrep, "fixed"), file = paste0("results/ss_lvb_temp/", sp, "_", best, "_sdrep.csv"))

plt.yrs = min(dat$yc[dat$age<2]):max(dat$year)
cairo_pdf(paste0("results/ss_lvb_temp/", sp, '_', best, '_k_estimates.pdf'), family = "Times", height = 7, width = 10)
par(mfrow = c(1,2), mar = c(1,1,0,1), oma = c(4,4,3,0))
for(i in 1:2) 
{
  plot.kaa.fn(i, get(best), ymax = round(get.ymax.kaa.fn(1:2, get(best), years = plt.yrs),2), do.xlab = TRUE, do.ylab = i == 1, do.nobs = i == 1, years = plt.yrs)
  mtext(side = 3, outer = FALSE, line = 1, cex = 1.5, ifelse(i==1, "Age 0", "Age > 0"))
}
mtext(side = 2, outer = TRUE, line = 2, "LVB k", cex = 1.5)
mtext(side = 1, outer = TRUE, line = 2, "Year", cex = 1.5)
dev.off()

plt.yrs = min(dat$yc):max(dat$year)
cairo_pdf(paste0("results/ss_lvb_temp/",sp,'_', best, '_laa_estimates.pdf'), family = "Times", height = 10, width = 10)
par(mfrow = c(3,3), mar = c(1,1,1,1), oma = c(4,4,0,0))
for(i in 1:9) plot.laa.fn(i, ymax = round(get.ymax.laa.fn(1:9,get(best), years = plt.yrs)), mod = get(best), do.xlabs = i %in% 7:9, do.ylabs = i %in% c(1,4,7), years = plt.yrs)
mtext(side = 2, outer = TRUE, line = 2, "Length (cm)", cex = 1.5)
mtext(side = 1, outer = TRUE, line = 2, "Year", cex = 1.5)
dev.off()
