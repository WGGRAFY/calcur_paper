library(TMB)
dll_name = "ss_lvb_temp"
source("code/ss_lvb_temp/helper_MV.r")
dyn.load(dynlib(paste0("code/ss_lvb_temp/",dll_name)))
comname = "petrale sole"
sp = "petsole"
ecov_year1 = 1950
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

best = "m2_long"

#cairo_pdf(paste0("results/ss_lvb_temp/", sp, '_', best, '_ecov_res.pdf'), family = "Times", height = 10, width = 10)
plot.Ecov.res.fn(get(best), years = modyears, ylab2 = colnames(tdat.obs))
#dev.off()

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

load("data/petrale_inputs.Rds")

years = 1969:2018
x = list(n_years = length(years)) 
#x$n_ages = n_a_dat
x$n_ages = 17 #this is the max age in the age comp observations. n_a_dat=40 is for weight at age, maturity at age, etc.
x$n_fleets = n_f #4
#the third GB index only has values in the first 4 years, so not useful for the joint model.
x$n_indices = 3 #CPUE for fleets 1,3 and survey fleet 7
x$n_ages_pop = n_a_dat #40

#remove selectivity stuff for 3rd gb index
x$n_selblocks = x$n_fleets + 1 #first two CPUE will use same selectivity as fleets 1,3
x$selblock_models = rep(2,x$n_selblocks) #logistic
x$selblock_pointer_fleets = t(matrix(1:x$n_fleets, x$n_fleets, x$n_years))
x$selblock_pointer_indices = t(matrix(c(1,3,5), x$n_indices, x$n_years))

x$age_comp_model_fleets = rep(1, x$n_fleets)
x$n_age_comp_pars_fleets = rep(1, x$n_fleets)
x$age_comp_model_indices = rep(1, x$n_indices)
x$n_age_comp_pars_indices = rep(1, x$n_indices)

x$fracyr_SSB = rep(spawning_timing_yr_fraction, x$n_years)
x$mature = t(matrix(maturity_at_age[1 + 1:x$n_ages,3], x$n_ages, x$n_years)) #starts at 0
x$waa_pointer_fleets = 1:x$n_fleets
x$waa_pointer_totcatch = x$n_fleets + 2 #6
x$waa_pointer_indices = c(1,3,5)
x$waa_pointer_ssb = x$n_fleets + 1 #doesn't matter
x$waa_pointer_jan1 = x$n_fleets + 1 #doesn't matter

#use average weight at age of ages 17:40 in plus group. males don't get very big. Would imply larger M for males.
weights  = exp(-(16:39)*natM[1]) #equilibrium survive to age 17:40 (no fishing)
weights = weights/sum(weights) #equilibrium proportion alive at age, give survive to age 17
waa = array(NA, dim = c(x$n_fleets+x$n_indices,x$n_years, x$n_ages))
for(i in 1:(x$n_fleets+x$n_indices)){
	temp = as.matrix(subset(waa_matrix, Fleet == i & Sex == 1 & Yr %in% years)[,-(1:7)])
	waa[i,,1:(x$n_ages-1)] = temp[,1:(x$n_ages-1)]
	waa[i,,x$n_ages] = apply(temp[,(x$n_ages):NCOL(temp)],1,function(x) sum(x*weights))
	temp = as.matrix(subset(waa_matrix, Fleet == i & Sex == 2 & Yr %in% years)[,-(1:7)])
	waa[i,,1:(x$n_ages-1)] = (waa[i,,1:(x$n_ages-1)] + temp[,1:(x$n_ages-1)])/2
	waa[i,,x$n_ages] = (waa[i,,x$n_ages] + apply(temp[,(x$n_ages):NCOL(temp)],1,function(x) sum(x*weights)))/2
}
x$waa = waa
x$Maa = matrix(natM[1], x$n_years, x$n_ages)

x$agg_catch = cbind(subset(catch_matrix, year %in% years)[,-1])
x$agg_catch_sigma = cbind(subset(catch_se_matrix, year %in% years)[,-1])

########################################
#load("data/rows.Rds")
source("code/ssm_temp/get_aref_fn.r")
x$catch_paa = matrix(NA, x$n_years*x$n_fleets, x$n_ages)
for(i in 1:x$n_fleets) {
	ind = age_comp_catch$rowind[,2] == i
	tyrs = age_comp_catch$rowind[ind,1]
	temp = age_comp_catch$age_comp[ind,]
	temp[tyrs==2009,]
	temp = age_comp_catch$Nsamp$Yr[!is.na(age_comp_catch$Nsamp[,i+1])]
	temp = age_comp_catch$age_comp[rows_catch_age_comp[[i]],]

	x$catch_paa[(i-1)*x$n_years + 1:x$n_years,] = age_comp_catch$age_comp[rows_catch_age_comp[[i]],]
x$use_catch_paa = cbind(gb$use_catch_paa)
x$catch_Neff = cbind(gb$catch_Neff)
x$catch_aref = cbind(gb$catch_aref)
########################################

x$units_indices = rep(1,3)
x$fracyr_indices = catch_timing_as_fraction[c(1:2,5)] #I think this is right?
temp = CPUE_data
temp = temp[match(years, temp$year),]
x$agg_indices = cbind(temp[,1+c(1:2,5)])
temp = matrix(1, NROW(x$agg_indices), NCOL(x$agg_indices))
temp[is.na(x$agg_indices)] = 0
x$use_indices = temp
temp = CPUE_se_log
temp = temp[match(years, temp$year),]
x$agg_index_sigma = cbind(temp[,1+c(1:2,5)])
x$units_index_paa = rep(2, x$n_indices)

x$index_paa = matrix(NA, x$n_years*x$n_indices, x$n_ages)
x$index_paa[2*x$n_years + 1:x$n_years,] = as.matrix(age_comp_survey$age_comp[match(years,age_comp_survey$Nsamp$Yr),])
temp = matrix(0, x$n_years, x$n_indices)
temp[match(age_comp_survey$Nsamp$Yr,years),3] = 1
x$use_index_paa = temp
temp = matrix(0, x$n_years, x$n_indices)
temp[match(age_comp_survey$Nsamp$Yr,years),3] = age_comp_survey$Nsamp[,2]
x$index_Neff = temp

temp = get_aref_fn(x$index_paa[2*x$n_years + which(x$index_Neff[1:x$n_years,3]>0),])
x$index_aref = matrix(17, x$n_years, x$n_indices)
x$index_aref[which(x$index_Neff[1:x$n_years,3]>0),3] = temp
##########################################

x$q_lower = rep(0, x$n_indices)
x$q_upper = rep(1000, x$n_indices)

x$n_estimated_selpars = 2* x$n_selblocks
x$n_other_selpars = 0
x$other_selpars = c()
x$estimated_selpar_pointers = 1:x$n_estimated_selpars
x$other_selpar_pointers = c()
x$selpars_lower = rep(0, x$n_estimated_selpars)
x$selpars_upper = rep(x$n_ages, x$n_estimated_selpars)

x$n_NAA_sigma = 2
x$NAA_sigma_pointers =  c(1,rep(2,x$n_ages_pop-1))
x$recruit_model = 2 # random about mean

x$use_growth_model = 1
x$percentSPR = 40

m2_long = readRDS(file = "results/ssm_temp/m2_long.RDS")


dat = m2_long$env$dat
x$Ecov_obs = dat$Ecov_obs
x$use_Ecov_obs = dat$use_Ecov_obs
#x$Ecov_obs[which(is.na(x$Ecov_obs))] = -999
x$Ecov_obs_sigma = dat$Ecov_obs_sigma
#x$Ecov_obs_sigma[which(is.na(x$Ecov_obs_sigma))] = -999

#which of the years in the environmental time series (Which may be longer) are the years of the assessment model
x$model_years = match(years, 1950:2018)

x$age_obs_growth = dat$age_obs
x$year_obs_growth = dat$year_obs
x$weight = dat$weight
x$iswt = dat$iswt
x$len = dat$len
x$islen = dat$islen
x$Ecov_maxages_b = dat$Ecov_maxages_b
x$Ecov_maxages_a = dat$Ecov_maxages_a
x$Ecov_maxages_Linf = dat$Ecov_maxages_Linf
x$fit_k_LVB = dat$fit_k

x$X_Linf = dat$X_Linf
x$obs_for_Linf_y = dat$obs_for_Linf_y


y = list(mean_rec_pars = 10)
y$logit_q = rep(-8, x$n_indices)
y$log_F1 = rep(-2, x$n_fleets)
y$F_devs = matrix(0, x$n_years-1, x$n_fleets)
y$log_N1 = 10 #random walk in the rest of the ages
#y$log_N1 = rep(10, x$n_ages) #x$N1_model = 0
y$log_NAA_sigma = rep(0, x$n_NAA_sigma)
y$estimated_selpars = rep(0, x$n_estimated_selpars)
n_catch_acomp_pars = c(1,1,1,3,1,2)[x$age_comp_model_fleets[which(apply(x$use_catch_paa,2,sum)>0)]]
n_index_acomp_pars = c(1,1,1,3,1,2)[x$age_comp_model_indices[which(apply(x$use_index_paa,2,sum)>0)]]
y$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
y$index_paa_pars = rep(0, sum(n_index_acomp_pars))
y$log_NAA = matrix(10, x$n_years-1, x$n_ages_pop)
y$N1_re = rep(10, x$n_ages_pop-1)
y$log_sig_N1 = 0

par = m2_long$parList
y$beta_b = par$beta_b
y$beta_a = par$beta_a
y$beta_k_LVB = par$beta_k
y$beta_Linf = par$beta_Linf
y$t0 = par$t0
y$beta_Ecov_b = par$beta_Ecov_b
y$beta_Ecov_a = par$beta_Ecov_a
y$beta_Ecov_k_LVB = par$beta_Ecov_k
y$beta_Ecov_Linf = par$beta_Ecov_Linf
#y$beta_Ecov_Linf = rep(0,max(1,x$Ecov_maxages_Linf+1))
y$k_LVB_AR_pars = par$k_AR_pars
y$k_LVB_re = par$k_re
y$beta_sig_Linf = par$beta_sig_Linf
y$beta_sig_L_obs = par$beta_sig_L_obs
y$beta_sig_W_obs = par$beta_sig_W_obs

y$Ecov_mu = par$Ecov_mu
y$Ecov_AR_pars = par$Ecov_AR_pars
y$Ecov_re = par$Ecov_re
y$log_Ecov_obs_sig_scale = par$log_Ecov_obs_sig_scale

temp = input
x = c(NROW(temp$par$beta_Ecov_k), NCOL(temp$par$beta_Ecov_k))
temp$map$beta_Ecov_k = matrix(NA, x[1],x[2])
temp$map$beta_Ecov_k[1,1] = 1
temp$map$beta_Ecov_k = factor(temp$map$beta_Ecov_k)
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
temp$par = m1_long$parList

ssm_input = list(dat = x, par = y)
ssm_input$map = list()
ssm_input$map$beta_k_LVB = temp$map$beta_k
ssm_input$map$beta_Ecov_k_LVB = temp$map$beta_Ecov_k
map.pars = c("beta_Ecov_b", "beta_Ecov_a", "log_Ecov_obs_sig_scale", "beta_Ecov_Linf", "beta_sig_L_obs")
ssm_input$map[map.pars] = temp$map[map.pars]

ssm_input$random = c("log_NAA", "N1_re", "Ecov_re", "k_re")

setwd("code/ssm_temp")
compile("ssm_temp.cpp", "-O0 -g")
setwd("../..")
ssm_mod <- MakeADFun(ssm_input$dat,ssm_input$par,DLL="ssm_temp", random = c("log_NAA", "Ecov_re", "k_LVB_re", "k_re"), 
  map = temp$map)

