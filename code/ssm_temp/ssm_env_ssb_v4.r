mat.dat <- read.csv('../../cod_gb_19702014s_0yr_mvavg.csv', header = TRUE)
mat.dat = mat.dat[which(mat.dat$yc %in% 1963:2014),]
mat.dat = mat.dat[which(mat.dat$AGE > 0),]
fall.temp.dat = read.csv("../../1963-2014-fall-temp-for-state-space.csv", header = TRUE)

source("~/work/retro/R/read.asap3.dat.fn.r")
source("~/work/retro/R/asap3_2_ssm_10_rev.dat.fn.r")
source("~/work/FATE/SSM/R/write.ssm_10_rev.dat.fn.r")
source("~/work/FATE/SSM/R/read.ssm_10_rev.dat.fn.r")
source("~/work/FATE/SSM/tmb/get_aref_fn.r")
source("~/work/SSM_M/R/ssm_10_rev_data_2_ssm_m_data_fn.r")
gbcod.asap3 = read.asap3.dat.fn("~/work/FATE/SSMSP/R/2015_COD_GB_MOD_asap_RUN6_DIFF_ZERO_2012.dat")
gbcod.ssm10 = asap3_2_ssm_10_rev.dat.fn(gbcod.asap3)
write.ssm_10_rev.dat.fn(gbcod.ssm10, "gbcod_ssm_10_rev.dat")
gb = ssm_10_rev_data_2_ssm_m_data_fn("gbcod_ssm_10_rev.dat")

#gomcod.asap3 = read.asap3.dat.fn("~/work/FATE/SSMSP/R/2015_COD_GM_MOD_ASAP_M02.dat")
#gomcod.ssm10 = asap3_2_ssm_10_rev.dat.fn(gomcod.asap3)
#write.ssm_10_rev.dat.fn(gomcod.ssm10, "gomcod_ssm_10_rev.dat")
#gom = ssm_10_rev_data_2_ssm_env_data_fn("gomcod_ssm_10_rev.dat")

x = list(n_years = gb$n_years) 
x$n_ages = gb$n_ages 
x$n_fleets = gb$n_fleets
#the third GB index only has values in the first 4 years, so not useful for the joint model.
x$n_indices = gb$n_indices-1
#dfo, autumn, (spr41_w), spr36_w

#remove selectivity stuff for 3rd gb index
x$n_selblocks = gb$n_selblocks - 1
x$selblock_models = c(gb$selblock_models[-5])
x$selblock_pointer_fleets = cbind(gb$selblock_pointer_fleets)
x$selblock_pointer_indices = cbind(gb$selblock_pointer_indices[,1:3])

x$age_comp_model_fleets = c(gb$age_comp_model_fleets)
x$n_age_comp_pars_fleets = c(gb$n_age_comp_pars_fleets)
x$age_comp_model_indices = c(gb$age_comp_model_indices[-3])
x$n_age_comp_pars_indices = c(gb$n_age_comp_pars_indices[-3])

x$fracseason_SSB = x$fracyr_SSB = cbind(gb$fracyr_SSB) #n_years x n_stocks
x$mature = gb$mature
x$waa_pointer_fleets = c(gb$waa_pointer_fleets)
x$waa_pointer_totcatch = c(gb$waa_pointer_totcatch) #n_regions
#the third index only has values in the first 4 years, so not very useful.
x$waa_pointer_indices = c(gb$waa_pointer_indices[-3])
x$waa_pointer_ssb = c(gb$waa_pointer_ssb) #n_stocks
x$waa_pointer_jan1 = c(gb$waa_pointer_jan1) #n_stocks
x$waa = gb$waa
x$Maa = cbind(gbcod.asap3$dat$M)

x$agg_catch = cbind(gb$agg_catch)
x$agg_catch_sigma = cbind(gb$agg_catch_sigma)
x$catch_paa = gb$catch_paa
x$use_catch_paa = cbind(gb$use_catch_paa)
x$catch_Neff = cbind(gb$catch_Neff)
x$catch_aref = cbind(gb$catch_aref)

x$units_indices = c(gb$units_indices[-3])
x$fracyr_indices = cbind(gb$fracyr_indices[,-3])
x$agg_indices = cbind(gb$agg_indices[,-3])
x$use_indices = cbind(gb$use_indices[,-3])
x$agg_index_sigma = cbind(gb$agg_index_sigma[,-3])
temp = as.matrix(read.csv("~/work/FATE/SSMSP/R/gbcod_survey_cvs.csv", header = TRUE)[,-1])/100
x$agg_index_sigma[,1] = sqrt(log(temp[,"DFO"]^2+1)) 
x$agg_index_sigma[,2] = sqrt(log(temp[,"Fall"]^2+1)) 
x$agg_index_sigma[,3] = sqrt(log(temp[,"Spring"]^2+1)) 
x$agg_index_sigma[which(x$use_indices == 0)] = 100
x$units_index_paa = c(gb$units_index_paa[-3])

temp = gb$index_paa[-(gb$n_years*2 + 1:gb$n_years),]
x$index_paa = cbind(temp)
x$use_index_paa = cbind(gb$use_index_paa[,-3])
x$index_Neff = cbind(gb$index_Neff[,-3])
x$index_aref = cbind(gb$index_aref[,-3])

x$q_lower = c(gb$q_lower[-3])
x$q_upper = c(gb$q_upper[-3])

x$n_estimated_selpars = gb$n_estimated_selpars - 2
x$n_other_selpars = gb$n_other_selpars
x$other_selpars = c(gb$other_selpars)
x$other_selpars[which(x$other_selpars<1e-15)] = 0.001
x$estimated_selpar_pointers = c(gb$estimated_selpar_pointers[1:8],gb$estimated_selpar_pointers[11:12] - 2)
x$other_selpar_pointers = c(gb$other_selpar_pointers)
x$selpars_lower = c(gb$selpars_lower[-(9:10)])
x$selpars_upper = c(gb$selpars_upper[-(9:10)])

x$n_NAA_sigma = gb$n_NAA_sigma
x$NAA_sigma_pointers =  gb$NAA_sigma_pointers
x$recruit_model = 2 # random about mean
x$use_mat_model = 1
x$use_growth_model = 1

temp = aggregate(as.integer(!(mat.dat$MATURITY == 'I')), by = list(yc = mat.dat$yc, age = mat.dat$AGE), FUN = sum)
temp = cbind.data.frame(temp, aggregate(as.integer(!(mat.dat$MATURITY == 'I')), by = list(yc = mat.dat$yc, age = mat.dat$AGE), FUN = length)[[3]])
names(temp)[3:4] = c("Y","N")
x$Y = temp$Y
x$N = temp$N
x$age_obs = temp$age
x$year_obs = temp$age + temp$yc - 1962
x$Ecov_maxages_k = rep(0,length(x$Y))
x$Ecov_maxages_a50 = rep(0,length(x$Y))
x$Ecov_obs = cbind(fall.temp.dat$anom_bt)
x$use_Ecov_obs = cbind(match(as.integer(!is.na(x$Ecov_obs)),c(0,1)) - 1)
x$Ecov_obs[which(is.na(x$Ecov_obs))] = -999
x$Ecov_obs_sigma = cbind(fall.temp.dat$sd1_bt)
x$Ecov_obs_sigma[which(is.na(x$Ecov_obs_sigma))] = -999
x$fit_k = 0
x$fit_a50 = 0
x$binomial = 0
x$model_years = match(1978:2014, 1963:2014)

#from fits of growth model by itself above
temp = input$dat

x$age_obs_growth = temp$age_obs
x$year_obs_growth = temp$year_obs
x$weight = temp$weight
x$iswt = temp$iswt
x$len = temp$len
x$islen = temp$islen
x$Ecov_maxages_b = temp$Ecov_maxages_b
x$Ecov_maxages_a = temp$Ecov_maxages_a
x$Ecov_maxages_Linf = temp$Ecov_maxages_Linf
x$fit_k_LVB = 0

#This must be consistent with x$Ecov_y used with the maturity data. Here we are using Ecov in the fall of the year that individuals are born
#x$Ecov_pmat = t(sapply(1977 + 1:x$n_years, function(y)  y - 1:x$n_ages - 1962))

y = list(mean_rec_pars = 10)
y$logit_q = rep(-8, x$n_indices)
y$log_F1 = rep(-2, x$n_fleets)
y$F_devs = matrix(0, x$n_years-1, x$n_fleets)
y$log_N1 = rep(10, x$n_ages) #x$N1_model = 0
y$log_NAA_sigma = rep(0, x$n_NAA_sigma)
y$estimated_selpars = rep(0, x$n_estimated_selpars)
n_catch_acomp_pars = c(1,1,1,3,1,2)[x$age_comp_model_fleets[which(apply(x$use_catch_paa,2,sum)>0)]]
n_index_acomp_pars = c(1,1,1,3,1,2)[x$age_comp_model_indices[which(apply(x$use_index_paa,2,sum)>0)]]
y$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
y$index_paa_pars = rep(0, sum(n_index_acomp_pars))
y$log_NAA = matrix(10, x$n_years-1, x$n_ages)

y$beta_k = 0
y$beta_a50 = 0
y$beta_Ecov_k = rep(0,max(1,x$Ecov_maxages_k+1))
y$beta_Ecov_a50 = rep(0,max(1,x$Ecov_maxages_a50+1))
y$Ecov_mu = 0
y$Ecov_AR_pars = rep(0, 2)
y$Ecov_re = cbind(c(x$Ecov_obs))
y$Ecov_re[which(x$use_Ecov_obs == 0)] = 0
y$beta_phi = 0
y$k_AR_pars = rep(0,2)
y$k_re = rep(0, length(y$Ecov_re))
y$a50_AR_pars = rep(0,2)
y$a50_re = rep(0, length(y$Ecov_re))

y$beta_b = log(3)
y$beta_a = -10
y$beta_k_LVB = rep(0, max(x$age_obs))
y$beta_Linf = log(100)
y$t0 = 0
y$beta_Ecov_b = rep(0,max(1,x$Ecov_maxages_b+1))
y$beta_Ecov_a = rep(0,max(1,x$Ecov_maxages_a+1))
y$beta_Ecov_k_LVB = rep(0,max(x$age_obs))
y$beta_Ecov_Linf = rep(0,max(x$age_obs))
#y$beta_Ecov_Linf = rep(0,max(1,x$Ecov_maxages_Linf+1))
y$k_LVB_AR_pars = rep(0,2)
y$k_LVB_re = rep(0, length(y$Ecov_re))
y$beta_sig_Linf = -1
y$beta_sig_L_obs = -10
y$beta_sig_W_obs = -1

gbcod.ss = list(dat = x, par = y)
gbcod.ss$dat$percentSPR = 40.0
gbcod.ss$map = list(
  catch_paa_pars = factor(rep(NA,length(gbcod.ss$par$catch_paa_pars))), 
  index_paa_pars = factor(rep(NA,length(gbcod.ss$par$index_paa_pars))),
  beta_Ecov_k = factor(rep(NA, length(gbcod.ss$par$beta_Ecov_k))),
  beta_Ecov_a50 = factor(rep(NA, length(gbcod.ss$par$beta_Ecov_a50))),
  #beta_phi = factor(rep(NA, length(gbcod.ss$par$beta_phi))),
  k_AR_pars = factor(rep(NA, length(gbcod.ss$par$k_AR_pars))),
  k_re = factor(rep(NA, length(gbcod.ss$par$k_re))),
  a50_AR_pars = factor(rep(NA, length(gbcod.ss$par$a50_AR_pars))),
  a50_re = factor(rep(NA, length(gbcod.ss$par$a50_re))),
  beta_k_LVB = factor(rep(1, length(gbcod.ss$par$beta_k_LVB))),
  beta_Ecov_b = factor(rep(NA, length(gbcod.ss$par$beta_Ecov_b))),
  beta_Ecov_a = factor(rep(NA, length(gbcod.ss$par$beta_Ecov_a))),
  beta_Ecov_k_LVB = factor(rep(NA, length(gbcod.ss$par$beta_Ecov_k_LVB))),
  k_LVB_AR_pars = factor(rep(NA, length(gbcod.ss$par$k_LVB_AR_pars))),
  k_LVB_re = factor(rep(NA, length(gbcod.ss$par$k_LVB_re))),
  beta_Ecov_Linf = factor(rep(NA, length(gbcod.ss$par$beta_Ecov_Linf))),
  beta_sig_L_obs = factor(NA)
  )

library(TMB)
dyn.unload("ssm_env_ssb_v4.so")
compile("ssm_env_ssb_v4.cpp","-O0 -g")
dyn.load("ssm_env_ssb_v4.so")
#gdbsource("ssm_env_ssb_v4.r")

#no temperature effects on maturity or growth, beta-binomial model for maturity
gbcod.ss.mod0 <- MakeADFun(gbcod.ss$dat,gbcod.ss$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re"), 
  map = gbcod.ss$map)
gbcod.ss.mod0$opt <- nlminb(gbcod.ss.mod0$par,gbcod.ss.mod0$fn,gbcod.ss.mod0$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod0$gr(gbcod.ss.mod0$opt$par)))
x = gbcod.ss.mod0
gbcod.ss.mod0$sdrep = sdreport(x)
gbcod.ss.mod0$rep = gbcod.ss.mod0$report()
z = summary(gbcod.ss.mod0$sdrep, "report")
range(z[rownames(z) == "log_SSB_E",2])
range(z[rownames(z) == "log_SSB40_E",2])
range(z[rownames(z) == "log_waa",2])

temp = summary(gbcod.ss.mod0$sdrep)
temp = temp[grep("log_SSB", rownames(temp)),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(1978:2014, exp(temp[1:gbcod.ss$dat$n_years,1]), type = 'n', ylim = c(0, max(exp(temp))), xlab = 'Year', ylab = 'SSB (mt)')
axis(1, lwd = 2)
axis(2, lwd = 2)
box(lwd = 2)
grid(col = gray(0.7), lwd = 2)
lines(1978:2014, exp(temp[1:gbcod.ss$dat$n_years,1]))
lines(1978:2014, exp(temp[gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years,1]), col = 'red')

temp = summary(gbcod.ss.mod0$sdrep)
temp = temp[grep("log_SSB", rownames(temp)),]
temp = temp[gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years,2]/temp[1:gbcod.ss$dat$n_years,2]

plot(1978:2014, 100*(temp-1), type = 'n', ylim = c(0, 100*max(temp-1)), xlab = 'Year', ylab = '% increase in CV')
axis(1, lwd = 2)
axis(2, lwd = 2)
box(lwd = 2)
grid(col = gray(0.7), lwd = 2)
lines(1978:2014, 100*(temp-1))

#constant beta-binomial maturity (BB_1), but growth from growth5 (G_6) fitted in ss_temp_growth_v4.r
temp = gbcod.ss
temp$par = gbcod.ss.mod0$env$parList()
temp$par$beta_k_LVB = growth5$parList$beta_k
temp$par$k_LVB_re = growth5$parList$k_re
temp$par$k_LVB_AR_pars  = growth5$parList$k_AR_pars
temp$par$beta_Ecov_k_LVB = growth5$parList$beta_Ecov_k
temp$par$Ecov_re = growth5$parList$Ecov_re
temp$par$Ecov_AR_pars = growth5$parList$Ecov_AR_pars
temp$par$Ecov_mu = growth5$parList$Ecov_mu
temp$par$beta_b = growth5$parList$beta_b
temp$par$beta_a = growth5$parList$beta_a
temp$par$beta_Linf = growth5$parList$beta_Linf
temp$par$t0 = growth5$parList$t0
temp$par$beta_sig_Linf = growth5$parList$beta_sig_Linf
temp$par$beta_sig_W_obs = growth5$parList$beta_sig_W_obs

temp$dat$fit_k_LVB = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_LVB_re","k_LVB_AR_pars"))]
temp$map$beta_Ecov_k_LVB = factor(c(1, rep(NA,length(temp$par$beta_Ecov_k_LVB)-1)))

gbcod.ss.mod.bb1.g6 <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re", "k_LVB_re"), 
  map = temp$map)
gbcod.ss.mod.bb1.g6$opt <- nlminb(gbcod.ss.mod.bb1.g6$par,gbcod.ss.mod.bb1.g6$fn,gbcod.ss.mod.bb1.g6$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod.bb1.g6$gr(gbcod.ss.mod.bb1.g6$opt$par)))
x = gbcod.ss.mod.bb1.g6
gbcod.ss.mod.bb1.g6$sdrep = sdreport(x)
gbcod.ss.mod.bb1.g6$rep = gbcod.ss.mod.bb1.g6$report()

#constant growth from growth0 (G_1), but maturity from mod4.bb.k.a50.AR (BB_5) fitted in ss_mat_paper_v4.r,
temp = gbcod.ss
temp$par = gbcod.ss.mod0$env$parList()
maxage = 2
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$par$beta_k = mod4.bb.k.a50.AR$parList$beta_k
temp$par$beta_a50 = mod4.bb.k.a50.AR$parList$beta_a50
temp$par$beta_phi = mod4.bb.k.a50.AR$parList$beta_phi
temp$par$k_AR_pars = mod4.bb.k.a50.AR$parList$k_AR_pars
temp$par$k_re = mod4.bb.k.a50.AR$parList$k_re
temp$par$a50_AR_pars = mod4.bb.k.a50.AR$parList$a50_AR_pars
temp$par$a50_re = mod4.bb.k.a50.AR$parList$a50_re

temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
gbcod.ss.mod.bb5.g1 <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re", "k_re", "a50_re"), 
  map = temp$map)
gbcod.ss.mod.bb5.g1$opt <- nlminb(gbcod.ss.mod.bb5.g1$par,gbcod.ss.mod.bb5.g1$fn,gbcod.ss.mod.bb5.g1$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod.bb5.g1$gr(gbcod.ss.mod.bb5.g1$opt$par)))
x = gbcod.ss.mod.bb5.g1
gbcod.ss.mod.bb5.g1$sdrep = sdreport(x)
gbcod.ss.mod.bb5.g1$rep = gbcod.ss.mod.bb5.g1$report()


#maturity from mod4.bb.k.a50.AR (BB_5) fitted in ss_mat_paper_v4.r, growth from growth5 (G_6) fitted in ss_temp_growth_v4.r
temp = gbcod.ss
temp$par = gbcod.ss.mod0$env$parList()
temp$par$beta_k_LVB = growth5$parList$beta_k
temp$par$k_LVB_re = growth5$parList$k_re
temp$par$k_LVB_AR_pars  = growth5$parList$k_AR_pars
temp$par$beta_Ecov_k_LVB = growth5$parList$beta_Ecov_k
temp$par$Ecov_re = growth5$parList$Ecov_re
temp$par$Ecov_AR_pars = growth5$parList$Ecov_AR_pars
temp$par$Ecov_mu = growth5$parList$Ecov_mu
temp$par$beta_b = growth5$parList$beta_b
temp$par$beta_a = growth5$parList$beta_a
temp$par$beta_Linf = growth5$parList$beta_Linf
temp$par$t0 = growth5$parList$t0
temp$par$beta_sig_Linf = growth5$parList$beta_sig_Linf
temp$par$beta_sig_W_obs = growth5$parList$beta_sig_W_obs

temp$dat$fit_k_LVB = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_LVB_re","k_LVB_AR_pars"))]
temp$map$beta_Ecov_k_LVB = factor(c(1, rep(NA,length(temp$par$beta_Ecov_k_LVB)-1)))

#mod4.bb.k.a50.AR from ss_mat_paper_v4.r
maxage = 2
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$par$beta_k = mod4.bb.k.a50.AR$parList$beta_k
temp$par$beta_a50 = mod4.bb.k.a50.AR$parList$beta_a50
temp$par$beta_phi = mod4.bb.k.a50.AR$parList$beta_phi
temp$par$k_AR_pars = mod4.bb.k.a50.AR$parList$k_AR_pars
temp$par$k_re = mod4.bb.k.a50.AR$parList$k_re
temp$par$a50_AR_pars = mod4.bb.k.a50.AR$parList$a50_AR_pars
temp$par$a50_re = mod4.bb.k.a50.AR$parList$a50_re

temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]

gbcod.ss.mod1 <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re", "k_LVB_re", "k_re", "a50_re"), 
  map = temp$map)
gbcod.ss.mod1$opt <- nlminb(gbcod.ss.mod1$par,gbcod.ss.mod1$fn,gbcod.ss.mod1$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod1$gr(gbcod.ss.mod1$opt$par)))
x = gbcod.ss.mod1
gbcod.ss.mod1$sdrep = sdreport(x)
gbcod.ss.mod1$rep = gbcod.ss.mod1$report()

x = summary(gbcod.ss.mod1$sdrep)
x[rownames(x) %in% "a50_re",1]
tail(gbcod.ss.mod1$rep$pmat)
x = summary(mod4.bb.k.a50.AR$sdrep)
x[rownames(x) %in% "a50_re",1]

gbcod.ss.mod1$report()$fix_pmat
gbcod.ss.mod1$env$parList()$a50_re

x <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re", "k_LVB_re", "k_re"), 
  map = temp$map)
x = x$report()$fix_pmat
x = gbcod.ss.mod1$report()$fix_pmat
k = exp(x[1] + x[4])
a50 = exp(x[5] + x[7])
1/(1+exp(-k*(1-a50)))
1/(1+exp(-x[8]))
k = exp(x[1])
a50 = exp(x[2])
1/(1+exp(-k*(1-a50)))
1/(1+exp(-x[9]*(1-x[10])))

tail(mod7.bb.k.a50.AR$rep$pmat)

z = summary(gbcod.ss.mod1$sdrep, "report")
range(z[rownames(z) == "log_SSB_E",2])
range(z[rownames(z) == "log_SSB40_E",2])
range(z[rownames(z) == "log_waa",2])

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[rownames(temp) %in% c("t0","Ecov_y","k_LVB_re","beta_Linf","beta_k_LVB","beta_Ecov_k_LVB"),1]
ind = sum(names(temp) == "Ecov_y")
ind = c(1:(sum(names(temp) == "Ecov_y")-2),ind)
temp = temp[-which(names(temp) == "Ecov_y")[ind]]
temp = temp[-which(names(temp) == "k_LVB_re")[ind]]
temp
k1 = exp(temp["beta_k_LVB"] + temp['beta_Ecov_k_LVB']*temp["Ecov_y"] + temp["k_LVB_re"])
k2 = 
exp(temp["beta_Linf"])*(1-exp(-k1*(1-temp["t0"])))
exp(temp["beta_Linf"])*(1-exp(-0.2421831*(1-temp["t0"])))

#gbcod.ss.mod1 with maturity known
temp = gbcod.ss
temp$dat = gbcod.ss.mod1$env$data
temp$par = gbcod.ss.mod1$env$parList()
temp$map = gbcod.ss.mod1$env$map
temp$map$beta_phi = factor(rep(NA,length(temp$par$beta_phi)))
temp$map$beta_k = factor(rep(NA,length(temp$par$beta_k)))
temp$map$beta_a50 = factor(rep(NA,length(temp$par$beta_a50)))
temp$map$a50_re = factor(rep(NA,length(temp$par$a50_re)))
temp$map$a50_AR_pars = factor(rep(NA,length(temp$par$a50_AR_pars)))
temp$map$k_re = factor(rep(NA,length(temp$par$k_re)))
temp$map$k_AR_pars = factor(rep(NA,length(temp$par$k_AR_pars)))
gbcod.ss.mod1.m_known <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re", "k_LVB_re"), 
  map = temp$map)
gbcod.ss.mod1.m_known$opt <- nlminb(gbcod.ss.mod1.m_known$par,gbcod.ss.mod1.m_known$fn,gbcod.ss.mod1.m_known$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod1.m_known$gr(gbcod.ss.mod1.m_known$opt$par)))
x = gbcod.ss.mod1.m_known
gbcod.ss.mod1.m_known$sdrep = sdreport(x)
gbcod.ss.mod1.m_known$rep = gbcod.ss.mod1.m_known$report()
z = summary(gbcod.ss.mod1.m_known$sdrep, "report")
range(z[rownames(z) == "log_SSB_E",2])
range(z[rownames(z) == "log_SSB40_E",2])
range(z[rownames(z) == "log_waa",2])


#gbcod.ss.mod1 with growth known
temp = gbcod.ss
temp$dat = gbcod.ss.mod1$env$data
temp$par = gbcod.ss.mod1$env$parList()
temp$map = gbcod.ss.mod1$env$map
temp$map$beta_b = factor(rep(NA,length(temp$par$beta_b)))
temp$map$beta_a = factor(rep(NA,length(temp$par$beta_a)))
temp$map$beta_k_LVB = factor(rep(NA,length(temp$par$beta_k_LVB)))
temp$map$beta_Linf = factor(rep(NA,length(temp$par$beta_Linf)))
temp$map$t0 = factor(rep(NA,length(temp$par$t0)))
temp$map$beta_sig_Linf = factor(rep(NA,length(temp$par$beta_sig_Linf)))
temp$map$beta_sig_W_obs = factor(rep(NA,length(temp$par$beta_sig_W_obs)))
temp$map$beta_Ecov_k_LVB = factor(rep(NA,length(temp$par$beta_Ecov_k_LVB)))
temp$map$k_LVB_re = factor(rep(NA,length(temp$par$k_LVB_re)))
temp$map$k_LVB_AR_pars = factor(rep(NA,length(temp$par$k_LVB_AR_pars)))
temp$map$Ecov_re = factor(rep(NA,length(temp$par$Ecov_re)))
temp$map$Ecov_AR_pars = factor(rep(NA,length(temp$par$Ecov_AR_pars)))
temp$map$Ecov_mu = factor(rep(NA,length(temp$par$Ecov_mu)))
gbcod.ss.mod1.g_known = MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "k_re", "a50_re"), 
  map = temp$map)
gbcod.ss.mod1.g_known$opt <- nlminb(gbcod.ss.mod1.g_known$par,gbcod.ss.mod1.g_known$fn,gbcod.ss.mod1.g_known$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod1.g_known$gr(gbcod.ss.mod1.g_known$opt$par)))
x = gbcod.ss.mod1.g_known
gbcod.ss.mod1.g_known$sdrep = sdreport(x)
gbcod.ss.mod1.g_known$rep = gbcod.ss.mod1.g_known$report()
z = summary(gbcod.ss.mod1.g_known$sdrep, "report")
range(z[rownames(z) == "log_SSB_E",2])
range(z[rownames(z) == "log_SSB40_E",2])
range(z[rownames(z) == "log_waa",2])

#gbcod.ss.mod1 with growth and maturity known
temp = gbcod.ss
temp$dat = gbcod.ss.mod1$env$data
temp$par = gbcod.ss.mod1$env$parList()
temp$map = gbcod.ss.mod1$env$map
temp$map$beta_phi = factor(rep(NA,length(temp$par$beta_phi)))
temp$map$beta_k = factor(rep(NA,length(temp$par$beta_k)))
temp$map$beta_a50 = factor(rep(NA,length(temp$par$beta_a50)))
temp$map$a50_re = factor(rep(NA,length(temp$par$a50_re)))
temp$map$a50_AR_pars = factor(rep(NA,length(temp$par$a50_AR_pars)))
temp$map$k_re = factor(rep(NA,length(temp$par$k_re)))
temp$map$k_AR_pars = factor(rep(NA,length(temp$par$k_AR_pars)))
temp$map$beta_b = factor(rep(NA,length(temp$par$beta_b)))
temp$map$beta_a = factor(rep(NA,length(temp$par$beta_a)))
temp$map$beta_k_LVB = factor(rep(NA,length(temp$par$beta_k_LVB)))
temp$map$beta_Linf = factor(rep(NA,length(temp$par$beta_Linf)))
temp$map$t0 = factor(rep(NA,length(temp$par$t0)))
temp$map$beta_sig_Linf = factor(rep(NA,length(temp$par$beta_sig_Linf)))
temp$map$beta_sig_W_obs = factor(rep(NA,length(temp$par$beta_sig_W_obs)))
temp$map$beta_Ecov_k_LVB = factor(rep(NA,length(temp$par$beta_Ecov_k_LVB)))
temp$map$k_LVB_re = factor(rep(NA,length(temp$par$k_LVB_re)))
temp$map$k_LVB_AR_pars = factor(rep(NA,length(temp$par$k_LVB_AR_pars)))
temp$map$Ecov_re = factor(rep(NA,length(temp$par$Ecov_re)))
temp$map$Ecov_AR_pars = factor(rep(NA,length(temp$par$Ecov_AR_pars)))
temp$map$Ecov_mu = factor(rep(NA,length(temp$par$Ecov_mu)))
gbcod.ss.mod1.known = MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA"), 
  map = temp$map)
gbcod.ss.mod1.known$opt <- nlminb(gbcod.ss.mod1.known$par,gbcod.ss.mod1.known$fn,gbcod.ss.mod1.known$gr, control = list(iter.max = 1000, eval.max = 1000))
max(abs(gbcod.ss.mod1.known$gr(gbcod.ss.mod1.known$opt$par)))
x = gbcod.ss.mod1.known
gbcod.ss.mod1.known$sdrep = sdreport(x)
gbcod.ss.mod1.known$rep = gbcod.ss.mod1.known$report()
z = summary(gbcod.ss.mod1.known$sdrep, "report")
range(z[rownames(z) == "log_SSB_E",2])
range(z[rownames(z) == "log_SSB40_E",2])
range(z[rownames(z) == "log_waa",2])

temp = summary(gbcod.mat.2.sdrep)
temp = temp[grep("log_SSB", rownames(temp)),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(1978:2014, exp(temp[1:gbcod$dat$n_years,1]), type = 'n', ylim = c(0, max(exp(temp))), xlab = 'Year', ylab = 'SSB (mt)')
axis(1, lwd = 2)
axis(2, lwd = 2)
box(lwd = 2)
grid(col = gray(0.7), lwd = 2)
lines(1978:2014, exp(temp[1:gbcod$dat$n_years,1]))
lines(1978:2014, exp(temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]), col = 'red')

z = summary(gbcod.ss.mod1$sdrep)
exp(z[rownames(z) == "log_SSB40_E",1])
plot(exp(z[rownames(z) == "log_SSB_E",1]), type = 'l')
lines(exp(z[rownames(z) == "log_SSB",1]), col = 'red')

GBcod.asap.res = read.table("../2015_COD_GB_MOD_asap_RUN6_DIFF_ZERO_2012.rep", skip = 953, nrows = 37)[,3]
GBcod.asap.res = cbind(SSB = GBcod.asap.res, F = read.table("../2015_COD_GB_MOD_asap_RUN6_DIFF_ZERO_2012.rep", skip = 874, nrows = 37)[,2])
GBcod.asap.res = cbind(GBcod.asap.res, R = read.table("../2015_COD_GB_MOD_asap_RUN6_DIFF_ZERO_2012.rep", skip = 913, nrows = 37)[,1])
tcol <- col2rgb('black')
tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')


cairo_pdf('~/work/cod/ss_maturity/tex/SSB_F_R_CJFAS_v2.pdf', family = "Times", height = 10, width = 5)
par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,5,0,0))
temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(1978:2014, temp[,1], type = 'n', axes = FALSE, ylim = c(0,max(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
lines(1978:2014, GBcod.asap.res[,1]/1000, lwd = 2, lty = 2)
lines(1978:2014, exp(temp[,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste("SSB (", 10^3, " mt)")), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[which(rownames(temp) == "log_F"),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ylim = c(0,max(exp(temp))), xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
lines(1978:2014, GBcod.asap.res[,2], lwd = 2, lty = 2)
lines(1978:2014, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(italic(F)), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod1$sdrep)
exp(temp[which(rownames(temp) == "log_N1"),1])
exp(temp[which(rownames(temp) == "log_NAA")[1:5],1])
ind = c(which(rownames(temp) == "log_N1"),t(matrix(which(rownames(temp) == "log_NAA"),36,10)))
x = list(matrix(temp[ind,1], 37, 10, byrow = TRUE))
x[2:3] = list(x[[1]] - matrix(qnorm(0.975)*temp[ind,2], 37, 10, byrow = TRUE),x[[1]] + matrix(qnorm(0.975)*temp[ind,2], 37, 10, byrow = TRUE))
plot(1978:2014, x[[1]][,1], type = 'n', axes = FALSE, ylim = c(0,max(exp(x[[3]][,1])))/1000, xlab = '', ylab = '')
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
lines(1978:2014, GBcod.asap.res[,3]/1000, lwd = 2, lty = 2)
lines(1978:2014, exp(x[[1]][,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(x[[2]][,1],rev(x[[3]][,1])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(italic(R), " (", 10^6,")")), line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

cairo_pdf('~/work/cod/ss_maturity/tex/BRPs_CJFAS_v2.pdf', family = "Times", height = 8, width = 6)
par(mfrow = c(2,1), mar = c(1,1,1,1), oma = c(4,4,0,0))

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[grep("log_SSB40", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
ind = gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years
temp = cbind(temp[ind,1], temp[ind,1] +qnorm(0.975)*cbind(-temp[ind,2],temp[ind,2]))
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = c(65,110))#ylim = c(range(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(1978:2014, exp(temp[,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(SSB[40], " (", 10^3, "mt)")), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[grep("log_F40", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
temp = cbind(temp[ind,1], temp[ind,1] +qnorm(0.975)*cbind(-temp[ind,2],temp[ind,2]))
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = c(0.17,0.23))#ylim = c(range(exp(temp))), xlab = '', ylab = '')
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(1978:2014, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(italic(F)[40]), line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[grep("SPR_0", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
ind = gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years
temp = cbind(log(temp[ind,1]), log(temp[ind,1]) +qnorm(0.975)*cbind(-temp[ind,2],temp[ind,2])/temp[ind,1])
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = range(exp(temp)))#ylim = c(range(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(1978:2014, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(SSB[textstyle("40%")], " (", 10^3, "mt)")), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[grep("SPR40", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
ind = gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years
temp = cbind(log(temp[ind,1]), log(temp[ind,1]) +qnorm(0.975)*cbind(-temp[ind,2],temp[ind,2])/temp[ind,1])
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = range(exp(temp)))#ylim = c(range(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(1978:2014, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(SSB[textstyle("40%")], " (", 10^3, "mt)")), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[grep("YPR40", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
ind = gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years
temp = cbind(log(temp[ind,1]), log(temp[ind,1]) +qnorm(0.975)*cbind(-temp[ind,2],temp[ind,2])/temp[ind,1])
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = range(exp(temp)))#ylim = c(range(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(1978:2014, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(SSB[textstyle("40%")], " (", 10^3, "mt)")), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[grep("log_Y40", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
ind = gbcod.ss$dat$n_years + 1:gbcod.ss$dat$n_years
temp = cbind(temp[ind,1], temp[ind,1] +qnorm(0.975)*cbind(-temp[ind,2],temp[ind,2]))
plot(1978:2014, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = range(exp(temp)))#ylim = c(range(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(1978:2014, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(1978:2014,2014:1978), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(SSB[textstyle("40%")], " (", 10^3, "mt)")), line = 3, cex = 1.5)


z = summary(gbcod.ss.mod1.known$sdrep, "report")
known.prec.res = cbind(z[rownames(z) == "log_SSB_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "SPR_0_E",2]/z[rownames(z) == "SPR_0_E",1])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "SPR40_E",2]/z[rownames(z) == "SPR40_E",1])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "YPR40_E",2]/z[rownames(z) == "YPR40_E",1])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_F40_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_SSB40_E",2])

z = summary(gbcod.ss.mod1.g_known$sdrep, "report")
g_known.prec.res = cbind(z[rownames(z) == "log_SSB_E",2])
g_known.prec.res = cbind(g_known.prec.res, z[rownames(z) == "SPR_0_E",2]/z[rownames(z) == "SPR_0_E",1])
g_known.prec.res = cbind(g_known.prec.res, z[rownames(z) == "SPR40_E",2]/z[rownames(z) == "SPR40_E",1])
g_known.prec.res = cbind(g_known.prec.res, z[rownames(z) == "YPR40_E",2]/z[rownames(z) == "YPR40_E",1])
g_known.prec.res = cbind(g_known.prec.res, z[rownames(z) == "log_F40_E",2])
g_known.prec.res = cbind(g_known.prec.res, z[rownames(z) == "log_SSB40_E",2])

z = summary(gbcod.ss.mod1.m_known$sdrep, "report")
m_known.prec.res = cbind(z[rownames(z) == "log_SSB_E",2])
m_known.prec.res = cbind(m_known.prec.res, z[rownames(z) == "SPR_0_E",2]/z[rownames(z) == "SPR_0_E",1])
m_known.prec.res = cbind(m_known.prec.res, z[rownames(z) == "SPR40_E",2]/z[rownames(z) == "SPR40_E",1])
m_known.prec.res = cbind(m_known.prec.res, z[rownames(z) == "YPR40_E",2]/z[rownames(z) == "YPR40_E",1])
m_known.prec.res = cbind(m_known.prec.res, z[rownames(z) == "log_F40_E",2])
m_known.prec.res = cbind(m_known.prec.res, z[rownames(z) == "log_SSB40_E",2])

z = summary(gbcod.ss.mod1$sdrep, "report")
prec.res = cbind(z[rownames(z) == "log_SSB_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "SPR_0_E",2]/z[rownames(z) == "SPR_0_E",1])
prec.res = cbind(prec.res, z[rownames(z) == "SPR40_E",2]/z[rownames(z) == "SPR40_E",1])
prec.res = cbind(prec.res, z[rownames(z) == "YPR40_E",2]/z[rownames(z) == "YPR40_E",1])
prec.res = cbind(prec.res, z[rownames(z) == "log_F40_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "log_SSB40_E",2])

cairo_pdf('~/work/cod/ss_maturity/tex/SSB_F40_YPR.pdf', family = "Times", height = 10, width = 5)
#png(filename = '~/work/cod/ss_maturity/tex/SSB_F40_YPR.png', width = 10*144, height = 7*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
plot(1978:2014, prec.res[,1], type = "n", ylim = c(0,20), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/known.prec.res)[,1]-1), lwd = 2)
lines(1978:2014, 100*((g_known.prec.res/known.prec.res)[,1]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((m_known.prec.res/known.prec.res)[,1]-1), lty = 3, lwd = 2)
#axis(1, lwd = 2, cex.axis = 1.5)
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
text(2007,18, "SSB", cex = 2)

mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)

plot(1978:2014, prec.res[,6], type = "n", ylim = c(0,20), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/known.prec.res)[,6]-1), lwd = 2)
lines(1978:2014, 100*((g_known.prec.res/known.prec.res)[,6]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((m_known.prec.res/known.prec.res)[,6]-1), lty = 3, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
axis(1, labels = FALSE, lwd = 2)
box(lwd = 2)
text(2007,18, expression(SSB[40]), cex = 2)

plot(1978:2014, prec.res[,5], type = "n", ylim = c(0,20), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/known.prec.res)[,5]-1), lwd = 2)
lines(1978:2014, 100*((g_known.prec.res/known.prec.res)[,5]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((m_known.prec.res/known.prec.res)[,5]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
text(2007,18, expression(italic(F)[40]), cex = 2)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

cairo_pdf('~/work/cod/ss_maturity/tex/SPR100.pdf', family = "Times", height = 5, width = 5)
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
plot(1978:2014, prec.res[,2], type = "n", ylim = c(0,500), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/g_known.prec.res)[,2]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((prec.res/m_known.prec.res)[,2]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
#text(2007,100, expression(SPR[100])), cex = 1.5)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)
dev.off()

plot(1978:2014, prec.res[,4], type = "n", ylim = c(0,20), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/known.prec.res)[,4]-1), lwd = 2)
lines(1978:2014, 100*((g_known.prec.res/known.prec.res)[,4]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((m_known.prec.res/known.prec.res)[,4]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, labels = FALSE, lwd = 2)
box(lwd = 2)
text(2007,18, expression(YPR(italic(F)[40])), cex = 1.5)

par(mfrow = c(1,2), mar = c(1,1,1,1), oma = c(4,4,0,0))
plot(1978:2014, prec.res[,2], type = "n", ylim = c(0,500), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/g_known.prec.res)[,2]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((prec.res/m_known.prec.res)[,2]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
text(2007,100, expression(SPR[100])), cex = 1.5)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)

plot(1978:2014, prec.res[,3], type = "n", ylim = c(0,500), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(1978:2014, 100*((prec.res/g_known.prec.res)[,3]-1), lty = 2, lwd = 2)
lines(1978:2014, 100*((prec.res/m_known.prec.res)[,3]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, labels = FALSE, lwd = 2)
box(lwd = 2)
text(2007,100, expression(SPR(italic(F)[40])), cex = 1.5)

plot(1978:2014, prec.res[,6], type = "n", ylim = c(0,0.2), ylab = "", xlab = "")
lines(1978:2014, (prec.res/known.prec.res)[,6]-1, lwd = 2)
lines(1978:2014, (g_known.prec.res/known.prec.res)[,6]-1, lty = 2, lwd = 2)
lines(1978:2014, (m_known.prec.res/known.prec.res)[,6]-1, lty = 3, lwd = 2)

m_known.prec.res/prec.res
prec.res/known.prec.res

z = summary(gbcod.ss.mod1$sdrep, "report")
temp = z[rownames(z) == "log_SSB_E",2]
z = summary(gbcod.ss.mod1.known$sdrep, "report")
temp = temp/z[rownames(z) == "log_SSB_E",2]

z = summary(gbcod.ss.mod1$sdrep, "report")
temp = z[rownames(z) == "log_SSB40_E",2]
z = summary(gbcod.ss.mod1.known$sdrep, "report")
temp = temp/z[rownames(z) == "log_SSB40_E",2]

z = summary(gbcod.ss.mod1$sdrep, "report")
exp(z[rownames(z) == "log_F40_E",1])
temp = z[rownames(z) == "log_F40_E",2]
z = summary(gbcod.ss.mod1.known$sdrep, "report")
exp(z[rownames(z) == "log_F40_E",1])
temp = temp/z[rownames(z) == "log_F40_E",2]

z = summary(gbcod.ss.mod1$sdrep, "report")
temp = z[rownames(z) == "SPR_0_E",2]/z[rownames(z) == "SPR_0_E",1]
z = summary(gbcod.ss.mod1.known$sdrep, "report")
temp = temp/z[rownames(z) == "SPR_0_E",2]

tcol <- col2rgb('black')
tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')



temp = summary(gbcod.ss.mod0$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
plot(1978:2014, exp(temp[,1]), type = 'l')
temp = summary(gbcod.ss.mod1$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
lines(1978:2014, exp(temp[,1]), col = 'red')

temp.fn = function(ylim. = c(0.7,1.4))
{
  plot(1978:2014, exp(x[,1]), type = 'n', ann = FALSE, axes = FALSE, ylim = ylim.)
  axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  box(lwd = 2)
  grid(col = gray(0.7), lwd = 1, lty = 2)
  #abline(h = 1, lwd = 2)
}

cairo_pdf('~/work/cod/ss_maturity/tex/SSB_BRP_ratio_CJFAS_v2.pdf', family = "Times", height = 10, width = 5)
par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
#par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = summary(gbcod.ss.mod0$sdrep)
x = temp[which(rownames(temp) == "log_SSB_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_SSB_E"),]
temp.fn()
lines(1978:2014, exp(x[,1]), lwd = 2)
temp = summary(gbcod.ss.mod.bb1.g6$sdrep)
x = temp[which(rownames(temp) == "log_SSB_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_SSB_E"),]
lines(1978:2014, exp(x[,1]), lwd = 2, lty = 2)
temp = summary(gbcod.ss.mod.bb5.g1$sdrep)
x = temp[which(rownames(temp) == "log_SSB_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_SSB_E"),]
lines(1978:2014, exp(x[,1]), lwd = 2, lty = 3)
mtext(side = 2, outer = FALSE, "SSB ratio", line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod0$sdrep)
x = temp[which(rownames(temp) == "log_SSB40_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_SSB40_E"),]
temp.fn(c(0.8,1.2))
lines(1978:2014, exp(x[,1]), lwd = 2)
temp = summary(gbcod.ss.mod.bb1.g6$sdrep)
x = temp[which(rownames(temp) == "log_SSB40_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_SSB40_E"),]
lines(1978:2014, exp(x[,1]), lwd = 2, lty = 2)
temp = summary(gbcod.ss.mod.bb5.g1$sdrep)
x = temp[which(rownames(temp) == "log_SSB40_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_SSB40_E"),]
lines(1978:2014, exp(x[,1]), lwd = 2, lty = 3)
mtext(side = 2, outer = FALSE, expression(paste(SSB[40], " ratio")), line = 3, cex = 1.5)

temp = summary(gbcod.ss.mod0$sdrep)
x = temp[which(rownames(temp) == "log_F40_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_F40_E"),]
temp.fn(c(0.9,1.2))
lines(1978:2014, exp(x[,1]), lwd = 2)
temp = summary(gbcod.ss.mod.bb1.g6$sdrep)
x = temp[which(rownames(temp) == "log_F40_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_F40_E"),]
lines(1978:2014, exp(x[,1]), lwd = 2, lty = 2)
temp = summary(gbcod.ss.mod.bb5.g1$sdrep)
x = temp[which(rownames(temp) == "log_F40_E"),]
temp = summary(gbcod.ss.mod1$sdrep)
x = x - temp[which(rownames(temp) == "log_F40_E"),]
lines(1978:2014, exp(x[,1]), lwd = 2, lty = 3)
axis(1, lwd = 2, cex.axis = 1.5)
mtext(side = 2, outer = FALSE, expression(paste(italic(F)[40], " ratio")), line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()


library(TMB)
dyn.load("ssm_env_ssb_v4.so")
peel.fit.fn = function(peel, model = gbcod.ss.mod1)
{
  
  print(peel)
  temp = list(dat = model$env$data, par = model$env$parList(), map = model$env$map)
  temp$dat$n_years = temp$dat$n_years - peel
  temp$dat$Maa = temp$dat$Maa[1:temp$dat$n_years,]
  ind = numeric()
  for(i in 1:temp$dat$n_indices) ind = c(ind, 1:temp$dat$n_years + (i-1)*model$env$data$n_years)
  print(model$env$data$n_years)
  print(ind)
  print(dim(temp$dat$index_paa))
  temp$dat$index_paa = temp$dat$index_paa[ind,]
  ind = numeric()
  for(i in 1:temp$dat$n_fleets) ind = c(ind, 1:temp$dat$n_years + (i-1)*model$env$data$n_years)
  temp$dat$catch_paa = temp$dat$catch_paa[ind,]
  log_NAA_na_ind = rbind(matrix(1:(temp$dat$n_ages*(temp$dat$n_years-1)), temp$dat$n_years-1, temp$dat$n_ages), matrix(rep(NA, peel*temp$dat$n_ages), peel, temp$dat$n_ages))
  F_devs_na_ind = rbind(matrix(1:(temp$dat$n_fleets * (temp$dat$n_years-1)), temp$dat$n_years-1, temp$dat$n_fleets), matrix(rep(NA, peel * temp$dat$n_fleets), peel, temp$dat$n_fleets))
  #Ecov_re_na_ind = c(1:(temp$dat$n_years_Ecov-1), rep(NA, peel))
  #if(model$env$random == "log_R") temp$map$log_R = factor(log_R_na_ind)
  #temp$map$Ecov_re = factor(Ecov_re_na_ind)
  temp$map$log_NAA = factor(log_NAA_na_ind)
  temp$map$F_devs = factor(F_devs_na_ind)
  temp$map$catch_paa_pars = factor(rep(NA,length(temp$par$catch_paa_pars))) 
  temp$map$index_paa_pars = factor(rep(NA,length(temp$par$index_paa_pars)))

  temp.mod <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA", "Ecov_re", "k_LVB_re", "k_re", "a50_re"), 
  map = temp$map)
  temp.opt = nlminb(temp.mod$par,temp.mod$fn,temp.mod$gr, control = list(iter.max = 1000, eval.max = 1000))
  return(list(opt = temp.opt, rep = temp.mod$report()))
}

#temp = peel.fit.fn(1)
#temp.mod <- MakeADFun(temp$dat,temp$par,DLL="ssm_env_ssb_v4", random = c("log_NAA","Ecov_re"), map = temp$map)
temp = list(peel.fit.fn(0))
  
gbcod.ss.mod1$rep$SSB
temp$peels = list(peel.fit.fn(1))
temp$peels[[2]] = peel.fit.fn(2)
temp$peels[[3]] = peel.fit.fn(3)
temp$peels[[4]] = peel.fit.fn(4)
temp$peels[[5]] = peel.fit.fn(5)
gbcod.ss.mod1$peels = temp$peels
gbcod.ss.mod1$peels[[6]] = peel.fit.fn(6)
gbcod.ss.mod1$peels[[7]] = peel.fit.fn(7)
gbcod.ss.mod1$peels[1:5] = temp$peels

mean(sapply(1:5, function(x) temp$peels[[x]]$rep$SSB[gbcod.ss$dat$n_years-x]/gbcod.ss.mod1$rep$SSB[gbcod.ss$dat$n_years-x] - 1))
mean(sapply(1:5, function(x) temp$peels[[x]]$rep$F[gbcod.ss$dat$n_years-x]/gbcod.ss.mod1$rep$F[gbcod.ss$dat$n_years-x] - 1))
mean(sapply(1:7, function(x) gbcod.ss.mod1$peels[[x]]$rep$SSB[gbcod.ss$dat$n_years-x]/gbcod.ss.mod1$rep$SSB[gbcod.ss$dat$n_years-x] - 1))
mean(sapply(1:7, function(x) gbcod.ss.mod1$peels[[x]]$rep$F[gbcod.ss$dat$n_years-x]/gbcod.ss.mod1$rep$F[gbcod.ss$dat$n_years-x] - 1))

x = latex(round(gbcod.mat.2.rep$pmat,2), file = '~/work/cod/ss_maturity/tex/mataa_table_M5.tex', 
  rowlabel = 'Year', rowname = 1978:2014, colheads = paste('Age ', 1:gbcod$n_ages, c(rep('',gbcod$n_ages-1),'+'), sep =''), 
  table.env = FALSE)
x = latex(round(gbcod.mat.2.rep$NAA/1000,2), file = '~/work/cod/ss_maturity/tex/NAA_table_M5.tex', 
  rowlabel = 'Year', rowname = 1978:2014, colheads = paste('Age ', 1:gbcod$n_ages, c(rep('',gbcod$n_ages-1),'+'), sep =''), 
  table.env = FALSE)
x = latex(round(gbcod$waa[gbcod$waa_pointer_fleets,,],2), file = '~/work/cod/ss_maturity/tex/WAA_table_catch.tex', 
  rowlabel = 'Year', rowname = 1978:2014, colheads = paste('Age ', 1:gbcod$n_ages, c(rep('',gbcod$n_ages-1),'+'), sep =''), 
  table.env = FALSE)
x = latex(round(gbcod$waa[gbcod$waa_pointer_ssb,,],2), file = '~/work/cod/ss_maturity/tex/WAA_table_ssb.tex', 
  rowlabel = 'Year', rowname = 1978:2014, colheads = paste('Age ', 1:gbcod$n_ages, c(rep('',gbcod$n_ages-1),'+'), sep =''), 
  table.env = FALSE)

