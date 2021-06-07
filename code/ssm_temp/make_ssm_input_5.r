#same as make_ssm_input_4.r except allow growth model to be estimated
load("data/petrale_inputs.Rds")
load(file = "data/temp_objects_for_ss_lvb_temp.RData")
ldata_years = min(tempdat$year):max(tempdat$year) #1977-2018, but there are gaps in data until 2003.
ecov_year1 = 1969-45 #just to make sure it goes back far enough
modyears = ecov_year1:max(ldata_years)

years = 1969:2018
x = list(n_years = length(years)) 
#x$n_ages = n_a_dat
x$n_ages = 17 #this is the max age in the age comp observations. n_a_dat=40 is for weight at age, maturity at age, etc.
x$n_fleets = 1 #4
#the third GB index only has values in the first 4 years, so not useful for the joint model.
x$n_indices = 2 #CPUE for fleets 1,3 and survey fleet 7
x$n_ages_pop = n_a_dat #40
#x$n_ages_pop = x$n_ages #40

#remove selectivity stuff for 3rd gb index
x$n_selblocks = x$n_fleets + 1 #first two CPUE will use same selectivity as fleets 1,3
x$selblock_models = rep(2,x$n_selblocks) #logistic
x$selblock_pointer_fleets = t(matrix(1:x$n_fleets, x$n_fleets, x$n_years))
x$selblock_pointer_indices = t(matrix(c(1,2), x$n_indices, x$n_years))

x$age_comp_model_fleets = rep(1, x$n_fleets)
x$n_age_comp_pars_fleets = rep(1, x$n_fleets)
x$age_comp_model_indices = rep(1, x$n_indices)
x$n_age_comp_pars_indices = rep(1, x$n_indices)

x$fracyr_SSB = rep(spawning_timing_yr_fraction, x$n_years)
x$mature = t(matrix(maturity_at_age[1 + 1:x$n_ages,3], x$n_ages, x$n_years)) #starts at 0
x$waa_pointer_fleets = 1:x$n_fleets
x$waa_pointer_totcatch = 1 #6
x$waa_pointer_indices = rep(1,x$n_indices)
x$waa_pointer_ssb = 1
x$waa_pointer_jan1 = 1 #doesn't matter

#use average weight at age of ages 17:40 in plus group. males don't get very big. Would imply larger M for males.
weights  = exp(-(16:39)*natM[1]) #equilibrium survive to age 17:40 (no fishing)
weights = weights/sum(weights) #equilibrium proportion alive at age, give survive to age 17
waa = array(NA, dim = c(1,x$n_years, x$n_ages))
for(i in 1){
	temp = as.matrix(subset(waa_matrix, Fleet == i & Sex == 1 & Yr %in% years)[,-(1:7)])
	waa[i,,1:(x$n_ages-1)] = temp[,1:(x$n_ages-1)]
	waa[i,,x$n_ages] = apply(temp[,(x$n_ages):NCOL(temp)],1,function(x) sum(x*weights))
	temp = as.matrix(subset(waa_matrix, Fleet == i & Sex == 2 & Yr %in% years)[,-(1:7)])
	waa[i,,1:(x$n_ages-1)] = (waa[i,,1:(x$n_ages-1)] + temp[,1:(x$n_ages-1)])/2
	waa[i,,x$n_ages] = (waa[i,,x$n_ages] + apply(temp[,(x$n_ages):NCOL(temp)],1,function(x) sum(x*weights)))/2
}
x$waa = waa
x$Maa = matrix(natM[1], x$n_years, x$n_ages)

x$agg_catch = as.matrix(cbind(subset(catch_matrix, year %in% years)[,-1]))
x$agg_catch = cbind(apply(x$agg_catch,1,sum))
x$agg_catch_sigma = as.matrix(cbind(subset(catch_se_matrix, year %in% years)[,2]))

source("code/ssm_temp/get_aref_fn.r")
x$catch_paa = matrix(NA, x$n_years*x$n_fleets, x$n_ages)
x$catch_aref = matrix(NA, x$n_years, x$n_fleets)
x$use_catch_paa = x$catch_Neff = matrix(0, x$n_years, x$n_fleets)
temp = list()
for(i in 1:4) {
	ind = age_comp_catch$rowind[,2] == i & age_comp_catch$rowind[,1] %in% years
	tyrs = age_comp_catch$rowind[ind,1]
	temp[[i]] = age_comp_catch$age_comp[ind,]
	temp[[i]] = t(apply(temp[[i]], 1, function(x) x/sum(x)))
	temp[[i]] = temp[[i]][match(years,tyrs),]
}
for(i in 1:x$n_years) for(j in 1:x$n_ages) x$catch_paa[i,j] = mean(sapply(temp, function(x) x[i,j]), na.rm = TRUE) 
x$catch_aref = cbind(get_aref_fn(x$catch_paa))
x$use_catch_paa[] = 1 
temp = matrix(0,x$n_years, 4)
for(i in 1:4) {
	ind = age_comp_catch$rowind[,2] == i & age_comp_catch$rowind[,1] %in% years
	tyrs = age_comp_catch$rowind[ind,1]
	temp[match(tyrs, years, nomatch= 0),i] = age_comp_catch$Nsamp[match(tyrs,age_comp_catch$Nsamp[,1]),i+1] 
}
x$catch_Neff = cbind(apply(temp,1,sum))

x$units_indices = rep(1,x$n_indices)
x$fracyr_indices = t(matrix(catch_timing_as_fraction[c(1,5)],x$n_indices,x$n_years)) #I think this is right?
temp = CPUE_data
temp = temp[match(years, temp$year),]
x$agg_indices = as.matrix(cbind(temp[,1+c(1,5)]))
temp = matrix(1, NROW(x$agg_indices), NCOL(x$agg_indices))
temp[is.na(x$agg_indices)] = 0
x$use_indices = temp
temp = CPUE_se_log
temp = temp[match(years, temp$year),]
x$agg_index_sigma = as.matrix(cbind(temp[,1+c(1,5)]))
x$units_index_paa = rep(2, x$n_indices)

x$index_paa = matrix(NA, x$n_years*x$n_indices, x$n_ages)
temp = as.matrix(age_comp_survey$age_comp[match(years,age_comp_survey$Nsamp$Yr),])
temp = t(apply(temp, 1, function(x) x/sum(x)))

x$index_paa[x$n_years + 1:x$n_years,] = temp
temp = matrix(0, x$n_years, x$n_indices)
temp[match(age_comp_survey$Nsamp$Yr,years),2] = 1
x$use_index_paa = temp
temp = matrix(0, x$n_years, x$n_indices)
temp[match(age_comp_survey$Nsamp$Yr,years),2] = age_comp_survey$Nsamp[,2]
x$index_Neff = temp

temp = as.matrix(age_comp_survey$age_comp)
temp = t(apply(temp, 1, function(x) x/sum(x)))
temp = get_aref_fn(temp)
x$index_aref = matrix(NA, x$n_years, x$n_indices)
x$index_aref[,2] = temp[match(years,age_comp_survey$Nsamp$Yr)]

x$q_lower = rep(0, x$n_indices)
x$q_upper = rep(1000, x$n_indices)

x$n_estimated_selpars = 2 * x$n_selblocks
x$n_other_selpars = 0
x$other_selpars = numeric(0)
x$estimated_selpar_pointers = 1:x$n_estimated_selpars
x$other_selpar_pointers = numeric(0)
x$selpars_lower = rep(0, x$n_estimated_selpars)
x$selpars_upper = rep(x$n_ages, x$n_estimated_selpars)

x$n_NAA_sigma = 2
x$NAA_sigma_pointers =  c(1,rep(2,x$n_ages_pop-1))
x$recruit_model = 2 # random about mean

x$use_growth_model = 1
x$percentSPR = 40

if(!exists("m2_long")) m2_long = readRDS(file = "results/ssm_temp/m2_long.RDS")


dat = m2_long$env$data
x$Ecov_obs = dat$Ecov_obs
x$use_Ecov_obs = dat$use_Ecov_obs
#x$Ecov_obs[which(is.na(x$Ecov_obs))] = -999
x$Ecov_obs_sigma = dat$Ecov_obs_sigma
#x$Ecov_obs_sigma[which(is.na(x$Ecov_obs_sigma))] = -999

#which of the years in the environmental time series (Which may be longer) are the years of the assessment model
x$model_years = match(years, modyears)

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
ssm_input = list(dat = x)



y = list(mean_rec_pars = 10)
y$logit_q = c(-8,log(1/(1000-1)))#rep(-8, x$n_indices)
y$log_F1 = rep(-3, x$n_fleets)
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

if(!exists("input")) input = readRDS(file = "code/ssm_temp/growth_input.RDS")

temp = input
z = c(NROW(temp$par$beta_Ecov_k), NCOL(temp$par$beta_Ecov_k))
temp$map$beta_Ecov_k = matrix(NA, z[1],z[2])
temp$map$beta_Ecov_k[1,1] = 1
temp$map$beta_Ecov_k = factor(temp$map$beta_Ecov_k)
temp$dat$fit_k = 1
temp$map = temp$map[-which(names(temp$map) %in% c("k_re","k_AR_pars"))]
temp$par = m2_long$parList

ssm_input$par = y
ssm_input$map = list()
ssm_input$map$index_paa_pars = as.factor(rep(NA, length(ssm_input$par$index_paa_pars)))
ssm_input$map$catch_paa_pars = as.factor(rep(NA, length(ssm_input$par$catch_paa_pars)))
ssm_input$map$beta_k_LVB = temp$map$beta_k
ssm_input$map$beta_Ecov_k_LVB = temp$map$beta_Ecov_k
map.pars = c("beta_Ecov_b", "beta_Ecov_a", "log_Ecov_obs_sig_scale", "beta_Ecov_Linf", "beta_sig_L_obs")
ssm_input$map[map.pars] = temp$map[map.pars]
#map.pars = names(m2_long$parList)
#map.pars[c(3,8,13,14)] = c("beta_k_LVB","beta_Ecov_k_LVB","k_LVB_AR_pars","k_LVB_re")
#for(i in map.pars) ssm_input$map[[i]] = as.factor(rep(NA,length(ssm_input$par[[i]])))

#ssm_input$random = c("log_NAA", "N1_re")#, "Ecov_re", "k_LVB_re")
ssm_input$random = c("log_NAA", "N1_re", "Ecov_re", "k_LVB_re")
# ssm_input$par$log_N1 = 5 #random walk in the rest of the ages
# ssm_input$par$N1_re[] = 5 #random walk in the rest of the ages
# ssm_input$par$log_NAA[] = 5
ssm_input$dat$do_growth = 1
#ssm_input$dat$fit_k_LVB = 0

saveRDS(ssm_input, file = "code/ssm_temp/ssm_input.RDS")
