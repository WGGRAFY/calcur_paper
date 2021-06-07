library(TMB)
#source("code/ssm_temp/make_ssm_input_1.r")
#source("code/ssm_temp/make_ssm_input_2.r")

ssm_input = readRDS("code/ssm_temp/ssm_input.RDS")
setwd("code/ssm_temp")
compile("ssm_temp.cpp", "-O0 -g")
setwd("../..")
dyn.load(dynlib(paste0("code/ssm_temp/ssm_temp")))



x = ssm_input
x$par$log_NAA_sigma[2] = log(0.05)
x$dat$catch_Neff = x$dat$catch_Neff*100
#x$map$estimated_selpars = factor(rep(NA, length(x$par$estimated_selpars)))
x$map$logit_q = factor(c(1,2,NA))#rep(NA, length(x$par$logit_q)))
#x$map$log_N1 = factor(rep(NA, length(x$par$log_N1)))
#x$map$mean_rec_pars = factor(rep(NA, length(x$par$mean_rec_pars)))
#x$map$log_F1 = factor(rep(NA, length(x$par$log_F1)))
#x$map$F_devs = factor(rep(NA, length(x$par$F_devs)))
x$map$log_NAA_sigma = factor(c(1,NA))
#x$map$log_sig_N1 = factor(rep(NA, length(x$par$log_sig_N1)))

source("code/ssm_temp/fit_tmb.r")
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp = fit_tmb(y, do.sdrep = FALSE)

tsel = matrix(0, x$dat$n_selblocks, x$dat$n_ages)
apply(temp$rep$selblocks,1,function(x) which(x == max(x))) # 7,7,7,6,9
#all catch_paa for first 3 fleets are 0 for first age, only 2 non-zero for 4th fleet
tsel[1:3,c(2:6,8:17)] = 1
tsel[4,c(2:5,7:17)] = 1
tsel[5,c(1:8,10:17)] = 1
tsel = t(tsel)
x$dat$n_estimated_selpars = sum(tsel)
x$dat$n_other_selpars = length(tsel)-sum(tsel)
x$dat$other_selpars = c(0,1,0,1,0,1,0,1,1)
x$dat$estimated_selpar_pointers = which(tsel == 1)
x$dat$other_selpar_pointers = which(tsel == 0)
x$dat$selpars_lower = rep(0, x$dat$n_estimated_selpars)
x$dat$selpars_upper = rep(1, x$dat$n_estimated_selpars)
x$par$estimated_selpars = rep(0, x$dat$n_estimated_selpars)

tpar = temp$parList
tpar = tpar[names(tpar) != "estimated_selpars"]
x$par[names(tpar)] = tpar

y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp = fit_tmb(y, do.sdrep = FALSE)
y$par = temp$opt$par
temp2 = fit_tmb(y, do.sdrep = FALSE)

tsel = matrix(0, x$dat$n_selblocks, x$dat$n_ages)
apply(temp$rep$selblocks,1,function(x) which(x == max(x))) # 7,7,7,6,9
#all catch_paa for first 3 fleets are 0 for first age, only 2 non-zero for 4th fleet
tsel[1:3,c(2:6,8:17)] = 1
tsel[4,c(2:5,7:17)] = 1
tsel[5,c(1:4,10:17)] = 1
tsel = t(tsel)
x$dat$n_estimated_selpars = sum(tsel)
x$dat$n_other_selpars = length(tsel)-sum(tsel)
x$dat$other_selpars = c(0,1,0,1,0,1,0,1,1,1,1,1,1)
x$dat$estimated_selpar_pointers = which(tsel == 1)
x$dat$other_selpar_pointers = which(tsel == 0)
x$dat$selpars_lower = rep(0, x$dat$n_estimated_selpars)
x$dat$selpars_upper = rep(1, x$dat$n_estimated_selpars)
x$par$estimated_selpars = rep(0, x$dat$n_estimated_selpars)

tpar = temp2$parList
tpar = tpar[names(tpar) != "estimated_selpars"]
x$par[names(tpar)] = tpar

y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp3 = fit_tmb(y, do.sdrep = FALSE)

x$par = temp3$parList
x$map$log_F1 = factor(rep(NA,length(x$par$log_F1)))
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp4 = fit_tmb(y, do.sdrep = FALSE)
temp4$sdrep = sdreport(temp4)
saveRDS(temp4, file = "code/ssm_temp/best_fit_so_far.RDS")

x$map$log_NAA_sigma = factor(c(1,NA,3))
x$par = temp4$parList
x$par$log_NAA_sigma = x$par$log_NAA_sigma[c(1,2,2)]
x$dat$NAA_sigma_pointers = c(1,2,2,rep(3,14))

y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp5 = fit_tmb(y, do.sdrep = FALSE)
#h <- stats::optimHess(temp5$opt$par, temp5$fn, temp5$gr)
Gr = temp5$gr(temp5$opt$par)
df <- data.frame(param = names(temp5$opt$par),
                 MLE = temp5$opt$par,
                 gr.at.MLE = Gr)
ind.hi <- which(Gr > 0.01)
temp5$badpar <- df[ind.hi,]

source("code/ssm_temp/make_ssm_input_3.r")
ssm_input = readRDS("code/ssm_temp/ssm_input.RDS")
source("code/ssm_temp/fit_tmb.r")
x = ssm_input
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp = fit_tmb(y, do.sdrep = FALSE)
temp$sdrep = sdreport(temp)
x = summary(temp$sdrep)
plot(years, exp(x[rownames(x) == "log_F40",1]), type = 'l', ylim = c(0,0.2))
lines(years, exp(x[rownames(x) == "log_F40_E",1]), col = 'red')

plot(years, exp(x[rownames(x) == "log_SSB",1]), type = 'l')
lines(years, exp(x[rownames(x) == "log_SSB_E",1]), col = 'red')

plot(years, exp(x[rownames(x) == "log_SSB40",1]), type = 'l')
lines(years, exp(x[rownames(x) == "log_SSB40_E",1]), col = 'red')


unique(rownames(x))

h = TMBhelper::check_estimability(temp4)
h <- stats::optimHess(temp4$opt$par, temp4$fn, temp4$gr)
v = solve(h)

check <- TMBhelper::check_estimability(y)
if(length(test$WhichBad) > 0){
  bad.par <- as.character(test$BadParams$Param[test$BadParams$Param_check=='Bad'])
  bad.par.grep <- grep(bad.par, test$BadParams$Param)
  model$badpar <- test$BadParams[bad.par.grep,]
  warning(paste("","Some fixed effect parameter(s) are not identifiable.",
    "Consider 1) removing them from the model by fixing input$par and input$map = NA, or",
    "2) changing your model configuration.","",
    paste(capture.output(print(test$BadParams[bad.par.grep,])), collapse = "\n"), sep="\n"))    
}


