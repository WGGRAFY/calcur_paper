library(TMB)
dyn.load(dynlib(paste0("code/ssm_temp/ssm_temp")))
#gdbsource("code/ssm_temp/temp.r", TRUE)
ssm_input = readRDS("code/ssm_temp/ssm_input.RDS")

x = ssm_input
x$dat$recruit_model = 3
x$par$mean_rec_pars = c(0,0)
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
