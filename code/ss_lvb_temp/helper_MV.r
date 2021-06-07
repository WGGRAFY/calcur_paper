library(plotrix)

fit.tmb.fn = function(model, n.newton=3)
{
  model$opt <- nlminb(model$par, model$fn,model$gr)
  if(n.newton) for(i in 1:n.newton) { # Take a few extra newton steps 
    g <- as.numeric(model$gr(model$opt$par))
    h <- optimHess(model$opt$par, model$fn, model$gr)
    model$opt$par <- model$opt$par - solve(h, g)
    model$opt$objective <- model$fn(model$opt$par)
  }
  model$rep <- model$report()
  if(is.null(model$env$random))  model$parList = model$env$parList(model$env$last.par)
  else model$parList = model$env$parList()

  model$sdrep <- try(sdreport(model,model$opt$par))
  return(model)
}

make_input = function(dat, edat)
{
  #fall.temp.dat = read.csv("1963-2014-fall-temp-for-state-space.csv", header = TRUE)
  x = list(weight = dat$weight_kg)
  x$iswt = as.integer(x$weight > 0 & !is.na(x$weight))
  x$len = dat$length_cm
  x$islen = as.integer(x$len > 0 & !is.na(x$len))
  x$age_obs = dat$age_years
  x$year_obs = x$age_obs + dat$yc - (min(modyears) - 1) #have to pay attention to this for when in the year the environmental covariate is affecting growth
  x$cy_obs = dat$yc - min(modyears) #have to pay attention to this for when in the year the environmental covariate is affecting growth
  x$Ecov_maxages_b = rep(0,length(x$len))
  x$Ecov_maxages_a = rep(0,length(x$len))
  x$Ecov_maxages_Linf = rep(0,length(x$len))
  #no covariates in likelihood yet
  #No covariates until we figure out what to use.
  if(missing(edat)) x$Ecov_obs = cbind(rep(0, nyears)) #cbind(fall.temp.dat$anom_bt)
  else x$Ecov_obs = edat[[1]] 
  x$use_Ecov_obs = matrix(NA, NROW(x$Ecov_obs), NCOL(x$Ecov_obs))
  x$use_Ecov_obs[] = match(as.integer(!is.na(x$Ecov_obs)),c(0,1)) - 1
  #let's see if NAs are allowed
  #x$Ecov_obs[which(is.na(x$Ecov_obs))] = -999
  if(missing(edat)) x$Ecov_obs_sigma = cbind(rep(1, length(x$Ecov_obs)))
  else x$Ecov_obs_sigma = edat[[2]]
  #x$Ecov_obs_sigma[which(is.na(x$Ecov_obs_sigma))] = -999
  if(missing(edat)) x$fit_Ecov = 0
  else x$fit_Ecov = 1
  if(missing(edat)) x$fit_Ecov_obs = 0
  else x$fit_Ecov_obs = 1
  x$fit_growth = 1
  x$fit_k = 0
  x$laa_matrix_max_age = max(x$age_obs, na.rm = TRUE)
  x$X_Linf = cbind(rep(1,length(x$len)))
  x$obs_for_Linf_y = 1
  
  y = list(beta_b = log(3))
  y$beta_a = -10
  y$beta_k = rep(0, max(x$age_obs))
  y$beta_Linf = log(100)
  y$t0 = 0
  y$beta_Ecov_b = rep(0,max(1,x$Ecov_maxages_b+1))
  y$beta_Ecov_b = matrix(y$beta_Ecov_b, length(y$beta_Ecov_b), NCOL(x$Ecov_obs))
  y$beta_Ecov_a = rep(0,max(1,x$Ecov_maxages_a+1))
  y$beta_Ecov_a = matrix(y$beta_Ecov_a, length(y$beta_Ecov_a), NCOL(x$Ecov_obs))
  y$beta_Ecov_k = rep(0,max(x$age_obs))
  y$beta_Ecov_k = matrix(y$beta_Ecov_k, length(y$beta_Ecov_k), NCOL(x$Ecov_obs))
  y$beta_Ecov_Linf = rep(0,max(x$age_obs))
  y$beta_Ecov_Linf = matrix(y$beta_Ecov_Linf, length(y$beta_Ecov_Linf), NCOL(x$Ecov_obs))
  y$Ecov_mu = rep(0, NCOL(x$Ecov_obs))
  y$Ecov_AR_pars = rep(0, 2)
  y$Ecov_AR_pars = matrix(y$Ecov_AR_pars, length(y$Ecov_AR_pars), NCOL(x$Ecov_obs))
  y$Ecov_re = x$Ecov_obs
  y$Ecov_re[which(x$use_Ecov_obs == 0)] = 0
  y$k_AR_pars = rep(0,2)
  y$k_re = rep(0, NROW(y$Ecov_re))
  y$log_Ecov_obs_sig_scale = rep(0, NCOL(x$Ecov_obs))
  y$beta_sig_Linf = -1
  y$beta_sig_L_obs = -10
  y$beta_sig_W_obs = -1
  out = list(dat = x, par = y)
  out$map = list(
    beta_k = factor(rep(1, length(out$par$beta_k))),
    beta_Ecov_b = factor(rep(NA, length(out$par$beta_Ecov_b))),
    beta_Ecov_a = factor(rep(NA, length(out$par$beta_Ecov_a))),
    beta_Ecov_k = factor(rep(NA, length(out$par$beta_Ecov_k))),
    k_AR_pars = factor(rep(NA, length(out$par$k_AR_pars))),
    k_re = factor(rep(NA, length(out$par$k_re))),
    log_Ecov_obs_sig_scale = factor(rep(NA, length(out$par$log_Ecov_obs_sig_scale))),
    beta_Ecov_Linf = factor(rep(NA, length(out$par$beta_Ecov_Linf))),
    beta_sig_L_obs = factor(NA)
  )
  if(x$fit_Ecov != 1){
    map$Ecov_mu = factor(rep(NA, length(out$par$Ecov_mu)))
    map$Ecov_AR_pars = factor(rep(NA, length(out$par$Ecov_AR_pars)))
    map$Ecov_re = factor(rep(NA, length(out$par$Ecov_re)))
  }
  return(out)
}

get.logkaa.fn = function(mod)
{
  temp = list(summary(mod$sdrep)[rownames(summary(mod$sdrep)) == "log_k_ya",])
  temp[[2]] = array(dim = c(3,nyears,max(mod$env$data$age_obs)))
  temp[[2]][1,,] = temp[[1]][,1]
  temp[[2]][2,,] = temp[[1]][,1] -qnorm(0.975)*temp[[1]][,2] 
  temp[[2]][3,,] = temp[[1]][,1] +qnorm(0.975)*temp[[1]][,2]
  temp = temp[[2]]
  return(temp)
} 

get.ymax.kaa.fn = function(ages, mod, years)
{
  ind = 1:length(modyears)
  if(!missing(years)) ind = which(modyears %in% years)
  temp = get.logkaa.fn(mod)
  #nyears = length(ind)
  tind = ind
  ymax = 0
  print(ind)
  print(ages)
  #ymax = max(exp(temp[2:3,tind, ages[1]]))
  for(i in ages) 
  {
    if(i %in% ind) tind = ind[ind > i]
    print(tind)
    ymax = max(ymax,exp(temp[2:3,tind,i]))
  }
  return(ymax)
}

plot.kaa.fn = function(a=1, mod = m3, ymax, do.xlabs = TRUE, do.ylabs = TRUE, do.nobs = FALSE, xlim, years)
{
  ind = 1:length(modyears)
  if(!missing(years)) ind = which(modyears %in% years)
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  yrs = modyears[-(1:a)]
  print(yrs)
  print(modyears)
  temp = get.logkaa.fn(mod)
  if(missing(ymax)) ymax = exp(max(temp[2:3,,2:9], na.rm = TRUE))
  temp = exp(temp[,(a+1):nyears,])
  if(missing(xlim)) xlim = range(modyears[ind])
  plot(yrs, temp[1,,a], type = 'n', axes = FALSE, xlab = "", ylab = "", ylim = c(0,ymax), xlim = xlim)
  #plot(modyears[-(1:a)], temp[1,,a], type = 'b', ylim = range(c(0,temp[2:3,,a])), xlim = range(modyears), xlab = "Year", ylab = "k")
  grid(col = gray(0.7), lwd = 1)
  if(do.xlabs) axis(1, lwd = 2, cex.axis = 1.5)
  else axis(1, labels = FALSE, lwd = 2)
  if(do.ylabs) axis(2, lwd = 2, cex.axis = 1.5)
  else axis(2, labels = FALSE, lwd = 2)
  box(lwd = 2)
  abline(h = exp(mod$parList$beta_k[a]))
  lines(yrs, temp[1,,a], lwd = 2)
  polygon(c(yrs,rev(yrs)), c(temp[2,,a], rev(temp[3,,a])), col = tcol, border = "transparent", lty = 2)
  yearborn = mod$env$data$cy_obs + min(modyears)
  yearobs = mod$env$data$age_obs + yearborn
  if(do.nobs) text(modyears,ymax*0.05, sapply(modyears, function(x) sum(yearobs == x)), srt = 90, adj = 1, cex = 0.8) 
}

get.ymax.laa.fn = function(ages, mod = m1, years)
{
  ind = 1:length(modyears)
  if(!missing(years)) ind = modyears %in% years
  dat = mod$env$data
  temp = summary(mod$sdrep)
  temp = temp[which(rownames(temp) == "log_laa"),]
  x = array(dim = c(3, length(modyears),max(dat$age_obs)))
  x[1,,] = temp[,1]
  x[2,,] = temp[,1]-qnorm(0.975)*temp[,2]
  x[3,,] = temp[,1]+qnorm(0.975)*temp[,2]
  ymax = max(exp(x[3,ind,ages]), na.rm = TRUE)
  for(age in ages)
  {
    x = t(sapply(modyears[ind], function(x) 
    {
      y = which(dat$year_obs == x - min(modyears) + 1 &  dat$age_obs == age)
      y = c(mean(dat$len[y]), sd(dat$len[y])/sqrt(length(y)))
      if(y[1] == 0 | is.na(y[2])) y = rep(NA,3)
      else y = c(y[1], exp(log(y[1]) + qnorm(0.975)*c(-1,1)*y[2]/y[1]))
      return(y)
    }))
    ymax = max(ymax, x[,3], na.rm = TRUE)
  }
  return(ymax)
}

get.laa.fn = function(age = 1, mod = m1)
{
  dat = mod$env$data
  temp = summary(mod$sdrep)
  temp = temp[which(rownames(temp) == "log_laa"),]
  x = array(dim = c(3, length(modyears),max(dat$age_obs)))
  x[1,,] = temp[,1]
  x[2,,] = temp[,1]-qnorm(0.975)*temp[,2]
  x[3,,] = temp[,1]+qnorm(0.975)*temp[,2]
  out = list(pred = exp(t(x[,,age])))
  x = t(sapply(modyears, function(x) 
  {
    y = which(dat$year_obs == x - min(modyears) + 1 &  dat$age_obs == age)
    y = c(mean(dat$len[y]), sd(dat$len[y])/sqrt(length(y)))
    if(y[1] == 0 | is.na(y[2])) y = rep(NA,3)
    else y = c(y[1], exp(log(y[1]) + qnorm(0.975)*c(-1,1)*y[2]/y[1]))
    return(y)
  }))
  out$obs = x
  return(out)
}

plot.laa.fn = function(x, ymax = 80,mod = m1, do.xlabs = TRUE, do.ylabs = TRUE, years)
{
  ind = 1:length(modyears)
  if(!missing(years)) ind = modyears %in% years
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  yrs = modyears[-(1:x)]
  res = get.laa.fn(x,mod)
  plot(yrs, res$pred[-(1:x),1], type = 'n', axes = FALSE, xlab = "", ylab = "", ylim = c(0,ymax), xlim = range(modyears[ind]))
  grid(col = gray(0.7), lwd = 1)
  if(do.xlabs) axis(1, lwd = 2, cex.axis = 1.5)
  else axis(1, labels = FALSE, lwd = 2)
  if(do.ylabs) axis(2, lwd = 2, cex.axis = 1.5)
  else axis(2, labels = FALSE, lwd = 2)
  box(lwd = 2)
  lines(yrs, res$pred[-(1:x),1], lwd = 2)
  polygon(c(yrs,rev(yrs)), c(res$pred[-(1:x),2],rev(res$pred[-(1:x),3])), col = tcol, border = "transparent", lty = 2)
  plotCI(yrs, res$obs[-(1:x),1], li = res$obs[-(1:x),2], ui = res$obs[-(1:x),3], add = TRUE, lwd = 2, slty = 2, sfrac = 0)
}

plot.laa.cohort.fn = function(cohort, ymax = 80,mod = m1, do.xlabs = TRUE, do.ylabs = TRUE)
{
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  ages = 1:20
  if(max(modyears)-cohort == 0) stop(paste0("cohort: ", cohort, " is at the end of the time series: ", max(modyears)))
  if(max(modyears)-cohort < 20) ages = 1:(max(modyears)-cohort)
  pred = t(sapply(ages, function(y) get.laa.fn(y,mod)$pred[cohort-modyears[1]+1+y,]))
  obs = t(sapply(ages, function(y) get.laa.fn(y,mod)$obs[cohort-modyears[1]+1+y,]))
  print(pred)
  print(obs)
  plot(ages, pred[,1], type = 'n', axes = FALSE, xlab = "", ylab = "", ylim = c(0,ymax), xlim = c(0,max(ages)))
  grid(col = gray(0.7), lwd = 1)
  if(do.xlabs) axis(1, lwd = 2, cex.axis = 1.5)
  else axis(1, labels = FALSE, lwd = 2)
  if(do.ylabs) axis(2, lwd = 2, cex.axis = 1.5)
  else axis(2, labels = FALSE, lwd = 2)
  box(lwd = 2)
  lines(ages, pred[,1], lwd = 2)
  polygon(c(ages,rev(ages)), c(pred[,2],rev(pred[,3])), col = tcol, border = "transparent", lty = 2)
  plotCI(ages, obs[,1], li = obs[,2], ui = obs[,3], add = TRUE, lwd = 2, slty = 2, sfrac = 0)
}

plot.laa.Ecov.fn = function(x, ymax = 80,mod = m1, do.xlabs = TRUE, do.ylabs = TRUE)
{
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  yrs = modyears[-(1:x)]
  Ecov = mod$rep$Ecov_y[-(1:x)]
  ord = order(Ecov)
  Ecov = Ecov[ord]
  res = get.laa.fn(x,mod)
  res$pred = res$pred[-(1:x),][ord,]
  res$obs = res$obs[-(1:x),][ord,]
  plot(Ecov, res$pred[,1], type = 'n', axes = FALSE, xlab = "", ylab = "", ylim = range(c(res$obs,res$pred),na.rm =TRUE))
  grid(col = gray(0.7), lwd = 1)
  if(do.xlabs) axis(1, lwd = 2, cex.axis = 1.5)
  else axis(1, labels = FALSE, lwd = 2)
  if(do.ylabs) axis(2, lwd = 2, cex.axis = 1.5)
  else axis(2, labels = FALSE, lwd = 2)
  box(lwd = 2)
  lines(Ecov, res$pred[,1], lwd = 2)
  polygon(c(Ecov,rev(Ecov)), c(res$pred[,2],rev(res$pred[,3])), col = tcol, border = "transparent", lty = 2)
  plotCI(Ecov, res$obs[,1], li = res$obs[,2], ui = res$obs[,3], add = TRUE, lwd = 2, slty = 2, sfrac = 0)
}

plot.Ecov.res.fn = function(model, yrange, years, ylab = "Bottom temperature anomaly", ylab2, xrange)
{
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  if(missing(xrange)) xrange = range(years)
  require(plotrix)
  n_Ecov = NCOL(model$env$data$Ecov_obs)
  print(n_Ecov)
  if(missing(years)) years = 1:NROW(model$env$data$Ecov_obs)
  Ecov = array(NA, dim = c(3, length(years), n_Ecov))
  temp = summary(model$sdrep)
  temp = temp[which(rownames(temp) == "Ecov_y"),]
  Ecov[1,,] = temp[,1]
  Ecov[2,,] = temp[,1] - qnorm(0.975)*temp[,2]
  Ecov[3,,] = temp[,1] + qnorm(0.975)*temp[,2]
  par(mfrow = c(n_Ecov,1), mar = c(1,1,1,3), oma = c(4,4,0,0))
  for(p in 1:n_Ecov) {
    yrange = range(Ecov[,,p], na.rm = TRUE)
    plot(years, Ecov[1,,p], type = 'n', axes = FALSE, xlab = "", ylab = "", ylim = yrange, xlim = xrange)
    grid(col = gray(0.7), lwd = 1)
    lines(years, Ecov[1,,p], lwd = 2)
    polygon(c(years,rev(years)), c(Ecov[2,,p],rev(Ecov[3,,p])), col = tcol, border = "transparent", lty = 2)
    obs.ind <- which(model$env$data$use_Ecov_obs[,p] == 1)
    plotCI(years[obs.ind], model$env$data$Ecov_obs[obs.ind,p], 
      li = (model$env$data$Ecov_obs - qnorm(0.975)*model$env$data$Ecov_obs_sigma)[obs.ind,p], 
      ui = (model$env$data$Ecov_obs + qnorm(0.975)*model$env$data$Ecov_obs_sigma)[obs.ind,p], add = TRUE, lwd = 2)
    if(p == n_Ecov) axis(1, lwd = 2, cex.axis = 1.5)
    else axis(1, labels = FALSE, lwd = 2)
    axis(2, lwd = 2, cex.axis = 1.5)
    box(lwd = 2)
    mtext(side = 4, outer = FALSE, line = 1, ylab2[p], cex = 1)
  }
  mtext(side = 2, outer = TRUE, line = 2, ylab, cex = 1.5)
  mtext(side = 1, outer = TRUE, line = 2, "Year", cex = 1.5)
}

mypalette = function(n, alpha=1){
  cols = c("dodgerblue","green","red")
  x = apply(sapply(cols, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
  colorRampPalette(x, alpha = TRUE)(n)
}
