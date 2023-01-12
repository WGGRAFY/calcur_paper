library(TMB)
dll_name = "ss_lvb_temp"
source("code/ss_lvb_temp/helper_MV.r")
dyn.load(dynlib(paste0("code/ss_lvb_temp/",dll_name)))
dat = readRDS("data/all_length_data.RDS")
modyears = min(dat$year):max(dat$year) 
nyears = length(modyears)

#get tempdat, regions, tdat.obs, tdat.sd
load(file = "data/temp_objects_for_ss_lvb_temp.RData")

#tempdat = read.csv("paul_temp_data_2020.csv", as.is = TRUE)
comnames = sort(unique(tempdat$species))
sps = c("dbrockfish", "lingcod", "phake", "psanddab", "petsole", "sablefish", "sbrockfish")

temp = cbind.data.frame(comname = comnames)
temp$first_age = 0
temp$area = NA
temp$depth = NA
for(i in 1:length(comnames)) 
{
  sp.tdat = subset(tempdat, species == comnames[i] & age < 6)
  temp$first_age[i] = min(sp.tdat$age)
  temp$area[i] = as.character(unique(sp.tdat$area[sp.tdat$age == temp$first_age[i]]))
  temp$depth[i] = as.character(unique(sp.tdat$depth_bin[sp.tdat$age == temp$first_age[i]]))
}
write.csv(temp, file = "results/ss_lvb_temp/first_age_area_depth.csv")

aic.tab = sapply(sps, function(x) 
  {
  print(x)
  read.csv(paste0("results/ss_lvb_temp/", x, "_aic.csv"), as.is=TRUE, header = TRUE)[,2]
})
aic.tab = apply(aic.tab,2, function(x) x - min(x))
rownames(aic.tab) = paste0("m",0:3)
write.csv(round(aic.tab, 2), file = "results/ss_lvb_temp/aic.csv")

best = sapply(sps, function(x) rownames(aic.tab)[aic.tab[,x] == 0])

temp = t(sapply(1:length(best), function(x) 
    {
      print(sps[x])
      sd = read.csv(paste0("results/ss_lvb_temp/", sps[x], "_", best[x], "_sdrep.csv"), as.is=TRUE)
      return(unlist(sd[sd[,1] == "beta_Ecov_k",2:3]))
  }))
rownames(temp) = comnames
colnames(temp) = c("beta_Ecov_k", "se")
write.csv(round(temp,2), file = "results/ss_lvb_temp/beta_Ecov_k_estimates.csv")

for(i in 1:length(sps)) {
  cairo_pdf(paste0('results/ss_lvb_temp/',sps[i], '_age_0_k_temperature.pdf'), family = "Times", height = 10, width = 10)
    temp = readRDS(paste0("results/ss_lvb_temp/", sps[i], "_", best[i], ".RDS"))
    ind = which(temp$parList$beta_Ecov_k[1,] != 0)
    plot(temp$rep$Ecov_y[,ind], exp(temp$rep$log_k_Ecov[,1]), xlim = c(-1,1), type = "n", axes = FALSE, ann = FALSE)
    grid(col = grey(0.7))
    axis(2, lwd = 2)
    axis(1, lwd = 2)
    box(lwd = 2)
    points(temp$rep$Ecov_y[,ind], exp(temp$rep$log_k_Ecov[,1]), pch = 19)
    #print(cbind(temp$rep$log_k_Ecov[,1], temp$opt$par["beta_k"] + temp$opt$par["beta_Ecov_k"] * temp$rep$Ecov_y[,ind]))
    text(-1, max(exp(temp$rep$log_k_Ecov[,1])), comnames[i], adj = c(0,1), cex = 1.5)
    mtext(side = 1, "Bottom Temperature Anomaly", line = 3)
    mtext(side = 2, "LVB k", line = 3)
  dev.off()
}

par(mfcol = c(4,2), oma = c(4,2,1,1), mar = c(1,3,1,1))
for(i in 1:length(sps)) {
  cairo_pdf(paste0('results/ss_lvb_temp/',sps[i], '_age_1_L_temperature.pdf'), family = "Times", height = 10, width = 10)
    temp = readRDS(paste0("results/ss_lvb_temp/", sps[i], "_", best[i], ".RDS"))
    L1 = exp(temp$parList$beta_Linf[1])*(1 - exp(-exp(temp$rep$log_k_Ecov[,1])*(1-temp$parList$t0)))
    dat = temp$env$data
    obs.L1 = t(sapply(modyears, function(x) 
    {
      y = which(dat$year_obs == x - min(modyears) + 1 &  dat$age_obs == 1)
      y = c(mean(dat$len[y]), sd(dat$len[y])/sqrt(length(y)))
      if(y[1] == 0 | is.na(y[2])) y = rep(NA,3)
      else y = c(y[1], exp(log(y[1]) + qnorm(0.975)*c(-1,1)*y[2]/y[1]))
      return(y)
    }))
    print(obs.L1)
    ind = which(temp$parList$beta_Ecov_k[1,] != 0)
    plot(temp$rep$Ecov_y[,ind], L1, xlim = c(-1,1), type = "n", axes = FALSE, ann = FALSE)
    #plot(m3$rep$Ecov_y, L1, xlim = c(-1,1), type = "n", axes = FALSE, ylim = range(c(L1,obs.L1[,1]), na.rm = TRUE))
    grid(col = grey(0.7))
    axis(2, lwd = 2)
    axis(1, lwd = 2)
    box(lwd = 2)
    points(temp$rep$Ecov_y[,ind], L1, pch = 19, col = "dodgerblue")
    #points(m3$rep$Ecov_y, obs.L1[,1], pch = 19)
    text(-1, max(L1), comnames[i], adj = c(0,1), cex = 1.5)
    mtext(side = 1, "Bottom Temperature Anomaly", line = 3)
    mtext(side = 2, "Length (cm) at age 1", line = 3)
  dev.off()
}
x = matrix(NA, 5, 4)
for(i in which(best == "m3")) {
  print(i)
  temp = readRDS(paste0("results/ss_lvb_temp/", sps[i], "_", best[i], ".RDS"))
  x[which(best == "m3")==i,] <- c(exp(temp$parList$beta_Linf), as.list(temp$sdrep, "Std. Error")$beta_Linf)
}
x[,3:4] = x[,1:2]*x[,3:4]
x = round(x,2)
x = cbind(paste0(x[,1], " (", x[,3], ")"), paste0(x[,2], " (", x[,4], ")"))
colnames(x) = c("Linf<2000", "Linf>1999")
rownames(x) = comnames[which(best == "m3")]
write.csv(x, file = "results/ss_lvb_temp/Linf_2000.csv")

n.tab = sapply(sps, function(x) 
  {
  print(x)
  temp = readRDS(paste0("results/ss_lvb_temp/", x, "_m0.RDS"))
  length(temp$env$data$islen)
})
names(n.tab) = comnames
write.csv(rbind(n.tab), file = "results/ss_lvb_temp/n_by_species.csv")

#get the number of fixed effects parameters
sapply(paste0("m",0:3), function(x) 
  {
  print(x)
  temp = readRDS(paste0("results/ss_lvb_temp/", sps[1], "_", x,".RDS"))
  length(temp$opt$par)
})

tcol <- col2rgb('black')
tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')

#png(filename = "calcur_results/temp_effect_by_species.png", width = 9*144, height = 12*144, res = 144, pointsize = 12, family = "Times")#,
cairo_pdf('results/ss_lvb_temp/temp_effect_by_species.pdf', family = "Times", height = 12, width = 9)
par(mfcol = c(4,2), oma = c(2,2,1,1), mar = c(3,3,1,1))
for(i in 1:length(sps)) {
    temp = readRDS(paste0("results/ss_lvb_temp/", sps[i], "_", best[i], ".RDS"))
    ind = which(!is.na(temp$env$data$Ecov_obs))
    predt = range(temp$env$data$Ecov_obs[ind[which(ind %in% temp$env$data$cy_obs)]])
    print(predt)
    predt = seq(predt[1],predt[2],0.01)
    X = cbind(1,predt)
    ind = c(which(rownames(temp$sdrep$cov.fixed) == "beta_k"), which(rownames(temp$sdrep$cov.fixed) == "beta_Ecov_k"))
    print(ind)
    cov = temp$sdrep$cov.fixed[ind,ind]
    beta = c(temp$parList$beta_k[1], temp$parList$beta_Ecov_k[1,1])
    print(beta)
    k = X %*% beta #predt * temp$parList$beta_Ecov_k[1,1]
    k = cbind(k, sqrt(diag(X %*% cov %*% t(X))))
    #kmult = cbind(kmult, abs(predt)*as.list(temp$sdrep, "Std. Error")$beta_Ecov_k[1,1])
    k = exp(cbind(k[,1], k[,1] - qnorm(0.975)*k[,2],k[,1] + qnorm(0.975)*k[,2]))
    plot(predt, k[,1], type = "n", axes = FALSE, ann = FALSE, ylim = range(k))
    grid(col = grey(0.7))
    axis(2, lwd = 2)
    axis(1, lwd = 2)
    box(lwd = 2)
    lines(predt, k[,1], lwd = 2)
    polygon(c(predt,rev(predt)), c(k[,2],rev(k[,3])), col = tcol, border = "transparent")
    legend("topleft", legend = comnames[i], adj = c(0,1), cex = 1.5, bty = "n")
    #legend("topleft", legend = comnames[i], adj = c(0,1), cex = 1.5, box.col = "white", bg = "white")
}
mtext(side = 1, "Bottom Temperature Anomaly", line = 0, outer = TRUE, cex = 1.5)
mtext(side = 2, "LVB k", line = 0, outer = TRUE, cex = 1.5)
dev.off()

