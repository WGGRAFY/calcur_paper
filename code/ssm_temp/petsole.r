library(TMB)
#source("code/ssm_temp/make_ssm_input_1.r")
#source("code/ssm_temp/make_ssm_input_2.r")
source("code/ssm_temp/make_ssm_input_3.r")

ssm_input = readRDS("code/ssm_temp/ssm_input.RDS")
setwd("code/ssm_temp")
compile("ssm_temp.cpp", "-O0 -g")
setwd("../..")
dyn.load(dynlib(paste0("code/ssm_temp/ssm_temp")))

source("code/ssm_temp/fit_tmb.r")
x = ssm_input
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp = fit_tmb(y, do.sdrep = FALSE)
temp$sdrep = sdreport(temp)
saveRDS(temp, file="code/ssm_temp/best_fit_so_far.RDS")
less_ages = readRDS(file="code/ssm_temp/best_fit_so_far.RDS")
x = summary(less$sdrep)
plot(years, exp(x[rownames(x) == "log_F40",1]), type = 'l', ylim = c(0,0.2))
lines(years, exp(x[rownames(x) == "log_F40_E",1]), col = 'red')

plot(years, exp(x[rownames(x) == "log_SSB",1]), type = 'l')
lines(years, exp(x[rownames(x) == "log_SSB_E",1]), col = 'red')

plot(years, exp(x[rownames(x) == "log_SSB40",1]), type = 'l')
lines(years, exp(x[rownames(x) == "log_SSB40_E",1]), col = 'red')


source("code/ssm_temp/make_ssm_input_4.r")
ssm_input = readRDS("code/ssm_temp/ssm_input.RDS")

x = ssm_input
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
more_ages = fit_tmb(y, do.sdrep = FALSE)
more_ages$sdrep = sdreport(more_ages)
saveRDS(more_ages, file="code/ssm_temp/more_ages.RDS")
x = summary(more_ages$sdrep)
plot(years, exp(x[rownames(x) == "log_F40",1]), type = 'l', ylim = c(0,0.2))
lines(years, exp(x[rownames(x) == "log_F40_E",1]), col = 'red')

plot(years, exp(x[rownames(x) == "log_SSB",1]), type = 'l')
lines(years, exp(x[rownames(x) == "log_SSB_E",1]), col = 'red')

plot(years, exp(x[rownames(x) == "log_SSB40",1]), type = 'l')
lines(years, exp(x[rownames(x) == "log_SSB40_E",1]), col = 'red')

#B-H implies steepnesses ~ 1
x = ssm_input
x$dat$recruit_model = 3
x$par$mean_rec_pars = c(0,0)
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
more_ages_sr = fit_tmb(y, do.sdrep = FALSE)
saveRDS(more_ages_sr, file = "code/ssm_temp/more_ages_sr.RDS")
more_ages_sr = readRDS(file = "code/ssm_temp/more_ages_sr.RDS")
x$par = more_ages_sr$parList
y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
temp = fit_tmb(y, do.sdrep = FALSE)

#Now include growth model estimation
source("code/ssm_temp/make_ssm_input_5.r")
ssm_input = readRDS("code/ssm_temp/ssm_input.RDS")
x = ssm_input
x$par = more_ages$parList

y <- MakeADFun(x$dat,x$par,DLL="ssm_temp", random = x$random, 
  map = x$map)
more_ages_growth_on = fit_tmb(y, do.sdrep = FALSE)
more_ages_growth_on$sdrep = sdreport(more_ages_growth_on)
saveRDS(more_ages_growth_on, file="code/ssm_temp/more_ages_growth_on.RDS")
#more_ages_growth_on = readRDS(file="code/ssm_temp/more_ages_growth_on.RDS")

#years = 1969:2018
tcol <- col2rgb('black')
tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')

#cairo_pdf('results/ssm_temp/weight_age_10.pdf', family = "Times", height = 10, width = 10)
png(filename = "results/ssm_temp/weight_age_10.png", width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = summary(more_ages_growth_on$sdrep)
temp = temp[which(rownames(temp) == "log_waa"),]
temp = list(matrix(temp[,1],ssm_input$dat$n_years, ssm_input$dat$n_ages_pop),matrix(temp[,2],ssm_input$dat$n_years, ssm_input$dat$n_ages_pop))
temp = cbind(temp[[1]][,10],temp[[2]][,10])
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(years, temp[,1], type = 'n', axes = FALSE, ylim = c(0,max(exp(temp))), xlab = '', ylab = '')
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
lines(years, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
lines(years, ssm_input$dat$waa[1,,10], lwd = 2)
polygon(c(years,rev(years)), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, "Mass (kg)", line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
text(max(years), max(exp(temp)), "Age 10", adj = c(1,1), cex = 2)
dev.off()

png(filename = "results/ssm_temp/petrale_SSB.png", width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = summary(more_ages_growth_on$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(years, temp[,1], type = 'n', axes = FALSE, ylim = c(0,max(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
#lines(years, GBcod.asap.res[,1]/1000, lwd = 2, lty = 2)
lines(years, exp(temp[,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(years,rev(years)), exp(c(temp[,2],rev(temp[,3])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste("SSB (", 10^3, " mt)")), line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

cairo_pdf('results/ssm_temp/petrale_SSB_F_R.pdf', family = "Times", height = 10, width = 5)
par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,5,0,0))
temp = summary(more_ages_growth_on$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(years, temp[,1], type = 'n', axes = FALSE, ylim = c(0,max(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
#lines(years, GBcod.asap.res[,1]/1000, lwd = 2, lty = 2)
lines(years, exp(temp[,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(years,rev(years)), exp(c(temp[,2],rev(temp[,3])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste("SSB (", 10^3, " mt)")), line = 3, cex = 1.5)

temp = summary(more_ages_growth_on$sdrep)
temp = temp[which(rownames(temp) == "log_F"),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(years, exp(temp[,1]), type = 'n', axes = FALSE, ylim = c(0,max(exp(temp))), xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
#lines(years, GBcod.asap.res[,2], lwd = 2, lty = 2)
lines(years, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(years,rev(years)), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(italic(F)), line = 3, cex = 1.5)

dat = ssm_input$dat
temp = summary(more_ages_growth_on$sdrep)
exp(temp[which(rownames(temp) == "log_N1"),1])
exp(temp[which(rownames(temp) == "N1_re"),1])
exp(temp[which(rownames(temp) == "log_NAA")[1:5],1])
ind = c(which(rownames(temp) == "log_N1"),which(rownames(temp) == "N1_re"),t(matrix(which(rownames(temp) == "log_NAA"),dat$n_years-1, dat$n_ages_pop)))
x = list(matrix(temp[ind,1], dat$n_years, dat$n_ages_pop, byrow = TRUE))
x[2:3] = list(x[[1]] - matrix(qnorm(0.975)*temp[ind,2], dat$n_years, dat$n_ages_pop, byrow = TRUE),x[[1]] + matrix(qnorm(0.975)*temp[ind,2], dat$n_years, dat$n_ages_pop, byrow = TRUE))
plot(years, x[[1]][,1], type = 'n', axes = FALSE, ylim = c(0,max(exp(x[[3]][,1])))/1000, xlab = '', ylab = '')
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lwd = 1, lty = 2)
#lines(years, GBcod.asap.res[,3]/1000, lwd = 2, lty = 2)
lines(years, exp(x[[1]][,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(years,rev(years)), exp(c(x[[2]][,1],rev(x[[3]][,1])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(italic(R), " (", 10^6,")")), line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

cairo_pdf('results/ssm_temp/petrale_BRPs.pdf', family = "Times", height = 8, width = 6)
par(mfrow = c(2,1), mar = c(1,1,1,1), oma = c(4,4,0,0))

temp = summary(more_ages_growth_on$sdrep)
temp = temp[grep("log_SSB40_E", rownames(temp)),]
temp = cbind(temp[,1], temp[,1] +qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(years, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = c(range(exp(temp)))/1000, xlab = '', ylab = '')
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(years, exp(temp[,1])/1000, lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(years,rev(years)), exp(c(temp[,2],rev(temp[,3])))/1000, col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(paste(SSB[40], " (", 10^3, "mt)")), line = 3, cex = 1.5)

temp = summary(more_ages_growth_on$sdrep)
temp = temp[grep("log_F40_E", rownames(temp)),]
#temp = temp[gbcod$dat$n_years + 1:gbcod$dat$n_years,1]/temp[1:gbcod$dat$n_years,1] # = 1
temp = cbind(temp[,1], temp[,1] +qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(years, exp(temp[,1]), type = 'n', axes = FALSE, ann = FALSE, ylim = c(range(exp(temp))), xlab = '', ylab = '')
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7), lty = 2)
lines(years, exp(temp[,1]), lwd = 2, type = 'b', pch = 19, cex = 1.5)
polygon(c(years,rev(years)), exp(c(temp[,2],rev(temp[,3]))), col = tcol, border = "transparent")
mtext(side = 2, outer = FALSE, expression(italic(F)[40]), line = 3, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

more_ages = readRDS(file="code/ssm_temp/more_ages.RDS")
z = summary(more_ages$sdrep, "report")
known.prec.res = cbind(z[rownames(z) == "log_SSB_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_SPR_0_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_SPR40_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_YPR40_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_F40_E",2])
known.prec.res = cbind(known.prec.res, z[rownames(z) == "log_SSB40_E",2])

z = summary(more_ages_growth_on$sdrep, "report")
prec.res = cbind(z[rownames(z) == "log_SSB_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "log_SPR_0_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "log_SPR40_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "log_YPR40_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "log_F40_E",2])
prec.res = cbind(prec.res, z[rownames(z) == "log_SSB40_E",2])


#cairo_pdf('results/ssm_temp/petrale_SSB_SSB40_F40_CV_ratio.pdf', family = "Times", height = 10, width = 5)
png(filename = 'results/ssm_temp/petrale_SSB_CV_ratio.png', width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
#par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = 100*((prec.res/known.prec.res)[,1]-1)
plot(years, temp, type = "n", ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, temp, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
#axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
text(max(years), max(temp), "SSB", adj = c(1,1), cex = 2)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

png(filename = 'results/ssm_temp/petrale_SSB40_CV_ratio.png', width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
#par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = 100*((prec.res/known.prec.res)[,6]-1)
plot(years, temp, type = "n", ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, temp, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
axis(1, lwd = 2, cex.axis = 1.5)
#axis(1, labels = FALSE, lwd = 2)
box(lwd = 2)
text(max(years), max(temp), expression(SSB[40]), adj = c(1,1), cex = 2)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

png(filename = 'results/ssm_temp/petrale_F40_CV_ratio.png', width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = 100*((prec.res/known.prec.res)[,5]-1)
plot(years, temp, type = "n", ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, temp, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
text(max(years), max(temp), expression(italic(F)[40]), adj = c(1,1), cex = 2)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

cairo_pdf('results/ssm_temp/petrale_SPR0_CV_ratio.pdf', family = "Times", height = 5, width = 5)
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
plot(years, 100*((prec.res/known.prec.res)[,2]-1), type = "n", ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, 100*((prec.res/known.prec.res)[,2]-1), lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
#text(2007,100, expression(SPR[100])), cex = 1.5)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)
dev.off()

plot(years, prec.res[,4], type = "n", ylim = c(0,20), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, 100*((prec.res/known.prec.res)[,4]-1), lwd = 2)
lines(years, 100*((g_known.prec.res/known.prec.res)[,4]-1), lty = 2, lwd = 2)
lines(years, 100*((m_known.prec.res/known.prec.res)[,4]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, labels = FALSE, lwd = 2)
box(lwd = 2)
text(2007,18, expression(YPR(italic(F)[40])), cex = 1.5)

par(mfrow = c(1,2), mar = c(1,1,1,1), oma = c(4,4,0,0))
plot(years, prec.res[,2], type = "n", ylim = c(0,500), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, 100*((prec.res/g_known.prec.res)[,2]-1), lty = 2, lwd = 2)
lines(years, 100*((prec.res/m_known.prec.res)[,2]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
text(2007,100, expression(SPR[100])), cex = 1.5)
mtext(side = 2, outer = TRUE, line = 2, "% increase in CV", cex = 1.5)

plot(years, prec.res[,3], type = "n", ylim = c(0,500), ylab = "", xlab = "", axes = FALSE)
grid(col = gray(0.7), lwd = 1)
lines(years, 100*((prec.res/g_known.prec.res)[,3]-1), lty = 2, lwd = 2)
lines(years, 100*((prec.res/m_known.prec.res)[,3]-1), lty = 3, lwd = 2)
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, labels = FALSE, lwd = 2)
box(lwd = 2)
text(2007,100, expression(SPR(italic(F)[40])), cex = 1.5)

plot(years, prec.res[,6], type = "n", ylim = c(0,0.2), ylab = "", xlab = "")
lines(years, (prec.res/known.prec.res)[,6]-1, lwd = 2)
lines(years, (g_known.prec.res/known.prec.res)[,6]-1, lty = 2, lwd = 2)
lines(years, (m_known.prec.res/known.prec.res)[,6]-1, lty = 3, lwd = 2)

m_known.prec.res/prec.res
prec.res/known.prec.res

z = summary(more_ages_growth_on$sdrep, "report")
temp = z[rownames(z) == "log_SSB_E",2]
z = summary(more_ages$sdrep, "report")
temp = temp/z[rownames(z) == "log_SSB_E",2]

z = summary(more_ages_growth_on$sdrep, "report")
temp = z[rownames(z) == "log_SSB40_E",2]
z = summary(more_ages$sdrep, "report")
temp = temp/z[rownames(z) == "log_SSB40_E",2]

z = summary(more_ages_growth_on$sdrep, "report")
exp(z[rownames(z) == "log_F40_E",1])
temp = z[rownames(z) == "log_F40_E",2]
z = summary(more_ages$sdrep, "report")
exp(z[rownames(z) == "log_F40_E",1])
temp = temp/z[rownames(z) == "log_F40_E",2]

z = summary(more_ages_growth_on$sdrep, "report")
temp = z[rownames(z) == "SPR_0_E",2]/z[rownames(z) == "SPR_0_E",1]
z = summary(more_ages$sdrep, "report")
temp = temp/z[rownames(z) == "SPR_0_E",2]

tcol <- col2rgb('black')
tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')



temp = summary(gbcod.ss.mod0$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
plot(years, exp(temp[,1]), type = 'l')
temp = summary(more_ages_growth_on$sdrep)
temp = temp[which(rownames(temp) == "log_SSB_E"),]
lines(years, exp(temp[,1]), col = 'red')

temp.fn = function(ylim. = c(0.7,1.4), labels = FALSE)
{
  plot(years, exp(x[,1]), type = 'n', ann = FALSE, axes = FALSE, ylim = ylim.)
  axis(1, labels = labels, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  box(lwd = 2)
  grid(col = gray(0.7), lwd = 1, lty = 2)
  #abline(h = 1, lwd = 2)
}

#cairo_pdf('results/ssm_temp/petrale_SSB_BRP_ratio.pdf', family = "Times", height = 10, width = 5)
png(filename = 'results/ssm_temp/petrale_SSB_ratio.png', width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
#par(mfrow = c(3,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = summary(more_ages_growth_on$sdrep)
x = temp[which(rownames(temp) == "log_SSB_E"),]
x = x - temp[which(rownames(temp) == "log_SSB"),]
temp.fn(c(0.8,1.25), labels = TRUE)
lines(years, exp(x[,1]), lwd = 2)
abline(h=1, lwd = 2, col = "red")
text(max(years), 1.25, "SSB", adj = c(1,1), cex = 2)
mtext(side = 2, "Ratio (Temperature/Constant)", outer = TRUE, line = 2, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off()

png(filename = 'results/ssm_temp/petrale_SSB40_ratio.png', width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = summary(more_ages_growth_on$sdrep)
x = temp[which(rownames(temp) == "log_SSB40_E"),]
x = x - temp[which(rownames(temp) == "log_SSB40"),]
temp.fn(c(0.8,1.25), labels = TRUE)
lines(years, exp(x[,1]), lwd = 2)
abline(h=1, lwd = 2, col = "red")
text(max(years), 1.25, expression(SSB[40]), adj = c(1,1), cex = 2)
mtext(side = 2, "Ratio (Temperature/Constant)", outer = TRUE, line = 2, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
dev.off() #mtext(side = 2, outer = FALSE, expression(paste(SSB[40], " ratio")), line = 3, cex = 1.5)

png(filename = 'results/ssm_temp/petrale_F40_ratio.png', width = 5*144, height = 5*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp = summary(more_ages_growth_on$sdrep)
x = temp[which(rownames(temp) == "log_F40_E"),]
x = x - temp[which(rownames(temp) == "log_F40"),]
temp.fn(c(0.8,1.25), labels = TRUE)
lines(years, exp(x[,1]), lwd = 2)
abline(h=1, lwd = 2, col = "red")
text(max(years), 1.25, expression(italic(F)[40]), adj = c(1,1), cex = 2)
mtext(side = 2, "Ratio (Temperature/Constant)", outer = TRUE, line = 2, cex = 1.5)
mtext(side = 1, outer = FALSE, "Year", line = 3, cex = 1.5)
#mtext(side = 2, outer = FALSE, expression(paste(italic(F)[40], " ratio")), line = 3, cex = 1.5)
dev.off()

############ stop here

library(TMB)
dyn.load("ssm_env_ssb_v4.so")
peel.fit.fn = function(peel, model = more_ages_growth_on)
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
  
more_ages_growth_on$rep$SSB
temp$peels = list(peel.fit.fn(1))
temp$peels[[2]] = peel.fit.fn(2)
temp$peels[[3]] = peel.fit.fn(3)
temp$peels[[4]] = peel.fit.fn(4)
temp$peels[[5]] = peel.fit.fn(5)
more_ages_growth_on$peels = temp$peels
more_ages_growth_on$peels[[6]] = peel.fit.fn(6)
more_ages_growth_on$peels[[7]] = peel.fit.fn(7)
more_ages_growth_on$peels[1:5] = temp$peels

mean(sapply(1:5, function(x) temp$peels[[x]]$rep$SSB[gbcod.ss$dat$n_years-x]/more_ages_growth_on$rep$SSB[gbcod.ss$dat$n_years-x] - 1))
mean(sapply(1:5, function(x) temp$peels[[x]]$rep$F[gbcod.ss$dat$n_years-x]/more_ages_growth_on$rep$F[gbcod.ss$dat$n_years-x] - 1))
mean(sapply(1:7, function(x) more_ages_growth_on$peels[[x]]$rep$SSB[gbcod.ss$dat$n_years-x]/more_ages_growth_on$rep$SSB[gbcod.ss$dat$n_years-x] - 1))
mean(sapply(1:7, function(x) more_ages_growth_on$peels[[x]]$rep$F[gbcod.ss$dat$n_years-x]/more_ages_growth_on$rep$F[gbcod.ss$dat$n_years-x] - 1))

x = latex(round(gbcod$waa[gbcod$waa_pointer_fleets,,],2), file = 'results/WAA_table_catch.tex', 
  rowlabel = 'Year', rowname = years, colheads = paste('Age ', 1:gbcod$n_ages, c(rep('',gbcod$n_ages-1),'+'), sep =''), 
  table.env = FALSE)
x = latex(round(gbcod$waa[gbcod$waa_pointer_ssb,,],2), file = 'results/WAA_table_ssb.tex', 
  rowlabel = 'Year', rowname = years, colheads = paste('Age ', 1:gbcod$n_ages, c(rep('',gbcod$n_ages-1),'+'), sep =''), 
  table.env = FALSE)

