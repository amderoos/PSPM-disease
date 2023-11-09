envvars <- ls()
rm(list = envvars[(envvars != "defplotpars") & 
                    (envvars != "Q075_MUI_Beta") & (envvars != "Q075_TI_Beta") & 
                    (envvars != "Q125_MUI_Beta") & (envvars != "Q125_TI_Beta")  & 
                    (envvars != "TI010_Beta_Qa") & (envvars != "TI010_Beta_Qb") & 
                    (envvars != "Q075_Beta_TIa") & (envvars != "Q075_Beta_TIb")])
rm("envvars")
basedir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../")
datadir <- paste0(basedir, "Figures/EBToutput/")
fname <- paste0(basedir, "Figures/Figure2.pdf")
setwd(basedir)

ToPdf <- T

if (ToPdf) pdf(file = fname, width = 7.0, height = 8.0)

library(rootSolve)
library(latex2exp)
parnames <- c("Rho", "Rmax", "Sb", "Sm", "M", "Q", "H", "SIGMA", 
              "TS", "TI", "TIJ", "TIA",
              "MUS", "MUI", "MUIJ", "MUIA", "BETA")
defpars <- c(0.1, 100.0, 0.1, 1.0, 1.0, 1.0, 3.0, 0.5, 
             0.1, 0.0, 0.0, 0.0, 
             0.015, 0.0, 0.0, 0.0, 1.0E-3)
names(defpars) <- parnames

if (!exists("defplotpars") ) defplotpars <- par(no.readonly = TRUE) # save default, for resetting

lmar   <- 6.0
rmar   <- 6.5
bmar   <- 5.5 
xlim   <- c(0.5, 1.5)
ylimTL <- c(0, 250)
ylimTR <- c(0.9, 2.2)
ylimBL <- c(0, 3.25)
ylimBR <- c(0.01, 0.15)

cexlab <- 1.6
cexaxs <- 1.6
cexttl <- 1.8
cexleg <- 1.4

ylablineL <- 4.0
ylablineR <- 5.0
xlabline  <- 3.0

axislwd <- 1
linelwd <- 3

parTI   <- 0.10
parMUI  <- 0.075

source("Figures/Compute-Equi-R0.R")

iszero <- 1.0E-8
params <- defpars
Requi <- uniroot.all(f = ReqFunc, interval = c(0.8, 1.0), tol = iszero, lower = 0.8, upper = 1.0, pars = params)
LRO <- ReqFunc(Requi, params) + 1
Bequi <- BeqFunc(Requi, params)
R0max <- R0Func(Requi, Bequi, params)

##############################################################################################################################
##### Compute data for the plot
##############################################################################################################################

qrange <- seq(0.5, 1.5, length.out = 1001)
Qres <- matrix(0, length(qrange), 12)
req <- 2.031
for (ind in (1:length(qrange))) {
  params <- defpars
  params["Q"]  <- qrange[ind]
  Qres[ind,  1] <- qrange[ind]
  Qres[ind,  2] <- uniroot.all(f = ReqFunc, tol = iszero, lower = 0.95*req, upper = 1.05*req, pars = params)
  Qres[ind,  3] <- BeqFunc(Qres[ind, 2], params)
  Qres[ind,  4] <- JnrFunc(Qres[ind, 2], Qres[ind, 3], params)
  Qres[ind,  5] <- AnrFunc(Qres[ind, 2], Qres[ind, 3], params)
  Qres[ind,  6] <- JbioFunc(Qres[ind, 2], Qres[ind, 3], params)
  Qres[ind,  7] <- AbioFunc(Qres[ind, 2], Qres[ind, 3], params)
  Qres[ind,  8] <- NtotFunc(Qres[ind, 2], Qres[ind, 3], params)
  Qres[ind,  9] <- ReqFunc(Qres[ind, 2], params) + 1
  Qres[ind, 10] <- R0Func(Qres[ind, 2], Qres[ind, 3], params)
  params["TI"] <- parTI
  Qres[ind, 11] <- R0Func(Qres[ind, 2], Qres[ind, 3], params)
  params["TI"] <- 0.0
  params["MUI"] <- parMUI
  Qres[ind, 12] <- R0Func(Qres[ind, 2], Qres[ind, 3], params)
  req <- Qres[ind, 2]
}
colnames(Qres) <- c("q", "R", "B", "Jnr", "Anr", "Jbio", "Abio", "Ntot", "LRO", "R0", "R0TI", "R0MUI")

##############################################################################################################################
##### Create the plot
##############################################################################################################################

par(defplotpars)
par(tcl = 0.5, mgp = c(3, 0.8, 0))

layout(matrix(c(1,2,3), 3, 1), heights = c(0.5, 0.3, 0.4))
# Top panel
par(mar=c(0, lmar, 0.1, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimTL)
axis(1, at=c(0.5, 0.75, 1.0, 1.25, 1.5), labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Juvenile-adult density", 2, line = ylablineL, cex = cexlab)

lines(Qres[,"q"], Qres[,"Jnr"], lwd = linelwd, col = "#0072B2")
lines(Qres[,"q"], Qres[,"Anr"], lwd = linelwd, col = "#D55E00")

legend("topleft", c("Juveniles", "Adults", "Resource"), col = c("#0072B2", "#D55E00", "#009E73"), lty = c(1, 1, 1), lwd = c(linelwd, linelwd, linelwd), cex = cexleg, horiz = T)

lblspos <- 1.0 + 0.2 * (0:7)
lbls <- c("0", "10", "20", "30", "40")
yusr <- par("usr")[3:4]
y2usr <- ylimTR + c(-0.04, 0.04) * (ylimTR[2] - ylimTR[1])
lines(Qres[,"q"], yusr[1] + (yusr[2] - yusr[1]) * (Qres[,"R"] - y2usr[1])/(y2usr[2] - y2usr[1]), lwd = linelwd, col = "#009E73")
axis(4, at = yusr[1] + (yusr[2] - yusr[1])*(lblspos - y2usr[1])/(y2usr[2] - y2usr[1]), labels = lblspos, cex.axis = cexaxs, 
     lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Resource density", 4, line = ylablineR, cex = cexlab)

# Middle panel
par(mar=c(0, lmar, 0, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimBL)

lines(Qres[,"q"], Qres[,"R0MUI"], lwd = linelwd, col = "black", lty = 1)
# lines(Qres[,"q"], Qres[,"R0"],    lwd = linelwd, col = "black", lty = 3)

axis(1, at=c(0.5, 0.75, 1.0, 1.25, 1.5), labels = F, 
     cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, at = (0:3), cex.axis = cexaxs, lwd = 0, lwd.ticks = 0, las = 2)
mtext(TeX("$R_0$"), 2, line = ylablineL - 0.4, cex = cexlab)

lblspos <- 0.02 + 0.02 * (0:7)
yusr <- par("usr")[3:4]
y2usr <- ylimBR + c(-0.04, 0.04) * (ylimBR[2] - ylimBR[1])

lines(xlim, yusr[1] + (yusr[2] - yusr[1]) * (rep(params["MUS"], 2) + parMUI - y2usr[1])/(y2usr[2] - y2usr[1]), 
      lwd = linelwd, col = "#0072B2")
lines(xlim, yusr[1] + (yusr[2] - yusr[1]) * (rep(params["MUS"], 2) + parMUI - y2usr[1])/(y2usr[2] - y2usr[1]), 
      lwd = linelwd, col = "#D55E00", lty = 2)

lines(xlim, yusr[1] + (yusr[2] - yusr[1]) * (rep(params["MUS"], 2) - y2usr[1])/(y2usr[2] - y2usr[1]), 
      lwd = 1, col = "black", lty = 2)

axis(4, at = yusr[1] + (yusr[2] - yusr[1])*(lblspos - y2usr[1])/(y2usr[2] - y2usr[1]), labels = F, 
     lwd = 0, lwd.ticks = axislwd)
axis(4, at = yusr[1] + (yusr[2] - yusr[1])*(0.02 + 0.04 * (0:3) - y2usr[1])/(y2usr[2] - y2usr[1]), labels = 0.02 + 0.04 * (0:3), cex.axis = cexaxs, 
     lwd = 0, lwd.ticks = 0, las = 2)
mtext("Mortality rate", 4, line = ylablineR, cex = cexlab)

legend("topleft", c(TeX("$R_0$"), "Infected juvenile mortality", "Infected adult mortality"), col = c("black", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = c(linelwd, linelwd, linelwd), cex = cexleg, horiz = T, text.width = NA)

# Bottom panel
par(mar=c(bmar, lmar, 0, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimBL)

lines(Qres[,"q"], Qres[,"R0TI"],  lwd = linelwd, col = "black", lty = 1)
# lines(Qres[,"q"], Qres[,"R0"],    lwd = linelwd, col = "black", lty = 3)

axis(1, at=c(0.5, 0.75, 1.0, 1.25, 1.5), labels = c("0.5", "0.75", "1", "1.25", "1.5"), 
     cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, at = (0:3), cex.axis = cexaxs, lwd = 0, lwd.ticks = 0, las = 2)
mtext("q", 1, line = xlabline, cex = cexlab)
mtext(TeX("$R_0$"), 2, line = ylablineL - 0.4, cex = cexlab)

lblspos <- 0.02 + 0.02 * (0:7)
yusr <- par("usr")[3:4]
y2usr <- ylimBR + c(-0.04, 0.04) * (ylimBR[2] - ylimBR[1])

funcresp <- params["SIGMA"] * params["M"] * Qres[,"R"] / (params["H"] + Qres[,"R"])
mortJ <- params["MUS"] - pmin((2.0 - Qres[,"q"]) * funcresp - (params["TS"] + parTI), 0)
mortA <- params["MUS"] - pmin((      Qres[,"q"]) * funcresp - (params["TS"] + parTI), 0)

lines(Qres[,"q"], yusr[1] + (yusr[2] - yusr[1]) * (mortJ - y2usr[1])/(y2usr[2] - y2usr[1]), 
      lwd = linelwd, col = "#0072B2")
lines(Qres[,"q"], yusr[1] + (yusr[2] - yusr[1]) * (mortA - y2usr[1])/(y2usr[2] - y2usr[1]), 
      lwd = linelwd, col = "#D55E00", lty = 2)

lines(xlim, yusr[1] + (yusr[2] - yusr[1]) * (rep(params["MUS"], 2) - y2usr[1])/(y2usr[2] - y2usr[1]), 
      lwd = 1, col = "black", lty = 2)

axis(4, at = yusr[1] + (yusr[2] - yusr[1])*(lblspos - y2usr[1])/(y2usr[2] - y2usr[1]), labels = F, 
     lwd = 0, lwd.ticks = axislwd)
axis(4, at = yusr[1] + (yusr[2] - yusr[1])*(0.02 + 0.04 * (0:3) - y2usr[1])/(y2usr[2] - y2usr[1]), labels = 0.02 + 0.04 * (0:3), cex.axis = cexaxs, 
     lwd = 0, lwd.ticks = 0, las = 2)

mtext("Mortality rate", 4, line = ylablineR, cex = cexlab)

legend("topleft", c(TeX("$R_0$"), "Infected juvenile mortality", "Infected adult mortality"), col = c("black", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = c(linelwd, linelwd, linelwd), cex = cexleg, horiz = T, text.width = NA)

par(defplotpars)

if (ToPdf) {
  dev.off()
  system(paste("open ", fname))
}
