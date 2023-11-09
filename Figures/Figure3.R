envvars <- ls()
rm(list = envvars[(envvars != "defplotpars") & 
                    (envvars != "Q075_MUI_Beta") & (envvars != "Q075_TI_Beta") & 
                    (envvars != "Q125_MUI_Beta") & (envvars != "Q125_TI_Beta")  & 
                    (envvars != "TI010_Beta_Qa") & (envvars != "TI010_Beta_Qb") & 
                    (envvars != "Q075_Beta_TIa") & (envvars != "Q075_Beta_TIb")])
rm("envvars")
basedir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../")
datadir <- paste0(basedir, "Figures/EBToutput/")
fname <- paste0(basedir, "Figures/Figure3.pdf")
setwd(basedir)

ToPdf <- T

if (ToPdf) pdf(file = fname, width = 6.0, height = 5.0)

library(rootSolve)
library(FindCurve)
library(latex2exp)
parnames <- c("Rho", "Rmax", "Sb", "Sm", "M", "Q", "H", "SIGMA", 
              "TS", "TI", "TIJ", "TIA",
              "MUS", "MUI", "MUIJ", "MUIA", "BETA")
defpars <- c(0.1, 100.0, 0.1, 1.0, 1.0, 1.0, 3.0, 0.5, 
             0.1, 0.0, 0.0, 0.0, 
             0.015, 0.0, 0.0, 0.0, 1.0E-3)
names(defpars) <- parnames

if (!exists("defplotpars") ) defplotpars <- par(no.readonly = TRUE) # save default, for resetting

lmar   <- 7.0
rmar   <- 0.1
bmar   <- 3.9 
xlim   <- c(0.5, 1.5)
ylim   <- c(0, 0.0025)
ylimTL <- c(0, 250)
ylimTR <- c(0.9, 2.2)
ylimBL <- c(0, 3.25)
ylimBR <- c(0.01, 0.15)

cexlab <- 1.8
cexaxs <- 1.3
cexttl <- 1.8
cexleg <- 1.5

ylablineL <- 6.0
ylablineR <- 5.0
xlabline  <- 2.5

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
Qres_MUI <- matrix(0, length(qrange), 12)
params <- defpars
params["MUI"] <- parMUI
req <- 2.031
for (ind in (1:length(qrange))) {
  params["Q"]  <- qrange[ind]
  Qres_MUI[ind,  1] <- qrange[ind]
  Qres_MUI[ind,  2] <- uniroot.all(f = ReqFunc, tol = iszero, lower = 0.95*req, upper = 1.05*req, pars = params)
  Qres_MUI[ind,  3] <- BeqFunc(Qres_MUI[ind, 2], params)
  Qres_MUI[ind,  4] <- JnrFunc(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  Qres_MUI[ind,  5] <- AnrFunc(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  Qres_MUI[ind,  6] <- JbioFunc(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  Qres_MUI[ind,  7] <- AbioFunc(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  Qres_MUI[ind,  8] <- NtotFunc(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  Qres_MUI[ind,  9] <- ReqFunc(Qres_MUI[ind, 2], params) + 1
  Qres_MUI[ind, 10] <- R0Func(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  Qres_MUI[ind, 11] <- params["BETA"] / Qres_MUI[ind, 10]
  oldbeta <- params["BETA"]
  params["BETA"] <- Qres_MUI[ind, 11]
  Qres_MUI[ind, 12] <- R0Func(Qres_MUI[ind, 2], Qres_MUI[ind, 3], params)
  params["BETA"] <- oldbeta
  req <- Qres_MUI[ind, 2]
}
colnames(Qres_MUI) <- c("q", "R", "B", "Jnr", "Anr", "Jbio", "Abio", "Ntot", "LRO", "R0", "BetaInvade", "TestR0")

qrange <- seq(0.5, 1.5, length.out = 1001)
Qres_TI <- matrix(0, length(qrange), 12)
params <- defpars
params["TI"] <- parTI
req <- 2.031
for (ind in (1:length(qrange))) {
  params["Q"]  <- qrange[ind]
  Qres_TI[ind,  1] <- qrange[ind]
  Qres_TI[ind,  2] <- uniroot.all(f = ReqFunc, tol = iszero, lower = 0.95*req, upper = 1.05*req, pars = params)
  Qres_TI[ind,  3] <- BeqFunc(Qres_TI[ind, 2], params)
  Qres_TI[ind,  4] <- JnrFunc(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  Qres_TI[ind,  5] <- AnrFunc(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  Qres_TI[ind,  6] <- JbioFunc(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  Qres_TI[ind,  7] <- AbioFunc(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  Qres_TI[ind,  8] <- NtotFunc(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  Qres_TI[ind,  9] <- ReqFunc(Qres_TI[ind, 2], params) + 1
  Qres_TI[ind, 10] <- R0Func(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  Qres_TI[ind, 11] <- params["BETA"] / Qres_TI[ind, 10]
  oldbeta <- params["BETA"]
  params["BETA"] <- Qres_TI[ind, 11]
  Qres_TI[ind, 12] <- R0Func(Qres_TI[ind, 2], Qres_TI[ind, 3], params)
  params["BETA"] <- oldbeta
  req <- Qres_TI[ind, 2]
}
colnames(Qres_TI) <- c("q", "R", "B", "Jnr", "Anr", "Jbio", "Abio", "Ntot", "LRO", "R0", "BetaInvade", "TestR0")

LPinit <- c(0.0008510851, 1.404672, 1.32286, 46.85942, 0.75)
if (!exists("TI010_Beta_Qa")) {
  pars <- defpars
  pars["TI"] <- parTI
  TI010_Beta_Qa <- FindCurve(modelname = "Equi/PSPM-Disease.h", 
                             curvetype = "LP", 
                             startpoint = LPinit, 
                             stepsize = -0.01, 
                             bounds   = c(0.0, 0.01, 0, 100, 0, 1000, 0, 1000, 0.676, 0.9), 
                             parameters = pars, 
                             options = c("par1", "16", "par2", "5", "noEXT"), 
                             clean = TRUE)
  colnames(TI010_Beta_Qa$curvepoints) <- c("par1", "R", "B", "Itot", "par2", "JnrS", "AnrS", "JnrI", "AnrI", "JbioS", "AbioS", 
                                             "JbioI", "AbioI", "PCbirthS", "PCbirthI", "TotMatS", "TotMatI", "RHS")
}

if (!exists("TI010_Beta_Qb")) {
  pars <- defpars
  pars["TI"] <- parTI
  TI010_Beta_Qb <- FindCurve(modelname = "Equi/PSPM-Disease.h", 
                             curvetype = "LP", 
                             startpoint = LPinit, 
                             stepsize = 0.01, 
                             bounds   = c(0.0, 0.006, 0, 100, 0, 1000, 1.0E-5, 1000, 0.676, 0.9), 
                             parameters = pars, 
                             options = c("par1", "16", "par2", "5", "noEXT"), 
                             clean = TRUE)
  colnames(TI010_Beta_Qb$curvepoints) <- c("par1", "R", "B", "Itot", "par2", "JnrS", "AnrS", "JnrI", "AnrI", "JbioS", "AbioS", 
                                             "JbioI", "AbioI", "PCbirthS", "PCbirthI", "TotMatS", "TotMatI", "RHS")
}

##############################################################################################################################
##### Create the plot
##############################################################################################################################

par(defplotpars)
par(tcl = 0.5, mgp = c(3, 0.6, 0))

# Plot the result

par(mar=c(bmar, lmar, 0.1, rmar))
plot(NULL, NULL, xlab = NA, ylab = NA, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")

LParea <- Qres_TI[(Qres_TI[, "q"] > 0.675) & (Qres_TI[, "q"] <= 0.9), c("q", "BetaInvade")]
LPbot  <- TI010_Beta_Qb$curvepoints[nrow(TI010_Beta_Qb$curvepoints):1, c("par2", "par1")]
LPbot  <- rbind(LPbot,  TI010_Beta_Qa$curvepoints[1:nrow(TI010_Beta_Qa$curvepoints), c("par2", "par1")])
LParea <- rbind(LParea, LPbot)

polygon(LParea[,1], LParea[,2], col = "grey", border = NA)
lines(LPbot[, 1], LPbot[, 2], lwd = 2, lty = 2, col = "red")

lines(Qres_TI[,"q"], Qres_TI[,"BetaInvade"], lwd = linelwd, col = "red")
lines(Qres_MUI[,"q"], Qres_MUI[,"BetaInvade"], lwd = linelwd, col = "black")

axis(1, at=c(0.5, 0.75, 1.0, 1.25, 1.5), 
     labels = c("0.5", "0.75", "1", "1.25", "1.5"), cex.axis = 1.4, lwd = 0, lwd.ticks = 1, cex.axis = cexaxs)
axis(2, cex.axis = 1.4, at = (0:5) * 5.0E-4, 
     labels = c(0, TeX(paste0('$', (1:5) * 0.5, '\\cdot 10^{-3}$'))), 
     lwd = 0, lwd.ticks = 1, las = 2, cex.axis = cexaxs)
mtext("q", 1.5, line = xlabline, cex = cexlab)
mtext(TeX("$\\beta$"), 2, line = ylablineL, cex = cexlab, las = 2)


par(defplotpars)

if (ToPdf) {
  dev.off()
  system(paste("open ", fname))
}
