envvars <- ls()
rm(list = envvars[(envvars != "defplotpars") & 
                    (envvars != "Q075_MUI_Beta") & (envvars != "Q075_TI_Beta") & 
                    (envvars != "Q125_MUI_Beta") & (envvars != "Q125_TI_Beta")  & 
                    (envvars != "TI010_Beta_Qa") & (envvars != "TI010_Beta_Qb") & 
                    (envvars != "Q075_Beta_TIa") & (envvars != "Q075_Beta_TIb")])
rm("envvars")
basedir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../")
datadir <- paste0(basedir, "Figures/EBToutput/")
fname <- paste0(basedir, "Figures/Figure5.pdf")
setwd(basedir)

ToPdf <- T
CompareEBT <- F

if (ToPdf) pdf(file = fname, width = 7.0, height = 8.0)

library(FindCurve)
library(latex2exp)
require(graphics)

parnames <- c("Rho", "Rmax", "Sb", "Sm", "M", "Q", "H", "SIGMA", 
              "TS", "TI", "TIJ", "TIA",
              "MUS", "MUI", "MUIJ", "MUIA", "BETA")
defpars <- c(0.1, 100.0, 0.1, 1.0, 1.0, 1.0, 3.0, 0.5, 
             0.1, 0.0, 0.0, 0.0, 
             0.015, 0.0, 0.0, 0.0, 1.0E-3)
names(defpars) <- parnames

if (!exists("defplotpars") ) defplotpars <- par(no.readonly = TRUE) # save default, for resetting

lmar   <- 6.5
rmar   <- 0.2
tmar   <- 0.2
bmar   <- 5.5 

cexlab <- 1.6
cexaxs <- 1.6
cexttl <- 1.8
cexleg <- 1.5

ylablineL <- 4.4
ylablineR <- 5.0
xlabline  <- 4.0

axislwd <- 1
linelwd <- 4
ebtlwd  <- 4

saddlelwd <- 3
saddlelty <- 2

hopflwd <- 3
hopflty <- 2

parTI   <- 0.10
parMUI  <- 0.075

xlim   <- c(0.0, 0.0035)
ylimTL <- c(0, 50)
ylimTR <- c(0.9, 2.2)
ylimML <- c(0, 55)
ylimBL <- c(0.01, 0.15)

##############################################################################################################################
##### Compute data for the plot
##############################################################################################################################

if (!exists("Q075_MUI_Beta")) {
  params <- defpars
  params["Q"] <- 0.75
  params["MUI"] <- parMUI
  init <- c(0.006, 1.67979513E+00, 4.31034951E+00, 4.53927322E+01)
  Q075_MUI_Beta <- FindCurve(modelname = "Equi/PSPM-Disease.h", 
                             curvetype = "EQ", 
                             startpoint = init, 
                             stepsize = -0.001, 
                             bounds   = c(0.0, 0.006, 0, 100, 0, 1000, 1.0E-5, 1000), 
                             parameters = params, 
                             options = c("par1", "16", "noEXT", "report", "10"), 
                             clean = TRUE)
  colnames(Q075_MUI_Beta$curvepoints) <- c("par1", "R", "B", "Itot", "par2", "JnrS", "AnrS", "JnrI", "AnrI", "JbioS", "AbioS", 
                                           "JbioI", "AbioI", "PCbirthS", "PCbirthI", "TotMatS", "TotMatI", "RHS")
}

if (!exists("Q075_TI_Beta")) {
  params <- defpars
  params["Q"] <- 0.75
  params["TI"] <- parTI
  init <- c(0.006, 3.077372, 0.6179408, 24.41046)
  Q075_TI_Beta <- FindCurve(modelname = "Equi/PSPM-Disease.h", 
                            curvetype = "EQ", 
                            startpoint = init, 
                            stepsize = -0.001, 
                            bounds   = c(0.0, 0.006, 0, 100, 0, 1000, 1.0E-5, 1000), 
                            parameters = params, 
                            options = c("par1", "16", "par2", "9", "noEXT", "report", "100"), 
                            clean = TRUE)
  colnames(Q075_TI_Beta$curvepoints) <- c("par1", "R", "B", "Itot", "par2", "JnrS", "AnrS", "JnrI", "AnrI", "JbioS", "AbioS", 
                                          "JbioI", "AbioI", "PCbirthS", "PCbirthI", "TotMatS", "TotMatI", "RHS")
  colnames(Q075_TI_Beta$bifpoints) <- c("par1", "R", "B", "Itot", "par2", "JnrS", "AnrS", "JnrI", "AnrI", "JbioS", "AbioS", 
                                        "JbioI", "AbioI", "PCbirthS", "PCbirthI", "TotMatS", "TotMatI", "RHS")
}

# Read the EBT data

EBT <- read.table(paste0(datadir, "/Q075_MUI0075_Beta_up.minmax.out"))
colnames(EBT) <- c("Time", "R",
                   "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                   "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                   "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                   "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim",
                   "par1", "Period", "CohNrS2", "CohNrI2")
oddrows  <- 2*(1:(nrow(EBT)/2)) - 1
evenrows <- 2*(1:(nrow(EBT)/2))
MUIminU <- EBT[oddrows, ]
MUImaxU <- EBT[evenrows,]

EBT <- read.table(paste0(datadir, "/Q075_MUI0075_Beta_dwn.minmax.out"))
colnames(EBT) <- c("Time", "R",
                   "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                   "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                   "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                   "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim",
                   "par1", "Period", "CohNrS2", "CohNrI2")
oddrows  <- 2*(1:(nrow(EBT)/2)) - 1
evenrows <- 2*(1:(nrow(EBT)/2))
MUIminD <- EBT[rev(oddrows), ]                                     # Reorder from low to high parameter values
MUImaxD <- EBT[rev(evenrows),]

EBT <- read.table(paste0(datadir, "/Q075_TI010_Beta_up.minmax.out"))
colnames(EBT) <- c("Time", "R",
                   "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                   "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                   "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                   "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim",
                   "par1", "Period", "CohNrS2", "CohNrI2")
oddrows  <- 2*(1:(nrow(EBT)/2)) - 1
evenrows <- 2*(1:(nrow(EBT)/2))
TIminU <- EBT[oddrows, ]
TImaxU <- EBT[evenrows,]

EBT <- read.table(paste0(datadir, "/Q075_TI010_Beta_dwn.minmax.out"))
colnames(EBT) <- c("Time", "R",
                   "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                   "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                   "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                   "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim",
                   "par1", "Period", "CohNrS2", "CohNrI2")
oddrows  <- 2*(1:(nrow(EBT)/2)) - 1
evenrows <- 2*(1:(nrow(EBT)/2))
TIminD <- EBT[rev(oddrows), ]                                     # Reorder from low to high parameter values
TImaxD <- EBT[rev(evenrows),]

##############################################################################################################################
##### Create the plot
##############################################################################################################################
##### Left column
##############################################################################################################################

# Extract the susceptible only cycle. Numerical simulations show this one is invadable from mui = 0.00154 and above
cyclesindx <- (1:nrow(MUImaxU))[(MUImaxU[, "Period"] > 0) & (MUImaxU[, "par1"] < 0.00154)]
cyclesU    <- split(cyclesindx,cumsum(c(1,diff(cyclesindx)!=1)))
cyclesU    <- cyclesU[lapply(cyclesU,length)>1]             # Remove elements of length 1

# Extract the susceptible-infected cycle
cyclesindx <- (1:nrow(MUImaxD))[(MUImaxD[, "Period"] > 0) & (MUImaxD[, "TnrI"] > 0.01)]
cyclesD    <- split(cyclesindx,cumsum(c(1,diff(cyclesindx)!=1)))
cyclesD    <- cyclesD[lapply(cyclesD,length)>1]             # Remove elements of length 1

MUIminC <- list(MUIminU[cyclesU[[1]],], MUIminD[cyclesD[[1]],])
MUImaxC <- list(MUImaxU[cyclesU[[1]],], MUImaxD[cyclesD[[1]],])

par(defplotpars)
par(tcl = 0.5, mgp = c(3, 0.8, 0))
layout(matrix((1:6), 3, 2), heights = c(0.5, 0.3, 0.4), widths = c(1.2, 1.0))

# Top panel
par(mar=c(0, lmar, tmar, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimTL)
axis(1, at=(0:6) * 1.0E-3, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Susceptible density", 2, line = ylablineL, cex = cexlab)

cyclesindx <- (1:nrow(Q075_MUI_Beta$curvepoints))[Q075_MUI_Beta$curvepoints[,"par1"] < max(MUIminC[[length(MUIminC)]][,"par1"])]
h1indx <- min(cyclesindx) - 1

lines(Q075_MUI_Beta$curvepoints[(1:h1indx),"par1"], Q075_MUI_Beta$curvepoints[(1:h1indx),"JnrS"], lwd = linelwd, col = "#0072B2")
lines(Q075_MUI_Beta$curvepoints[(1:h1indx),"par1"], Q075_MUI_Beta$curvepoints[(1:h1indx),"AnrS"], lwd = linelwd, col = "#D55E00")

lines(Q075_MUI_Beta$curvepoints[cyclesindx,"par1"], Q075_MUI_Beta$curvepoints[cyclesindx,"JnrS"], lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(Q075_MUI_Beta$curvepoints[cyclesindx,"par1"], Q075_MUI_Beta$curvepoints[cyclesindx,"AnrS"], lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(xlim[1], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(Q075_MUI_Beta$curvepoints, 1)[,"JnrS"], 2),
      lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(c(xlim[1], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(Q075_MUI_Beta$curvepoints, 1)[,"AnrS"], 2),
      lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(par("usr")[2], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(Q075_MUI_Beta$curvepoints, 1)[,"JnrS"], 2),
      lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(c(par("usr")[2], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(Q075_MUI_Beta$curvepoints, 1)[,"AnrS"], 2),
      lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

if (CompareEBT) {
  points(MUIminD[, "par1"], MUIminD[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  
  points(MUIminD[, "par1"], MUIminD[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(MUIminC))) {
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "JnrS")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "AnrS")]), lwd = ebtlwd, col = "#D55E00")
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "JnrS")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "AnrS")]), lwd = ebtlwd, col = "#D55E00")
  }
}

# Middle panel
par(mar=c(0, lmar, 0, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimML)
axis(1, at=(0:6) * 1.0E-3, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Infected density", 2, line = ylablineL, cex = cexlab)

lines(Q075_MUI_Beta$curvepoints[(1:h1indx),"par1"], Q075_MUI_Beta$curvepoints[(1:h1indx),"JnrI"], lwd = linelwd, col = "#0072B2")
lines(Q075_MUI_Beta$curvepoints[(1:h1indx),"par1"], Q075_MUI_Beta$curvepoints[(1:h1indx),"AnrI"], lwd = linelwd, col = "#D55E00")

lines(Q075_MUI_Beta$curvepoints[cyclesindx,"par1"], Q075_MUI_Beta$curvepoints[cyclesindx,"JnrI"], lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(Q075_MUI_Beta$curvepoints[cyclesindx,"par1"], Q075_MUI_Beta$curvepoints[cyclesindx,"AnrI"], lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(xlim[1], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(0, 2), lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(c(xlim[1], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(0, 2), lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(par("usr")[2], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(0, 2), lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(c(par("usr")[2], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(0, 2), lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

if (CompareEBT) {
  points(MUIminD[, "par1"], MUIminD[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  
  points(MUIminD[, "par1"], MUIminD[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(MUIminC))) {
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "JnrI")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "AnrI")]), lwd = ebtlwd, col = "#D55E00")
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "JnrI")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "AnrI")]), lwd = ebtlwd, col = "#D55E00")
  }
}

legend("topleft", c("Juveniles", "Adults"), col = c("#0072B2", "#D55E00"), lty = c(1, 1), lwd = c(linelwd, linelwd), cex = cexleg, horiz = T)

# Bottom panel
par(mar=c(bmar, lmar, 0, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimBL)

axis(1, at=(0:6) * 1.0E-3, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, at = 0.02 + 0.04 * (0:3), cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext(TeX("$\\beta$"), 1, line = xlabline, cex = cexlab)
mtext("Mortality rate", 2, line = ylablineL, cex = cexlab)

lines(xlim, rep(defpars["MUS"], 2) + parMUI, lwd = linelwd, col = "#0072B2")
lines(xlim, rep(defpars["MUS"], 2) + parMUI, lwd = linelwd, col = "#D55E00", lty = 2)

lines(xlim, rep(defpars["MUS"], 2), lwd = 1, col = "black", lty = 2)

if (CompareEBT) {
  points(MUIminD[, "par1"], MUIminD[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  
  points(MUIminD[, "par1"], MUIminD[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(MUIminC))) {
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "mortIJ")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "mortIA")]), lwd = ebtlwd, col = "#D55E00")
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "mortIJ")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "mortIA")]), lwd = ebtlwd, col = "#D55E00")
  }
}

##############################################################################################################################
##### Right column
##############################################################################################################################
# Extract the susceptible only cycle. Numerical simulations show this one is invadable from TI = 0.00118 and above
cyclesindx <- (1:nrow(TImaxU))[(TImaxU[, "Period"] > 0)  & (TImaxU[, "par1"] < 0.00118)]
cyclesU    <- split(cyclesindx,cumsum(c(1,diff(cyclesindx)!=1)))
cyclesU    <- cyclesU[lapply(cyclesU,length)>1]             # Remove elements of length 1

# Extract the susceptible-infected cycle
cyclesindx <- (1:nrow(TImaxD))[(TImaxD[, "Period"] > 0) & (TImaxD[, "TnrI"] > 0.01)]
cyclesD    <- split(cyclesindx,cumsum(c(1,diff(cyclesindx)!=1)))
cyclesD    <- cyclesD[lapply(cyclesD,length)>1]             # Remove elements of length 1

TIminC <- list(TIminU[head(cyclesU[[1]], -1),], TIminD[cyclesD[[1]],])
TImaxC <- list(TImaxU[head(cyclesU[[1]], -1),], TImaxD[cyclesD[[1]],])

# Top panel
par(mar=c(0, 0, tmar, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimTL)
axis(1, at=(0:6) * 1.0E-3, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)

lpindx <- which.min(Q075_TI_Beta$curvepoints[,"par1"])
lastindx <- nrow(Q075_TI_Beta$curvepoints)

cyclesindx <- (1:nrow(Q075_TI_Beta$curvepoints))[(Q075_TI_Beta$curvepoints[,"par1"] > min(TIminC[[2]][,"par1"])) & (Q075_TI_Beta$curvepoints[,"par1"] < max(TIminC[[2]][,"par1"]))]
h1indx <- min(cyclesindx) - 1
h2indx <- max(cyclesindx) + 1

lines(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], Q075_TI_Beta$curvepoints[(1:h1indx),"JnrS"], lwd = linelwd, col = "#0072B2")
lines(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], Q075_TI_Beta$curvepoints[(1:h1indx),"AnrS"], lwd = linelwd, col = "#D55E00")

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx),c("par1", "JnrS")]), lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx),c("par1", "AnrS")]), lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], Q075_TI_Beta$curvepoints[(h2indx:lpindx),"JnrS"], lwd = linelwd, col = "#0072B2")
lines(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], Q075_TI_Beta$curvepoints[(h2indx:lpindx),"AnrS"], lwd = linelwd, col = "#D55E00")

lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx),c("par1", "JnrS")]), lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx),c("par1", "AnrS")]), lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(Q075_TI_Beta$curvepoints[lastindx, "JnrS"], 2),
      lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(Q075_TI_Beta$curvepoints[lastindx, "AnrS"], 2),
      lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(Q075_TI_Beta$curvepoints[lastindx, "JnrS"], 2),
      lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(Q075_TI_Beta$curvepoints[lastindx, "AnrS"], 2),
      lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

if (CompareEBT) {
  points(TIminD[, "par1"], TIminD[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "JnrS"], col = "#56B4E9", pch = 19, cex = 0.2)
  
  points(TIminD[, "par1"], TIminD[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "AnrS"], col = "#F0E442", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(TIminC))) {
    lines(smooth.spline(TIminC[[ii]][, c("par1", "JnrS")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(TIminC[[ii]][, c("par1", "AnrS")]), lwd = ebtlwd, col = "#D55E00")
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "JnrS")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "AnrS")]), lwd = ebtlwd, col = "#D55E00")
  }
}

# Middle panel
par(mar=c(0, 0, 0, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimML)
axis(1, at=(0:6) * 1.0E-3, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)

lines(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], Q075_TI_Beta$curvepoints[(1:h1indx),"JnrI"], lwd = linelwd, col = "#0072B2")
lines(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], Q075_TI_Beta$curvepoints[(1:h1indx),"AnrI"], lwd = linelwd, col = "#D55E00")

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx),c("par1", "JnrI")]), lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx),c("par1", "AnrI")]), lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], Q075_TI_Beta$curvepoints[(h2indx:lpindx),"JnrI"], lwd = linelwd, col = "#0072B2")
lines(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], Q075_TI_Beta$curvepoints[(h2indx:lpindx),"AnrI"], lwd = linelwd, col = "#D55E00")

lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx),c("par1", "JnrI")]), lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx),c("par1", "AnrI")]), lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(0, 2), lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(0, 2), lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(0, 2), lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(0, 2), lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

if (CompareEBT) {
  points(TIminD[, "par1"], TIminD[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "JnrI"], col = "#56B4E9", pch = 19, cex = 0.2)
  
  points(TIminD[, "par1"], TIminD[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "AnrI"], col = "#F0E442", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(TIminC))) {
    lines(smooth.spline(TIminC[[ii]][, c("par1", "JnrI")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(TIminC[[ii]][, c("par1", "AnrI")]), lwd = ebtlwd, col = "#D55E00")
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "JnrI")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "AnrI")]), lwd = ebtlwd, col = "#D55E00")
  }
}

# Bottom panel
par(mar=c(bmar, 0, 0, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimBL)

axis(1, at=(0:6) * 1.0E-3, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, at = 0.02 + 0.04 * (0:3), labels = F, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext(TeX("$\\beta$"), 1, line = xlabline, cex = cexlab)

funcresp <- defpars["SIGMA"] * defpars["M"] * Q075_TI_Beta$curvepoints[,"R"] / (defpars["H"] + Q075_TI_Beta$curvepoints[,"R"])
mortJ <- defpars["MUS"] - pmin((2.0 - 0.75) * funcresp - (defpars["TS"] + parTI), 0)
mortA <- defpars["MUS"] - pmin((      0.75) * funcresp - (defpars["TS"] + parTI), 0)

lines(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], mortJ[(1:h1indx)], lwd = linelwd, col = "#0072B2")
lines(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], mortA[(1:h1indx)], lwd = linelwd, col = "#D55E00")

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx), "par1"], mortJ[(h1indx:h2indx)]), lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx), "par1"], mortA[(h1indx:h2indx)]), lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], mortJ[(h2indx:lpindx)], lwd = linelwd, col = "#0072B2")
lines(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], mortA[(h2indx:lpindx)], lwd = linelwd, col = "#D55E00")

lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx), "par1"], mortJ[(lpindx:lastindx)]), lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx), "par1"], mortA[(lpindx:lastindx)]), lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(mortJ[length(mortJ)], 2), lwd = hopflwd, lty = hopflty, col = "#0072B2")
lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(mortA[length(mortA)], 2), lwd = hopflwd, lty = hopflty, col = "#D55E00")

lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(mortJ[length(mortJ)], 2), lwd = saddlelwd, lty = saddlelty, col = "#0072B2")
lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(mortA[length(mortA)], 2), lwd = saddlelwd, lty = saddlelty, col = "#D55E00")

lines(xlim, rep(defpars["MUS"], 2), lwd = 1, col = "black", lty = 2)

if (CompareEBT) {
  points(TIminD[, "par1"], TIminD[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "mortIJ"], col = "#56B4E9", pch = 19, cex = 0.2)
  
  points(TIminD[, "par1"], TIminD[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "mortIA"], col = "#F0E442", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(TIminC))) {
    lines(smooth.spline(TIminC[[ii]][, c("par1", "mortIJ")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(TIminC[[ii]][, c("par1", "mortIA")]), lwd = ebtlwd, col = "#D55E00")
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "mortIJ")]), lwd = ebtlwd, col = "#0072B2")
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "mortIA")]), lwd = ebtlwd, col = "#D55E00")
  }
}

par(defplotpars)

if (ToPdf) {
  dev.off()
  system(paste("open ", fname))
}
