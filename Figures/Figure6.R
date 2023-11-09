envvars <- ls()
rm(list = envvars[(envvars != "defplotpars") & 
                    (envvars != "Q075_MUI_Beta") & (envvars != "Q075_TI_Beta") & 
                    (envvars != "Q125_MUI_Beta") & (envvars != "Q125_TI_Beta")  & 
                    (envvars != "TI010_Beta_Qa") & (envvars != "TI010_Beta_Qb") & 
                    (envvars != "Q075_Beta_TIa") & (envvars != "Q075_Beta_TIb")])
rm("envvars")
basedir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../")
datadir <- paste0(basedir, "Figures/EBToutput/")
fname <- paste0(basedir, "Figures/Figure6.pdf")
setwd(basedir)

ToPdf <- T
CompareEBT <- F

if (ToPdf) pdf(file = fname, width = 8.0, height = 6.0)

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
linecol <- "black"
rescol  <- "#009E73"

ebtlwd  <- 4
saddlecol <- "darkgrey"
saddlelwd <- 2
saddlelty <- 2

hopfcol <- "darkgrey"
hopflwd <- 2
hopflty <- 2

parTI   <- 0.10
parMUI  <- 0.075

xlim   <- c(0.0, 0.0035)
ylimTL <- c(38, 81)
ylimTR <- c(0.9, 2.2)
ylimBL <- c(0.8, 2.25)

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
                             options = c("par1", "16", "noEXT", "report", "100"), 
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
                            options = c("par1", "16", "par2", "9", "noEXT", "report", "10"), 
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
MUIavgU <- read.table(paste0(datadir, "/Q075_MUI0075_Beta_up.avg.out"))
colnames(MUIavgU) <- c("Time", "R",
                       "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                       "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                       "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                       "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim",
                       "par1", "Period", "CohNrS2", "CohNrI2")

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
TIavgU <- EBT <- read.table(paste0(datadir, "/Q075_TI010_Beta_up.avg.out"))
colnames(TIavgU) <- c("Time", "R",
                      "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                      "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                      "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                      "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim",
                      "par1", "Period", "CohNrS2", "CohNrI2")

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
MUIavgS    <- mean(MUIavgU[(MUImaxU[, "Period"] > 0) & (MUImaxU[, "par1"] < 0.00154), "Tnr"])

# Extract the susceptible-infected cycle
cyclesindx <- (1:nrow(MUImaxD))[(MUImaxD[, "Period"] > 0) & (MUImaxD[, "TnrI"] > 0.01)]
cyclesD    <- split(cyclesindx,cumsum(c(1,diff(cyclesindx)!=1)))
cyclesD    <- cyclesD[lapply(cyclesD,length)>1]             # Remove elements of length 1

MUIminC <- list(MUIminU[cyclesU[[1]],], MUIminD[cyclesD[[1]],])
MUImaxC <- list(MUImaxU[cyclesU[[1]],], MUImaxD[cyclesD[[1]],])

par(defplotpars)
par(tcl = 0.5, mgp = c(3, 0.8, 0))
layout(matrix((1:4), 2, 2), heights = c(0.5, 0.4), widths = c(1.2, 1.0))

# Top panel
par(mar=c(0, lmar, tmar, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimTL)
axis(1, at=(0:6) * 1.0E-3, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Total density", 2, line = ylablineL, cex = cexlab)

cyclesindx <- (1:nrow(Q075_MUI_Beta$curvepoints))[Q075_MUI_Beta$curvepoints[,"par1"] < max(MUIminC[[length(MUIminC)]][,"par1"])]
h1indx <- min(cyclesindx) - 1

totdens <- rowSums(Q075_MUI_Beta$curvepoints[,c("JnrS", "AnrS", "JnrI", "AnrI")])
lines(smooth.spline(Q075_MUI_Beta$curvepoints[(1:h1indx),"par1"], totdens[(1:h1indx)]), lwd = linelwd, col = linecol)

lines(smooth.spline(Q075_MUI_Beta$curvepoints[cyclesindx,"par1"], totdens[cyclesindx]), lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(c(xlim[1], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(totdens, 1), 2),
      lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(range(MUIavgU[(MUImaxU[, "Period"] > 0)  & (MUImaxU[, "par1"] < 0.00154), "par1"]), rep(MUIavgS, 2), lwd = 2, lty = 1, col = linecol)

lines(c(par("usr")[2], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(totdens, 1), 2),
      lwd = saddlelwd, lty = saddlelty, col = saddlecol)

if (CompareEBT) {
  points(MUIminD[, "par1"], MUIminD[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxD[, "par1"], MUImaxD[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUIminU[, "par1"], MUIminU[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(MUImaxU[, "par1"], MUImaxU[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(MUIminC))) {
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "Tnr")]), lwd = ebtlwd, col = linecol)
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "Tnr")]), lwd = ebtlwd, col = linecol)
  }
}

# Bottom panel
par(mar=c(bmar, lmar, 0, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimBL)

axis(1, at=(0:6) * 1.0E-3, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext(TeX("$\\beta$"), 1, line = xlabline, cex = cexlab)
mtext("Resource density", 2, line = ylablineL, cex = cexlab)

lines(smooth.spline(Q075_MUI_Beta$curvepoints[(1:h1indx), c("par1", "R")]), lwd = linelwd, col = rescol)
lines(smooth.spline(Q075_MUI_Beta$curvepoints[cyclesindx, c("par1", "R")]), lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(c(xlim[1], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(Q075_MUI_Beta$curvepoints, 1)[,"R"], 2),
      lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(c(par("usr")[2], tail(Q075_MUI_Beta$curvepoints, 1)[,"par1"]), rep(tail(Q075_MUI_Beta$curvepoints, 1)[,"R"], 2),
      lwd = saddlelwd, lty = saddlelty, col = saddlecol)

params <- defpars
params["MUI"] <- 0.075
params["Q"]  <- 0.75

# Rcrit <- params["TS"] * params["H"] / (params["SIGMA"] * (2 - params["Q"]) * params["M"] - params["TS"])
# lines(par("usr")[1:2], rep(Rcrit, 2), lty = 2, lwd = 1, col = "black")
# 
Rcrit <- params["TS"] * params["H"] / (params["SIGMA"] *      params["Q"]  * params["M"] - params["TS"])
lines(par("usr")[1:2], rep(Rcrit, 2), lty = 5, lwd = 2, col = "#D55E00")

if (CompareEBT) {
  points(TIminD[, "par1"], MUIminD[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], MUImaxD[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], MUIminU[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], MUImaxU[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(MUIminC))) {
    lines(smooth.spline(MUIminC[[ii]][, c("par1", "R")]), lwd = ebtlwd, col = rescol)
    lines(smooth.spline(MUImaxC[[ii]][, c("par1", "R")]), lwd = ebtlwd, col = rescol)
  }
}

##############################################################################################################################
##### Right column
##############################################################################################################################
# Extract the susceptible only cycle. Numerical simulations show this one is invadable from TI = 0.00118 and above
cyclesindx <- (1:nrow(TImaxU))[(TImaxU[, "Period"] > 0)  & (TImaxU[, "par1"] < 0.00118)]
cyclesU    <- split(cyclesindx,cumsum(c(1,diff(cyclesindx)!=1)))
cyclesU    <- cyclesU[lapply(cyclesU,length)>1]             # Remove elements of length 1
TIavgS     <- mean(TIavgU[(TImaxU[, "Period"] > 0)  & (TImaxU[, "par1"] < 0.00118), "Tnr"])

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

totdens <- rowSums(Q075_TI_Beta$curvepoints[,c("JnrS", "AnrS", "JnrI", "AnrI")])
lines(smooth.spline(Q075_TI_Beta$curvepoints[(1:h1indx),"par1"], totdens[(1:h1indx)]), lwd = linelwd, col = linecol)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx), "par1"], totdens[(h1indx:h2indx)]), lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx), "par1"], totdens[(lpindx:lastindx)]), lwd = saddlelwd, lty = saddlelty, col = saddlecol)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h2indx:lpindx),"par1"], totdens[(h2indx:lpindx)]), lwd = linelwd, col = linecol)

lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(totdens[lastindx], 2),
      lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(range(TIavgU[(TImaxU[, "Period"] > 0)  & (TImaxU[, "par1"] < 0.00118), "par1"]), rep(TIavgS, 2), lwd = 2, lty = 1, col = linecol)

lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(totdens[lastindx], 2),
      lwd = saddlelwd, lty = saddlelty, col = saddlecol)

if (CompareEBT) {
  points(TIminD[, "par1"], TIminD[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "Tnr"], col = "#56B4E9", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(TIminC))) {
    lines(smooth.spline(TIminC[[ii]][, c("par1", "Tnr")]), lwd = ebtlwd, col = linecol)
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "Tnr")]), lwd = ebtlwd, col = linecol)
  }
}

# Bottom panel
par(mar=c(bmar, 0, 0, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimBL)

axis(1, at=(0:6) * 1.0E-3, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext(TeX("$\\beta$"), 1, line = xlabline, cex = cexlab)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(1:h1indx), c("par1", "R")]), lwd = linelwd, col = rescol)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h1indx:h2indx), c("par1", "R")]), lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(h2indx:lpindx), c("par1", "R")]), lwd = linelwd, col = rescol)

lines(smooth.spline(Q075_TI_Beta$curvepoints[(lpindx:lastindx), c("par1", "R")]), lwd = saddlelwd, lty = saddlelty, col = saddlecol)

lines(c(xlim[1], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(tail(Q075_TI_Beta$curvepoints, 1)[, "R"], 2), lwd = hopflwd, lty = hopflty, col = hopfcol)

lines(c(par("usr")[2], Q075_TI_Beta$curvepoints[lastindx, "par1"]), rep(tail(Q075_TI_Beta$curvepoints, 1)[, "R"], 2), lwd = saddlelwd, lty = saddlelty, col = saddlecol)

params <- defpars
params["TI"] <- 0.10
params["Q"]  <- 0.75

# Rcrit <- params["TS"] * params["H"] / (params["SIGMA"] * (2 - params["Q"]) * params["M"] - params["TS"])
# lines(par("usr")[1:2], rep(Rcrit, 2), lty = 2, lwd = 1, col = "black")
# 
Rcrit <- params["TS"] * params["H"] / (params["SIGMA"] *      params["Q"]  * params["M"] - params["TS"])
lines(par("usr")[1:2], rep(Rcrit, 2), lty = 5, lwd = 2, col = "#D55E00")
 
Rcrit <- (params["TS"] + params["TI"] + params["TIJ"]) * params["H"] / (params["SIGMA"] * (2 - params["Q"]) * params["M"] - (params["TS"] + params["TI"] + params["TIJ"]))
lines(par("usr")[1:2], rep(Rcrit, 2), lty = 5, lwd = 2, col = "#0072B2")

# Rcrit <- (params["TS"] + params["TI"] + params["TIA"]) * params["H"] / (params["SIGMA"] *      params["Q"] *  params["M"] - (params["TS"] + params["TI"] + params["TIA"]))
# lines(par("usr")[1:2], rep(Rcrit, 2), lty = 3, lwd = 1, col = "#D55E00")

if (CompareEBT) {
  points(TIminD[, "par1"], TIminD[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxD[, "par1"], TImaxD[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TIminU[, "par1"], TIminU[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
  points(TImaxU[, "par1"], TImaxU[, "R"], col = "#56B4E9", pch = 19, cex = 0.2)
} else {
  for (ii in (1:length(TIminC))) {
    lines(smooth.spline(TIminC[[ii]][, c("par1", "R")]), lwd = ebtlwd, col = rescol)
    lines(smooth.spline(TImaxC[[ii]][, c("par1", "R")]), lwd = ebtlwd, col = rescol)
  }
}

par(defplotpars)

if (ToPdf) {
  dev.off()
  system(paste("open ", fname))
}
