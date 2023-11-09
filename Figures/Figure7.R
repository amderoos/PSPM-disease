envvars <- ls()
rm(list = envvars[(envvars != "defplotpars") & 
                    (envvars != "Q075_MUI_Beta") & (envvars != "Q075_TI_Beta") & 
                    (envvars != "Q125_MUI_Beta") & (envvars != "Q125_TI_Beta")  & 
                    (envvars != "TI010_Beta_Qa") & (envvars != "TI010_Beta_Qb") & 
                    (envvars != "Q075_Beta_TIa") & (envvars != "Q075_Beta_TIb")])
rm("envvars")
basedir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../")
datadir <- paste0(basedir, "Figures/EBToutput/")
fname <- paste0(basedir, "Figures/Figure7.pdf")
setwd(basedir)

ToPdf <- T

if (ToPdf) pdf(file = fname, width = 8.0, height = 6.0)

if (!exists("defplotpars") ) defplotpars <- par(no.readonly = TRUE) # save default, for resetting

lmar   <- 5.0
rmar   <- 0.7
tmar   <- 0.2
bmar   <- 4.0

cexlab <- 1.6
cexaxs <- 1.6
cexttl <- 1.8
cexleg <- 1.3

ylablineL <- 3.5
ylablineR <- 4.0
xlabline  <- 3.0

axislwd <- 1
linelwd <- 2

xlim   <- c(0.0, 1300)
ylimTL <- c(0, 45)
ylimML <- c(0, 45)

##############################################################################################################################
##### Read data for the plot
##############################################################################################################################

EBText <- read.table(paste0(datadir, "/Q075_TI010_Beta0001_IJ14.out"))
EBTinv <- read.table(paste0(datadir, "/Q075_TI010_Beta0001_IJ15.out"))

colnames(EBText) <- c("Time", "R",
                      "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                      "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                      "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                      "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim")
colnames(EBTinv) <- c("Time", "R",
                      "JbioS", "JbioI", "JBioT", "AbioS", "AbioI", "ABioT", "TbioS", "TbioI", "TBio",
                      "JnrS",  "JnrI",  "JnrT",  "AnrS",  "AnrI",  "AnrT",  "TnrS",  "TnrI",  "Tnr",
                      "Ingest", "nuJ", "nuIJ", "nuA", "nuIA", "mortJ", "mortIJ", "mortA", "mortIA", "FecS", "FecI",
                      "BRbioS", "BRbioI", "CohNrS", "CohNrI", "CohLim")

##############################################################################################################################
##### Create the plot
##############################################################################################################################
##### Left column
##############################################################################################################################

par(defplotpars)
par(tcl = 0.5, mgp = c(3, 0.8, 0))

layout(matrix((1:4), 2, 2), heights = c(0.4, 0.5), widths = c(1.2, 1.0))
# Top panel
par(mar=c(0, lmar, tmar, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimTL)
axis(1, at=(0:6) * 200, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Susceptible density", 2, line = ylablineL, cex = cexlab)

lines(EBText[,"Time"], EBText[,"JnrS"], lwd = linelwd, col = "#0072B2")
lines(EBText[,"Time"], EBText[,"AnrS"], lwd = linelwd, col = "#D55E00")
# lines(EBText[,"Time"], EBText[,"TnrS"], lwd = linelwd, col = "black")

lblspos <- 1.0 + 0.2 * (0:7)
lbls <- c("0", "10", "20", "30", "40")
# Bottom panel
par(mar=c(bmar, lmar, 0, 0))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimML)
axis(1, at=(0:6) * 200, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)
mtext("Infected density", 2, line = ylablineL, cex = cexlab)

lines(EBText[,"Time"], EBText[,"JnrI"], lwd = linelwd, col = "#0072B2")
lines(EBText[,"Time"], EBText[,"AnrI"], lwd = linelwd, col = "#D55E00")
# lines(EBText[,"Time"], EBText[,"TnrI"], lwd = linelwd, col = "black")

mtext("Time", 1, line = xlabline, cex = cexlab)

legend("topright", c("Juveniles", "Adults"), col = c("#0072B2", "#D55E00"), lty = c(1, 1), lwd = c(linelwd, linelwd), cex = cexleg, horiz = T)

##############################################################################################################################
##### Right column
##############################################################################################################################

# Top panel
par(mar=c(0, 0, tmar, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimTL)
axis(1, at=(0:6) * 200, labels = F, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)

lines(EBTinv[,"Time"], EBTinv[,"JnrS"], lwd = linelwd, col = "#0072B2")
lines(EBTinv[,"Time"], EBTinv[,"AnrS"], lwd = linelwd, col = "#D55E00")
# lines(EBTinv[,"Time"], EBTinv[,"TnrS"], lwd = linelwd, col = "black")

# Bottom panel
par(mar=c(bmar, 0, 0, rmar))
plot(NULL, NULL, type = "l",
     xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylimML)

axis(1, at=(0:6) * 200, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd)
axis(2, labels = F, cex.axis = cexaxs, lwd = 0, lwd.ticks = axislwd, las = 2)

lines(EBTinv[,"Time"], EBTinv[,"JnrI"], lwd = linelwd, col = "#0072B2")
lines(EBTinv[,"Time"], EBTinv[,"AnrI"], lwd = linelwd, col = "#D55E00")
# lines(EBTinv[,"Time"], EBTinv[,"TnrI"], lwd = linelwd, col = "black")

mtext("Time", 1, line = xlabline, cex = cexlab)

par(defplotpars)

if (ToPdf) {
  dev.off()
  system(paste("open ", fname))
}
