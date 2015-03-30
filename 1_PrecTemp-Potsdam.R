
# extreme precipitation and temperature
# Berry Boessenkool, 2014-2015, berry-b@gmx.de
# Script for data analysis and creation of graphics in main document
# Each section should run independent of the previous sections 
# (if the Rdata files are in the working directory)
# Sections are numbered roughly according to figure numbers in masters thesis

# Get the CRAN-releases of packages used in this script:
install.packages("lmomco") # for linear moments and distribution functions
install.packages("pbapply") # for progress bars
# get the current packages extremeStat and berryFunctions from github:
install.packages("devtools")
devtools::install_github("BerryBoessenkool/berryFunctions")
devtools::install_github("BerryBoessenkool/extremeStat") 
devtools::install_github("tpoisot/digitize")


# 0. format original data
# 1. P-T-Relationships from literature
# 2. Empirical precipitation quantiles in temperature bins
# 3. Distribution fitting, effect of cutoff point
# 4. Quantile ~ sample size
# 5. medians quantile ~ n
# 6. temperature dependent simulation
# 7. parametric quantiles
# 8 + 9. cutoff effect in section 3
# 10. gofProp effect
# 11. Analysis for each station



# 0. format original data ------------------------------------------------------
# Original data (PREC):
files <- dir("Original_data_Potsdam/KL03342", full=T)
# read first file:
datP <- read.fwf(files[1], widths=c(6,16,8,4), dec=",",
               col.names=c("stat","time","prec","n_5min_prec"), as.is=T)
# Add other files:
for(j in 2:length(files)) datP <- rbind(datP,
           read.fwf(files[j], widths=c(6,16,8,4), na.strings=c(" -9,99", "-99"),
              col.names=c("stat","time","prec","n_5min_prec"),dec=",", as.is=T))
# replace NAs with actual NAs
#datP[datP[,"prec"]==-9.99,"prec"] <- NA
#datP[datP[,"n_5min_prec"]==-99,"n_5min_prec"] <- NA
# remove completely-NA rows:
datP <- datP[ !is.na(datP$prec) & !is.na(datP$n_5min_prec) , ]

# check Station numbers:
if(length(unique(datP$stat))>1) stop("Several Station numbers in datP")
datP$time[duplicated(datP$time)]
datP <- datP[,-1]
# write to file:
write.table(datP, file="Original_data_Potsdam/Potsdam-Prec-formatted.txt",
            row.names=F, quote=F, sep="\t")
# actual time class:
datP$time <- as.POSIXct(strptime(datP$time, format="%d.%m.%Y %H:%M", tz="CET"))


# Original data (TEMP):
datT <- read.fwf("Original_data_Potsdam/kh03342.tem", widths=c(5,8,4),
                 col.names=c("stat","date","temp"), as.is=T)
# Temperature in °C:
datT$temp <- datT$temp/10
# reformat dates:
datT$date <- as.Date(as.character(datT$date), "%Y%m%d")
# NAs:
# datT[datT$temp==-99.9, "temp"] <- NA
if(length(unique(datT$stat))>1) stop("Several Station numbers in datT")
datT <- datT[,-1]
# write to disc:
write.table(datT, file="Original_data_Potsdam/Potsdam-Temp-formatted.txt",
            row.names=F, quote=F, sep="\t")


# merge prec and temp:
datP$date <- as.Date(datP$time)
dat <- merge(datP, datT, by="date", all=TRUE)
# Save as R-Objekt, to make starting from here faster:
save(dat, file="1_dat.Rdata")



# 1. P-T-Relationships from literature -----------------------------------------
library("digitize")  # x1,x2, y1,y2
setwd("PT_Literature_Plots")
files <- dir(pattern="*.JPG")
digitizePT <- function(fn){
 calpoints <- ReadAndCal(fn)
 write.table(as.data.frame(calpoints), paste0("calpoints/",fn,".txt"), row.names=F )
 logy <- readline("log y axis? y/n : ")
 vals <- readline("x1,x2, y1,y2, separated by comma: ")
 vals <- as.numeric(strsplit(vals, ",")[[1]])
 if(logy=="y") vals[3:4] <- log10(vals[3:4])
 write.table(vals, paste0("calvals/",fn,".txt"), row.names=F, col.names=F )
 data <- DigitData()
 dev.copy(png, paste0("digitalizedPlot/",fn, ".png")) ; dev.off()
 data2 <- Calibrate(data, calpoints, x1=vals[1], vals[2], vals[3], vals[4])
 if(logy=="y") data2$y <- 10^data2$y
 write.table(data2, paste0("PTnumbers/",fn,".txt"), row.names=F, col.names=F )
 }
digitizePT(files[24])


# read stuff back in:
#library(RColorBrewer)
source("1_PrecTemp-CC.R") # transp, col, probs, intervals, mid, cc_lines
setwd("PT_Literature_Plots")
files <- dir(pattern="*.JPG")
ptlit <- lapply(files, function(fn) read.table(paste0("PTnumbers/",fn,".txt")))
setwd("..")
names(ptlit) <- files

# scale mm/day intensities to hourly
#for(i in  1:12) ptlit[[i]]$V2 <- ptlit[[i]]$V2/24   # Berg_seas_2_01-12
for(i in 20:22) ptlit[[i]]$V2 <- ptlit[[i]]$V2/24   # Utsumi_3b-d
#
litfac <- as.factor(substr(files, 1, 7))
litcol <- brewer.pal(8, "Set1")
litcol <- litcol[c(6, 2:5, 1, 7:8)]
#display.brewer.all(n=NULL, type="qual", select=NULL, exact.n=TRUE)
#dput(levels(litfac))
litnames <- c("Berg (2009)", "Berg (2013a)", "Berg (2013)", "Hardwick Jones (2010)",
             "Lenderink (2008)", "Mishra (2012)", "Utsumi (2011)", "Westra (2014)")
monthcol <- colorRampPalette(c("blue","grey","red","grey","blue"))(12)
#ccstarts <- lapply(1:8, function(x) NA)
#ccstarts[[2]] <- c(9)
#ccstarts[[3]] <- c(2,8)
#ccstarts[[4]] <- c(7.5,21)
#ccstarts[[7]] <- c(2.5,7.5)
#ccstarts[[8]] <- c(4)

pdf("figure/1_PT_Literature.pdf", height=4, pointsize=9)
par(mfrow=c(2,3), mar=c(2,2,0.2,0.2), oma=c(1.5,1.5,0.2,0), mgp=c(1.6,0.7,0), cex=1)
for(i in c(1:4,7:8) ){       #  1){#
  inpanel <- which( as.numeric(litfac)==i )
  trange <- range(do.call(rbind, ptlit[inpanel])[,1], na.rm=TRUE)
  prange <- range(do.call(rbind, ptlit[inpanel])[,2], na.rm=TRUE)
  if(prange[1] < 1) prange[1] <- 1
  plot(1, log="y", xlab="", ylab="", ylim=prange, xlim=trange, type="n", las=1)
  # cc_lines(ccstarts[[i]], maincc=T)
  if(i!=1) cc_lines(NA, mainargs=c(col=2, lty=2))
  axis(2, logVals()$all, labels=FALSE)
  if(i!=1) for(k in inpanel) lines(ptlit[[k]])
  if(i==1) text(-10, 48, "mm/day", adj=c(0,1))
  if(i==1) for(k in inpanel) lines(ptlit[[k]], col=monthcol[k])
  if(i==1) colPointsLegend(z=1:12, x1=20, y1=45, y2=25, nbins=12, col=monthcol,
                           title="month of year", dens=F)
  }
  title(ylab="Precipitation  [mm/h]", line=0.1, outer=TRUE, cex=1)
  title(xlab="Temperature  [°C]",     line=0.1, outer=TRUE, cex=1)
dev.off()



# 2. Empirical precipitation quantiles in temperature bins ---------------------
rm(list=ls()) # clear workspace to ensure every section is reproducible
library(berryFunctions) # for quantileMean, logAxis
load("1_dat.Rdata")
source("1_PrecTemp-CC.R") # transp, col, probs, intervals, mid, cc_lines

# Run only once:
dat <- dat[dat$prec > 0.5,]    # 7667,  25840 with !=0
# quantiles (average of all 9 methods included in R) for each bin:
binQ <- sapply( tapply( dat$prec,  cut(dat$temp, intervals),
         quantileMean, probs=probs) , I)
# number of values per bin:
binN <- tapply( dat$prec,  cut(dat$temp, intervals), length)
## control:
## sum(binN, na.rm=T)  ;  sum(datP$prec!=0)
## 6 values are missing now (4 to 17 at other Stations)
## sum(is.na(dat$temp)) # =6
# save results:
save(dat, binN, binQ, intervals, mid, file="2_dat_binN_binQ_intervals_mid.Rdata")
# End run only once

load("2_dat_binN_binQ_intervals_mid.Rdata")

pdf("figure/2_empiricalquantile_log.pdf", height=4, width=3.5, pointsize=9)
  par(mar=c(3,2.8,0.2,0.4), mgp=c(1.8,0.5,0))
  plot(dat$temp, dat$prec, ylab="Precipitation [mm/h]", xlab="", yaxt="n", xaxt="n",
       yaxs="i", ylim=c(1.6, 45), type="n", xlim=c(0,30), xaxs="i", log="y")##, cex.axis=0.8, tcl=-0.2)
  title(xlab="Daily average temperature  [°C]", mgp=c(1.7,1,0) )
  logAxis(side=2, mgp=c(2.3,0.7,0) )
  axis(1,c(0,10,20,30))
  points(dat$temp, dat$prec, col=transp, pch=16, cex=1)
  abline(v=intervals, col=8)
  cc_lines(c(1.5,3.5))
  ## plot quantiles only for bins with more than 100 values:
  Sel <- binN > 100  # n-based Selection
  for(q in 1:4) points(mid[Sel], binQ[q,Sel], col=col[q], pch=q+14, type="o",
                       cex=1, lwd=2)
  # Legend, Annotations
  lg <- legend("topleft", c("99.99 %Q","99.9","99","90", "VPsat [hPa]", "CC-rate"),
         pch=c(4:1+14, NA, NA), col=c(rev(col),1,1), bg="white",
         lty=c(rep(1,5),3), inset=c(0, 0), cex=1)
  #legend(lg$rect$left, lg$rect$top - lg$rect$h,
  #      "Number of data points\nper temperature bin",
  #      pch=110, cex=0.7, bg="white", adj=c(0,0.3)) # c(110, NA)
  # write number of values per bin
  textField(mid, c(1.8, 1.95), binN, cex=0.8, quiet=T)
  box()
dev.off()



# 3. Distribution fitting, effect of cutoff point ------------------------------
rm(list=ls()) # clear workspace to ensure every section is reproducible
load("1_dat.Rdata")
source("1_PrecTemp-CC.R") # transp, col, probs, intervals, mid, cc_lines
library(berryFunctions) # quantileMean, logAxis, functions in extremeStat
library(extremeStat) # distLfit, distLplot
library(lmomco) # qlmomco
library(pbapply)
nbest <- 8

logPrec <- log10(dat[dat$prec > 0.5 , "prec"])
dlf <- distLfit(logPrec, gofProp=0.1, plot=F)

pdf("figure/3_FitDist.pdf", height=3, pointsize=9)
  par(mfrow=c(1,3), mar=c(3,3.7,0.4,0), mgp=c(1.7,0.6,0), cex=1)
  # select best fits only by RMSE of upper 5% of data
  distLplot(dlf, log=TRUE, nbest=nbest, main="", xlab="Precipitation [mm/h]",
             ylab="", legargs=list(bg="white", cex=1), percentargs=c(lwd=2, lty=2),
             logargs=c(allticks=TRUE, lcol=NA), percentline=FALSE )
  title(ylab="Probability Density Function (PDF)", mgp=c(2.3, 0.7, 0))
  rect(xleft=0.6, ybottom=0, xright=2, ytop=0.09, border="purple", col=addAlpha("purple"))
  box()
  # focus on tail
  par(mar=c(3,3,0.4,0.7))
  tquan <- distLquantile(dlf=dlf, probs=0.99, plot=TRUE, breaks=20, xlim=c(0.6,2),
           ylim=c(0,0.09), nbest=nbest, log=TRUE, main="", percentline=FALSE,
           legargs=list(bg="white", cex=1), percentargs=c(lwd=2, lty=2),
           linargs=c(lwd=2), logargs=c(allticks=TRUE, lcol=NA), ylab="", xlab="")
  title(xlab="Purple box precipitation range")
  box()#col="purple")
  # truncated data
  distLquantile(dlf=dlf, truncate=0.9, probs=0.99, plot=TRUE, log=TRUE, ylab="",
           nbest=nbest, legargs=list(bg="white", cex=1), xlab="Truncated to top 10%",
           linargs=c(lwd=2), logargs=c(allticks=TRUE, lcol=NA), main="")
  box()
dev.off()


dlf1 <- distLfit(log10(dat[dat$prec > 0.1 , "prec"]), plot=F)
dlf2 <- distLfit(log10(dat[dat$prec > 2,    "prec"]), plot=F)

pdf("figure/3_cutoffeffect.pdf", height=2, pointsize=9)
  par(mfrow=c(1,3), mar=c(3,3.7,0.4,0), mgp=c(1.7,0.6,0), cex=1)
  distLplot(dlf1, log=TRUE, nbest=nbest,main="", xlab="Precipitation [mm/h]",
                  legargs=list(bg="white", cex=0.8), percentargs=c(lwd=2), ylab="" )
  title(ylab="Probability Density Function (PDF)", mgp=c(2.3, 0.7, 0))
  rect(xleft=0.6, ybottom=0, xright=2, ytop=0.25, border="purple", col=addAlpha("purple"))
  box()
  # focus on tail
  par(mar=c(3,3,0.4,0.7))
  distLplot(dlf1, breaks=20, xlim=c(0.6,2), ylim=c(0,0.1), nbest=nbest, log=TRUE,
            main="", legargs=list(bg="white", cex=0.8), ylab="", xlab="")
  title(xlab=paste("cutoff point (prec > cutoff): 0.1 mm/h") )
  # different cutoff
  distLplot(dlf2, breaks=10, xlim=c(0.6,2), ylim=c(0,0.1), nbest=nbest, log=TRUE,
            main="", legargs=list(bg="white", cex=0.8), ylab="", xlab="")
  title(xlab=paste("cutoff point (prec > cutoff): 2 mm/h") )
dev.off()



# 4. Quantile ~ sample size ----------------------------------------------------
rm(list=ls()) # clear workspace to ensure every section is reproducible
library(berryFunctions) # for logAxis, quantileBands
library(extremeStat) # for distLfit, distLquantile
library(pbapply) # for progress bar in sapply loop
library(RColorBrewer) # for brewer.pal
library(lmomco) # for rlmomco, qlmomco
load("1_dat.Rdata")
probs <- c(0.9, 0.95, 0.99, 0.999, 0.9999) # different from other probs!
logPrec <- log10(dat[dat$prec>0.5, "prec"])
dlf <- distLfit(logPrec, plot=F, ks=F)
samplesize <- 25:1000 # samplesize <- seq(25, 1000, length=20)


# execute this only once - time consuming! 200 runs: 1.5 - 3.4 h per dist type
system.time(qn_wak <- pbreplicate(n=500, sapply(samplesize, function(nn) # 500: 8.4 h
   distLquantile(rlmomco(nn, para=dlf$parameter$wak), probs=probs, quiet=TRUE,
                 order=FALSE, truncate=0.8, type="wak")  ), simplify = FALSE))
save(qn_wak, file="5_qn_wak_500.Rdata")

system.time(qn_wei <- pbreplicate(n=200, sapply(samplesize, function(nn)
   distLquantile(rlmomco(nn, para=dlf$parameter$wei), probs=probs, quiet=TRUE,
   order=FALSE, truncate=0.8, type="wei")), simplify = FALSE))
save(qn_wei, file="5_qn_wei_200.Rdata")

system.time(qn_kap <- pbreplicate(n=200, sapply(samplesize, function(nn)
   distLquantile(rlmomco(nn, para=dlf$parameter$kap), probs=probs, quiet=TRUE,
   order=FALSE, truncate=0.8, type="kap")), simplify = FALSE))
save(qn_kap, file="5_qn_kap_200.Rdata")

system.time(qn_realdat <- pbreplicate(n=200, sapply(samplesize, function(nn) # 3.8 h
   distLquantile(sample(logPrec, nn), probs=probs, quiet=TRUE,
   order=FALSE, truncate=0.8, type=c("wak","wei","kap","gpa"))), simplify = FALSE))
save(qn_realdat, file="5_qn_realdat_200.Rdata")


# continue here:
load("5_qn_wak_200.Rdata")
load("5_qn_kap_200.Rdata")
load("5_qn_wei_200.Rdata")
load("5_qn_realdat_200.Rdata")

plotqn <- function(
  qn=qn_wak,         # object
  dist="wak",        # distribution abbreviation
  p=4,               # integer corresponding to percentile in probs
  ylim=c(0.3, 2),  # y axis range
  smooth=15          # smoothing of quantileBands
  )
  {
  par(mar=c(3,3.6,0.2,0.2), mgp=c(1.8,0.7,0), mfrow=c(1,2))
  # Empirical quantiles
  plot(samplesize, qn[[1]][p,], type="n", ylim=ylim, xlim=lim0(980), las=1, ylab="",
      xlab="sample size n", yaxt="n", main="")
  title(ylab=paste0("random sample ", probs[p]*100, "% quantile [mm/h]"), mgp=c(2.5,1,0) )
  logAxis(side=2)
  mat <- sapply(qn, function(x) x[p+5,])
  quantileBands(t(mat), smooth=smooth, probs=0:10/10, meanargs=list(col=2), border=NA,
                x=samplesize, col=brewer.pal(9,"BuGn")[3:7], txi=NA, add=TRUE, na.rm=TRUE)
  lines(samplesize, movAv(apply(mat, 1, median, na.rm=TRUE), smooth))
  tquan <- qlmomco(probs[p], para=dlf$parameter[[dist]])
  abline(h=tquan, lwd=1, col="purple")
  box()
  # Legend
  legend("bottomright", c(paste("smoothing bandwidth =", smooth),
    paste("distribution quantile:", round(10^tquan, 1), "mm/h"), "median of 200 simulations",
    "simulation mean (not smoothed)"), bg="white", lwd=c(NA,1,1,1), col=c(NA, "purple", 1, 2), cex=1)
  # Parametrical quantiles
  plot(samplesize, qn[[1]][p,], type="n", ylim=ylim, xlim=lim0(980), las=1, ylab="",
      xlab="sample size n", yaxt="n", main="")
  title(ylab="parametrical quantile of random sample", mgp=c(2,1,0) )
  logAxis(side=2)
  mat <- sapply(qn, function(x) x[p,])
  quantileBands(t(mat), smooth=smooth, probs=0:10/10, meanargs=list(col=2), border=NA,
                x=samplesize, col=brewer.pal(9,"BuGn")[3:7], txi=NA, add=TRUE, na.rm=TRUE)
  lines(samplesize, movAv(apply(mat, 1, median, na.rm=TRUE), smooth))
  abline(h=tquan, lwd=1, col="purple")
  box()
  }

# Actual plots
pdf(paste0("figure/4_quantile_n_wak.pdf"), height=2.9, pointsize=9)
  for(i in 1:5) plotqn(qn_wak, "wak", i, ylim=c(0.3, 2), smooth=15)
dev.off()

pdf(paste0("figure/4_quantile_n_kap.pdf"), height=4)
  for(i in 1:5) plotqn(qn_kap, "kap", i, ylim=c(0.3, 2), smooth=15)
dev.off()

pdf(paste0("figure/4_quantile_n_wei.pdf"), height=4)
  for(i in 1:5) plotqn(qn_wei, "wei", i, ylim=c(0.3, 2), smooth=15)
dev.off()


# samples from observed intensities (real data)
pdf(paste0("figure/4_quantile_n_realdata.pdf"), height=2.5, width=3.5, pointsize=9)
  type=c(2,1,3:5)                # distribution type: wak, wei, kap, gpa, empirical
  p=4                            # integer corresponding to percentile in probs
  ylim=c(1,1.5)                  # y axis range
  smooth=15                      # smoothing of quantileBands
  col=col=brewer.pal(5, "Dark2") # colors for bands
  bandrange=c(0.45, 0.55)        # limits of band range
  #
  par(mar=c(3,3,0.2,0.2), mgp=c(1.8,0.7,0) )
  dn <- c("wak", "wei", "kap", "gpa", "empirical")
  plot(samplesize, qn_realdat[[1]][p,], type="n", ylim=ylim, xlim=lim0(980), las=1, ylab="",
      xlab="sample size n", yaxt="n", main="")
  title(ylab=paste("subsample", probs[p]*100, "% quantile [mm/h]"), mgp=c(1.9,1,0) )
  logAxis(side=2)
  axis(2, log10(30), 30, las=1)
  for(typ in type) # wak, wei, kap, gpa, empirical
  {
  mat <- sapply(qn_realdat, function(x) x[5*typ-5+p,])
  top <- movAv(apply(mat, 1, quantileMean, probs=bandrange[2], na.rm=T), smooth)
  bot <- movAv(apply(mat, 1, quantileMean, probs=bandrange[1], na.rm=T), smooth)
  sel <- !is.na(top)
  polygon(x=c(samplesize[sel], rev(samplesize[sel])), y=c(bot[sel], rev(top[sel])),
          col=addAlpha(col[typ]), border=NA)
  lines(samplesize, movAv(apply(mat, 1, median, na.rm=TRUE), smooth), col=col[typ])
  }
  # legend:
  legend("bottomright", dn[type], lty=1, col=col[type], bg="white")
  abline(h=quantileMean(logPrec, probs=probs[p]))
dev.off()


# 5. medians quantile ~ n ------------------------------------------------------
# continue from section 4
medqn <- function(qn, p)
   {
   out <- movAv(apply(sapply(qn, function(x) x[p,]), MARGIN=1, median), 15)
   out[out==0] <- NA
   out
   }

pdf(paste0("figure/5_quantile_n_median.pdf"), height=2.1, pointsize=9)
  col <- c(NA,NA,"forestgreen","darkblue","red")
  layout(matrix(1:2, ncol=2), widths=c(6,4.8))
  par(mar=c(3,3.2,0.2,2.0), mgp=c(1.8,0.8,0), lend=1)
  # empirical quantiles
  plot(1, type="n", ylim=c(0.95, 1.58), xlim=lim0(980), las=1, ylab="",
      xlab="sample size n", yaxt="n")
  title(ylab="Median quantile of simulations", mgp=c(2.1,1,0) )
  logAxis(side=2)
  axis(2, log10(30), 30, las=1)
  for(i in 4:3) lines(samplesize, medqn(qn_wei, i+5), col=col[i], lty=1, lwd=2)
  for(i in 4:3) lines(samplesize, medqn(qn_wak, i+5), col=col[i], lty=1, lwd=2)
  for(i in 4:3) lines(samplesize, medqn(qn_kap, i+5), col=col[i], lty=1, lwd=2)
  q_wei <- qlmomco(probs[3:4], para=dlf$parameter[["wei"]])
  q_wak <- qlmomco(probs[3:4], para=dlf$parameter[["wak"]])
  q_kap <- qlmomco(probs[3:4], para=dlf$parameter[["kap"]])
  points(rep(970,2), q_wei, col=col[3:4], pch=16)
  points(rep(970,2), q_wak, col=col[3:4], pch=16)
  points(rep(970,2), q_kap, col=col[3:4], pch=16)
  text(rep(990,2), c(1.10, 1.52), rep("wei",2), adj=0, xpd=TRUE, col=col[3:4])
  text(rep(990,2), c(1.05, 1.37), rep("wak",2), adj=0, xpd=TRUE, col=col[3:4])
  text(rep(990,2), c(1.00, 1.25), rep("kap",2), adj=0, xpd=TRUE, col=col[3:4])
  text(200, 1.44, "99.9 % Quantile", col=col[4])
  text(500, 0.97, "99 % Quantile",   col=col[3])
  # theoretical quantiles
  par(mar=c(3,1.6,0.2,0.2) )
  plot(1, type="n", ylim=c(0.95, 1.58), xlim=lim0(980), las=1, ylab="",
      xlab="sample size n", yaxt="n")
  logAxis(side=2)
  for(i in 3:4) lines(samplesize, medqn(qn_wei, i), col=col[i], lty=1, lwd=2)
  for(i in 3:4) lines(samplesize, medqn(qn_wak, i), col=col[i], lty=1, lwd=2)
  for(i in 3:4) lines(samplesize, medqn(qn_kap, i), col=col[i], lty=1, lwd=2)
dev.off()



# 6. temperature dependent simulation ------------------------------------------
rm(list=ls()) # clear workspace to ensure every section is reproducible
load("1_dat.Rdata")
source("1_PrecTemp-CC.R") # transp, col, probs, intervals, mid, cc_lines
library("lmomco")
library("pbapply")
library("berryFunctions")


# Only once:
load("Lange-Zeitreihen-DWD/PT_data.Rdata")
PT <- PT_data[[1]]
for(i in 2:12) PT <- rbind(PT, PT_data[[i]])
PT <- PT[PT$prec > 0.5,]
p <- parwak(lmoms(log10(PT$prec)))
p$para
# sample sizes:
intervals <- -10:40
binN <- tapply(PT$prec, cut(PT$temp, intervals), length)
intervals <- intervals[binN>=25 & !is.na(binN)]
binN <- binN[binN>=25 & !is.na(binN)]
mid <- intervals[-length(intervals)]+0.5
n <- length(mid)
# calculate parameters of the Wakeby distribution across temperature bins:
wakpar <- sapply( tapply( PT$prec,  cut(PT$temp, intervals),
         function(x) parwak(lmoms(log10(x)))$para), I)
# Quantile across temperature bins:
tq999 <- sapply(1:n, function(i) quawak(0.999, list(type="wak", para=wakpar[,i])))
# empirical Quantile:
eq999 <- tapply(log10(PT$prec), cut(PT$temp, intervals), quantileMean, probs=0.999)
# parameters of Wakeby distribution increasing at CC rate:
linpar <- function(temp)
   {
   para <- sapply(1:5, function(i)
         predict(lm(I(wakpar[i,])~mid), newdata=data.frame(mid=temp)))
   list(type="wak", para=para)
   }
# 99.9% Quantile of distribution
qlin <- sapply(mid, function(t) quawak(0.999, linpar(t)))
# Save stuff:
save(wakpar, PT, p, n, mid, binN, intervals, tq999, eq999, linpar, qlin, file="4_wakpar_etc.Rdata")

# Empirical + theoretical quantile from random numbers sample of original size
system.time(
verif <- pbreplicate(n=500, sapply(1:length(mid), function(i)
                     {rand <- rlmomco(n=binN[i], para=linpar(mid[i]))
                     randtrunc <- sort(rand)[ -1:-(0.8*length(rand)) ]
                     equant <- quantileMean(rand, prob=0.999)
                     tquant <- quawak(0.995, parwak(lmoms(randtrunc), checklmom=F))
                     c(equant, tquant)}))
) # 40 min
save(verif, file="7_tempdep_500.Rdata")


# start here later:
load("4_wakpar_etc.Rdata")
load("7_tempdep_500.Rdata")

pdf("figure/6_tempdep_500.pdf", height=3, pointsize=9)
  layout(matrix(c(1:5, rep(6,5)), ncol=2), width=c(4,6))
  # Show parameter of the distribution across temperature bins:
  par(mar=c(0,2.8,1.5,1), oma=c(3,0,0,0), mgp=c(3,0.5,0), cex=1)
  for(i in 1:5)
    {
    plot(mid, wakpar[i,], type="l", las=1, xaxt="n", yaxt="n")
    title(main=names(p$para)[i], cex.main=1, font.main=1, line=0.2)
    axis(2, pretty2(par("usr")[3:4], n=2, min.n=2, force=T), las=1)
    #abline(lm(I(wakpar[i,])~mid), col=2)
    #abline(lm(I(wakpar[i,2:25])~mid2), col=3)
    lines(mid, linpar(mid)$para[,i], col=2)
    if(i==5) {axis(1) ; title(xlab="Temperature  [°C]", mgp=c(1.5,1,0), xpd=NA)}
    }
# Show 99.9% quantiles from random samples:
  par(mar=c(0,3.5,1,1), mgp=c(1.8,0.5,0))
  plot(mid, qlin, type="l", col=2, lwd=3, ylim=c(0.6, 1.7), yaxt="n", xaxs="i",
       ylab="Precipitation  [mm/h]", xlab="")
  title(xlab="Temperature  [°C]", mgp=c(1.5,1,0), xpd=NA )
  logAxis(2)
  ### Empirical from samples
  ##textpos <- rep(22,11); textpos[c(3,4,6,8,9)] <- NA
  ##quantileBands(t(verif), x=mid, add=TRUE, txi=textpos,
  ##               col=rgb(0,0,1, alpha=1:5/7), probs=0:10/10)
  # Empirical:
  equant <- apply(verif[1:n*2-1,], MARGIN=1, quantileMean, probs=c(0.1, 0.5, 0.9))
  polygon(x=c(mid, rev(mid)), y=c(equant[1,], rev(equant[3,])), col=addAlpha(4), border=NA)
  # Theoretical:
  tquant <- apply(verif[1:n*2,], MARGIN=1, quantileMean, probs=c(0.1, 0.5, 0.9))
  polygon(x=c(mid, rev(mid)), y=c(tquant[1,], rev(tquant[3,])), col=addAlpha(3), border=NA)
  # CC rate:
  temp <- seq(-8,32, len=100)
  lines(temp, log10( 6.1094*exp(17.625*temp/(temp+243.04))), col=1)
  # Distribution theoretical quantiles again (on top)
  lines(mid, qlin, col=2, lwd=3)
  # Simulation medians:
  lines(mid, equant[2,], col=4)
  lines(mid, tquant[2,], col="darkgreen")
  # write sample sizes:
  text(mid, rep(c(0.58, 0.61, 0.64), len=length(mid)), binN, cex=0.7, xpd=T)
  box()
dev.off()



# 7. parametric quantiles -----------------------------------------------------
rm(list=ls()) # clear workspace to ensure every section is reproducible
load("1_dat.Rdata")
source("1_PrecTemp-CC.R") # transp, col, probs, intervals, mid, cc_lines
library(berryFunctions) # l2df
library(extremeStat) # distLquantile
library(lmomco) # for qlmomco, par*** or lmom2par, lmoms
dat <- dat[dat$prec>0.5,]
#intervals <- intervals[2:14] # last bin has only 1 value, first one 3
#mid <- intervals[-1] - diff(intervals)/2
#cex <- 0.7

binN <- tapply( dat$prec, cut(dat$temp, intervals), length)
binQ <- tapply( log10(dat$prec),  cut(dat$temp, intervals),
     distLquantile, probs=probs, truncate=0.8, quiet=T, type=c("wei","wak","kap"))
binQ <- sapply(binQ, I)

pdf("figure/7_theoreticalquantiles_compare.pdf", height=3.5, pointsize=9)
  col <- c("orange","purple", "forestgreen", "black")
  ##layout(matrix(1:2, ncol=2), widths=c(5,5))
  par(mar=c(0.5,2.5,0,0), mgp=c(3,0.8,0), oma=c(2.5,1.0,0.3,0.2), mfrow=c(2,2), cex=1)
  for(q in 1:4)
  {
  plot(dat$temp, dat$prec, ylab="", xlab="", yaxt="n", #yaxs="i",
       ylim=10^c(min(binQ[4*1:4-4+q,], na.rm=T),pmin(max(binQ[4*1:4-4+q,], na.rm=T), 2.7)),
       type="n", log="y", xaxt=if(q %in% 1:2) "n" else "s")
  logAxis(side=2, mgp=c(2.3,0.7,0) )
  points(dat$temp, dat$prec, col=transp, pch=16, cex=1)
  cc_lines(mainargs=c(lty=2, lwd=2, col=2), startpunkte=1000)
  for(i in 1:4) lines(mid, 10^binQ[4*i-4+q,], col=col[i])
  # legend:
  legend("topleft", paste(probs[q]*100, "%Q"), bty="n")
  if(q==2)legend("bottomright", c("Weibull","Wakeby","Kappa","Empirical"), col=col,
         lty=1, bg="white", cex=1, inset=c(0.1, 0))
  if(q==1) textField(mid, c(1.8, 2.2), binN, cex=0.8, mar=0.3, quiet=T)    #  0.5   c(0.4, 1.0)
  }
title(ylab="Precipitation  [mm/h]", line=0, outer=T)
title(xlab="Daily average temperature  [°C]", line=1.2, outer=T)
dev.off()



# 8 + 9. cutoff effect in section 3


# 10. gofProp effect -----------------------------------------------------------
rm(list=ls()) # clear workspace to ensure every section is reproducible
load("1_dat.Rdata")
library(extremeStat) # distLfit, distLplot
library(berryFunctions) # rainbow2, quantileMean
library(lmomco) # qlmomco

# only use entries with prec >0.5 mm/h:
LogPrec <- log10(dat[dat$prec > 0.5 , "prec"]  )
# Fit distributions (via linear moments):
dlf <- distLfit(LogPrec, plot=FALSE, gofProp=0.05, ks=FALSE)
# effect of gofProp with prec >0.5  (see oldcode.R for >3mm)
gps <- c(20:2*5, 9:1, 0.5)

# Execute only once:
rmse_0.5 <- pbapply::pblapply(gps/100, function(gp)
           distLfit(LogPrec, plot=F, gofProp=gp, progbars=F, time=F, quiet=T)$gof)
save(rmse_0.5, file="rmse_0.5.Rdata")

# continue here
load("rmse_0.5.Rdata")
ranks_0.5 <- sapply(rmse_0.5, rownames)
dn <- unique(as.vector(sapply(rmse_0.5, rownames)))
set.seed(3); mycol=sample(rainbow2(16))

pdf("figure/10_fitdist_gofProp-effect.pdf", height=3, pointsize=9)
  layout(matrix(1:2, ncol=2), widths=c(6, 5))
  par(mar=c(3.5,4,0.5,0.5), cex=1)
  plot(1, type="n", ylim=c(0.0005,0.4), xlim=c(110,-10), log="y", mgp=c(3,0.8,0),
      xlab="", ylab="", las=1, yaxt="n")
  title(ylab="RMSE of distributions", mgp=c(2.9,1,0))
  logAxis(2, mgp=c(3,0.6,0), base=1)
  for(i in 16:1)
  lines(gps, sapply(rmse_0.5, function(x) x[dn[i],"RMSE"]), col=mycol[i])
  text(rep(100,16), rmse_0.5[[1]]$RMSE,         dn, cex=0.7, adj= 1.05, col=mycol)
  text(rep(  0,16), rmse_0.5[[29]][dn, "RMSE"], dn, cex=0.7, adj=-0.05, col=mycol)
  # set.seed(1); rs <- sample(3:27, 17) # random location selection
  rs <- c(19,17,16,14,13,13,8,7,16,13,11,18,12,15,14,18)
  for(i in 1:16)
  textField(gps[rs[i]], rmse_0.5[[rs[i]]][dn[i], "RMSE"], dn[i], cex=0.7,
             col=mycol[i], margin=0.1, field="round", fill="white")
  #
  par(mar=c(3.5,1.5,0.5,0.5))
  plot(1, type="n", ylim=c(17,1), xlim=c(110,-10), xaxt="s", yaxt="n", mgp=c(3,0.8,0),
      xlab="")
  title(ylab="Rank of distributions", mgp=c(0.5,1,0))
  ranks <- apply(ranks_0.5, 2, function(x) match(dn,x))
  for(i in 16:1) lines(gps, ranks[i,], col=mycol[i])
  text(rep(100,16), ranks[, 1], dn, cex=0.7, adj= 1.05, col=mycol)
  textField(rep( 50,17), ranks[,11], dn, cex=0.7, adj= 0.5 , col=mycol, margin=0, field="round")
  text(rep(  1,16), ranks[,29], dn, cex=0.7, adj=-0.05, col=mycol)
  mtext("Proportion  [%]  of highest data used to calculate RMSE as indicator of goodness of fit.",
      side=1, line=-1.5, outer=TRUE)
  #mtext("Focus on tail influences goodness of fit  and distribution selection",
  #    side=3, line=-2, outer=TRUE, cex=1.2, font=2)
dev.off()






# 11. Analysis for each station ------------------------------------------------
load("Lange-Zeitreihen-DWD/PT_data.Rdata")
library("lmomco")
library("berryFunctions")
intervals <- seq(-8,32, by=2) 
mid <- intervals[-1] - diff(intervals)/2
transp <- rgb(.3, .3, .3, alpha=0.3)
col <- c("green3","cyan2","darkblue","red")
col2 <- c("orange","purple", "forestgreen")
cc_lines <- function()  {  temp <- seq(-8,32, len=100)
                           cc <- function(temp) 6.1094*exp(17.625*temp/(temp+243.04))
                           changerate <- diff(cc(-8:33))/cc(-8:32) + 1
                           startpunkte <- if(par("ylog")) 10^(seq( -1.95, 2, len=13)) else 1.1^(seq(-15,40, len=14))
                           for(startp in startpunkte)lines(-8:33, cumprod(c(startp, changerate)), lty=3)
                           lines(temp, cc(temp) )  }

pdf("figure/9_allstations_tquant.pdf")
par(mar=c(1.8,2,0,0), mgp=c(3,0.5,0), cex.main=1, oma=c(1.1,1.3,2,0.1), mfrow=c(2,2), lend=1)
for(i in 1:14)
{
message(i)
PT <- PT_data[[i]]
PT <- PT[PT$prec>0.5,]
# distribution quantiles per temperature bin:
binQt <- tapply(PT$prec, cut(PT$temp, intervals), function(x){
  if(length(x)>5)
    sapply(c("kap", "wak", "wei"), function(d){
      param <- lmom2par(lmoms(log10(x)), type=d)
      if(is.null(param)) quants <- rep(NA, 4) else
         quants <- 10^qlmomco(f=c(0.9, 0.99, 0.999, 0.9999), param)
      if(length(quants)==0) quants <- rep(NA, 4)
      quants
      })
  else NULL
  })
binQt <- lapply(binQt, function(x) if(is.null(x)) matrix(NA, ncol=3, nrow=4) else x)
binQe <- tapply(PT$prec, cut(PT$temp, intervals),
         quantileMean, probs=c(0.9, 0.99, 0.999, 0.9999))
binQe <- lapply(binQe, function(x) if(is.null(x)) rep(NA,4) else x)
binN <- tapply(PT$prec, cut(PT$temp, intervals), length)
binN[is.na(binN)] <- 0
binQe <- lapply(1:length(binQe), function(j) if(binN[j]<5) rep(NA,4) else binQe[[j]])
# plot
for(q in 1:4)
  {
  plot(PT$temp, PT$prec, ylab="", xlab="", yaxt="n", yaxs="i", ylim=c(1, 99),
       type="n", log="y", xlim=c(-6,31))
  logAxis(side=2, mgp=c(2.3,0.7,0) )
  points(PT$temp, PT$prec, col=transp, pch=16)
  cc_lines()
  lines(mid, sapply(binQe, "[", q), lwd=2)#, type="o", pch=q+14, col=col[q])
  for(d in 1:3)
     lines(mid, sapply(binQt, "[", q, d), col=col2[d], lwd=2)#, type="o", pch=q+14, cex=cex
  legend("topleft", c("90.00 %", "99.00 %", "99.90 %", "99.99 %")[q], border=NA, bg="white")
  if(q==1) legend("topright", c("Kappa","Wakeby","Weibull","Empirical"), col=c(col2, 1),
       lty=1, bg="white", lwd=2)
  if(q==4) text(mid, rep(c(1.2, 1.5), len=length(mid)), binN, cex=0.7)
  }
title(main=names(PT_data)[i], outer=T)
title(ylab="Precipitation  [mm/h]", mgp=c(0.1,1,0), outer=T)
title(xlab="Daily average temperature  [°C]", mgp=c(-0.1,-1,0), outer=T )
}
dev.off()

