# Analysis of temperature-dependency of rainfall intensity
# With 142 freely available hourly time series (>15 years) from DWD FTP-server
# Berry Boessenkool, 2016
browseURL("http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183")

# Each section should run independently in a clean R session
# Frequently needed auxiliary functions are in Code_aid.R
# some code is commented out to avoid accidental calling (e.g. if time consuming)
# computing times are noted in minutes (e.g. # 5 min)
# Historically, section 3 came before section 2, but this is organized by content

# Contents ---------------------------------------------------------------------

# 0. Packages

# 1. Data 
# 1.1. Select DWD stations
# 1.2. Meta data weather stations
# 1.3. Download DWD station data
# 1.4. Dewpoint temperature
# 1.5. Dew point temperature of previous 5 hours

# 2. Sample size dependency
# 2.1. SSD all distributions
# 2.2. Distribution weights
# 2.3. SSD GPD 
# 2.4. Truncation dependency 
# 2.5. tempdep Wakeby distribution
# 2.6. Sample size bias lin vs log

# 3. Hourly Prec-Temp relationship
# 3.1. Raw data visualisation
# 3.2. PT-quantiles computation
# 3.3. PT-quantiles visualization
# 3.4. PTQ per station


# 0. Packages ------------------------------------------------------------------
if(FALSE){ # You need to download and install the packages only once
packinst <- function(n) if(!requireNamespace(n, quietly=TRUE)) install.packages(n)
sapply(c("berryFunctions", "extremeStat", "pblapply", "maps", "gplots", "gtools",
         "mapdata", "OSMscale", "RCurl"), packinst) 
berryFunctions::instGit("brry/berryFunctions")# must be version >= 1.13.10(2017-01-03)
berryFunctions::instGit("brry/extremeStat")   # must be version >= 1.2.12 (2017-01-04)
berryFunctions::instGit("brry/rdwd")  # not yet on CRAN must be >= 0.5.4  (2016-11-25)
berryFunctions::instGit("brry/OSMscale")      # must be version >= 0.3.12 (2016-11-24)
}
library(berryFunctions); library(pbapply) # potentially needed in every section



# 1. Data ----------------------------------------------------------------------

# 1.1. Select DWD stations -----------------------------------------------------

# Filenames for suitable stations:
library("rdwd")
#files <- selectDWD(res="hourly", var=c("air_temperature", "precipitation"),
#                   per="h", current=TRUE)
#save(files, file="dataprep/files.Rdata")
load("dataprep/files.Rdata")

# select stations with file available for both prec and temp:
temp_id <- as.integer(substr(files[[1]], 109,113))
prec_id <- as.integer(substr(files[[2]], 107,111))
files[[1]] <- files[[1]][temp_id %in% prec_id]
files[[2]] <- files[[2]][prec_id %in% temp_id]

# select stations with the longest time records (according to file name):
T_dur <- as.numeric(substr(files[[1]],124,131))/1e4 -
         as.numeric(substr(files[[1]],115,122))/1e4 
P_dur <- as.numeric(substr(files[[2]],122,129))/1e4 -
         as.numeric(substr(files[[2]],113,120))/1e4
m_dur <- pmin(T_dur,P_dur) 

# CHECK: m_dur = minimun duration
hist(m_dur, breaks=30, col="salmon", main="Common hourly prec and temp time series length [years]")
sum(m_dur>15) # 142 stations

files[[1]] <- files[[1]][m_dur>15]
files[[2]] <- files[[2]][m_dur>15]
temp_id <- as.integer(substr(files[[1]], 109,113))
prec_id <- as.integer(substr(files[[2]], 107,111))

# CHECK: station IDs
stopifnot(all(temp_id == prec_id))
headtail(data.frame(TEMP=substr(files[[1]], 106,131), PREC=substr(files[[2]], 104,129)), 3)


# 1.2. Metadata ----------------------------------------------------------------

data("metaIndex", package="rdwd")
meta <- metaIndex[metaIndex$res=="hourly" & 
                  metaIndex$var=="air_temperature" &
                  metaIndex$per=="historical" &
                  metaIndex$Stations_id %in% prec_id &
                  metaIndex$Stations_id %in% temp_id, ] 
meta$id <- meta$Stations_id
meta$name <- meta$Stationsname
meta$state <- meta$Bundesland
meta$lat <- meta$geoBreite
meta$long <- meta$geoLaenge
meta$ele <- meta$Stationshoehe

meta$T_dur <- T_dur[m_dur>15]
meta$P_dur <- P_dur[m_dur>15]
meta$m_dur <- m_dur[m_dur>15]
meta <- sortDF(meta, "id", decreasing = FALSE)
meta <- meta[,-(1:12)]
rownames(meta) <- NULL
save(meta, file="dataprods/meta.Rdata")

rm(T_dur, P_dur, m_dur, prec_id, temp_id, metaIndex)

# CHECK: coordinates #   View(meta)
colPoints(geoLaenge, geoBreite, m_dur, data=meta, add=F, asp=1.5)


# 1.3. Download data -----------------------------------------------------------

# Actually download files:
#fnamesT <- dataDWD(files[[1]], read=FALSE, sleep=30) # 35 min
#fnamesP <- dataDWD(files[[2]], read=FALSE, sleep=20) # 22 min
#save(fnamesT, file="dataprep/fnamesT.Rdata")
#save(fnamesP, file="dataprep/fnamesP.Rdata")
rm(files)

# read files:
load("dataprep/fnamesT.Rdata"); load("dataprep/fnamesP.Rdata")
library("rdwd")
#temp <- readDWD(fnamesT) # 10-27 min (because of unzipping etc)
#prec <- readDWD(fnamesP) #  4- 8 min (smaller files, shorter time series)
#save(temp, file="dataprep/temp.Rdata")
#save(prec, file="dataprep/prec.Rdata")
rm(fnamesP, fnamesT)

# After R restart, analysis can be continued here
load("dataprep/temp.Rdata"); load("dataprep/prec.Rdata")

# CHECK: number of non NA entries:
nnNAp <- sapply(prec, function(x) sum(!is.na(x$NIEDERSCHLAGSHOEHE)))
nnNAt <- sapply(temp, function(x) sum(!is.na(x$LUFTTEMPERATUR)))
hist(nnNAp/1000, breaks=20, col="blue")
hist(nnNAt/1000, breaks=20, col="red")
rm(nnNAp, nnNAt)

# Note: missing entries are replaced with NAs, so temperatures for 5 hour average don't get lost
# merge prec and temp (# 2 min):
PT_all <- pblapply(seq_along(temp), function(i) 
  {
  x <- prec[[i]]$MESS_DATUM
  hours <- seq.POSIXt(from=min(x,na.rm=T), to=max(x,na.rm=T), by="hour")
  x <- merge(prec[[i]], data.frame(MESS_DATUM=hours), all=TRUE)
  merge(x, temp[[i]], by="MESS_DATUM")
  }) 
rm(prec, temp)
names(PT_all) <- sapply(PT_all, "[", 1, 2)
save(PT_all, file="dataprep/PT_all.Rdata")


# 1.4. Dewpoint temperature ----------------------------------------------------

load("dataprep/PT_all.Rdata")
# CHECK: dew point temperature:
pdf("fig/Potsdam_Temp-Hum.pdf", height=5)
source("Code_aid.R")
x <- PT_all[[104]] # Potsdam 104, ID 3987
hist(x$REL_FEUCHTE, col="orange", breaks=40 )
library(hexbin) # gplot.hexbin
plot(hexbin(x$LUFTTEMPERATUR, x$REL_FEUCHTE, xbins=80), colramp=seqPal,
             xlab="hourly air temperature  [?C", ylab="Relative Humidity  [%]",
             main="Potsdam 1995-2015")
# attr(methods(class=class(hexbin(x$temp, x$hum))), "info")
x$dewtemp <- aid$dewtemp(x$LUFTTEMPERATUR, x$REL_FEUCHTE) 
p <- plot(hexbin(x$LUFTTEMPERATUR, x$dewtemp, xbins=80), colramp=seqPal) 
pushHexport(p$plot.vp) ; grid::grid.abline(0,1) ; rm(p)
# german wikipedia formula:
x$dewtemp2 <- aid$dewtemp2(x$LUFTTEMPERATUR, x$REL_FEUCHTE) 
hist(x$dewtemp-x$dewtemp2, breaks=40, col="purple") 
plot(x$LUFTTEMPERATUR, x$NIEDERSCHLAGSHOEHE, log="y", pch=16, col=addAlpha(1)) 
plot(x$dewtemp,        x$NIEDERSCHLAGSHOEHE, log="y", pch=16, col=addAlpha(1)) 
rm(x)
dev.off()


# 1.5. Dew point temperature of previous 5 hours -------------------------------

PT <- pblapply(PT_all, function(x) {                 # 10 seconds
  # Reduce filesize: column selection
  x <- x[, c("MESS_DATUM", "LUFTTEMPERATUR", "REL_FEUCHTE", "NIEDERSCHLAGSHOEHE")]
  colnames(x) <- c("date", "temp", "hum", "prec")
  # dew point temperature
  x$dewtemp <- aid$dewtemp(x$temp, x$hum) 
  # fill single missing NAs in dewtemp
  x$dtna <- approx(x$dewtemp, n=nrow(x))$y
  x$dtna[is.na(x$dewtemp[-1]) & is.na(x$dewtemp[-nrow(x)])] <- NA
  # dew point temperature of the preceding 5 hours (excluding the hour of rainfall): ## ToDo: note in paper
  x$temp5 <- filter(x$dtna, c(0,rep(1/5,5)), sides=1)
  x$dtna <- NULL
  ### x$temp51 <- pbsapply(1:nrow(x), function(i) mean(x[pmax(1,i-5):(i-1),"dewtemp"], na.rm=TRUE)) # ca 20-30 secs per run!
  ### x$temp52 <- filter(x$dewtemp, c(0,rep(1/5,5)), sides=1) # ca 0.01 sec
  # Reduce filesize: keep only rainfall values >0.5 mm/h:
  x <- x[!is.na(x$prec) & !is.na(x$temp5) & x$prec>=0.5, ]
  x})
names(PT) <- names(PT_all)
rm(PT_all)
save(PT, file="dataprods/PT.Rdata")



# 2. SSD: Sample size dependency -----------------------------------------------

load("dataprods/PT.Rdata")
# All 136k rainfall records between 10 and 12 degrees event dewpoint temperature # ToDo: note in paper
PREC <- unlist(lapply(PT, function(x) x[x$temp5>10 & x$temp5<12, "prec"] ))
PREC <- unname(PREC)
save(PREC, file="dataprods/PREC.Rdata")

hist(PREC, breaks=50, col="deepskyblue1")
logHist(PREC)
logHist(PREC, breaks=80)
binprec <- function(temp, ...)
  {
  bin <- temp+c(-1,1)
  PREC <- unlist(lapply(PT, function(x) x[x$temp5>bin[1] & x$temp5<bin[2], "prec"] ))
  logHist(PREC, breaks=80, main=paste0("Histogramm of all hourly rainfall ",
          "records at 142 stations in ", formatC(round(temp,1), format="f", digits=1),
          "\U{00B0}C bin"), xlab="Precipitation  [mm/h]", ...)
  title(sub=format(length(PREC), big.mark="'"), adj=1)
  d <- density(log10(PREC))
  lines(d$x, d$y*4)
  invisible(PREC)
  }

pdf("fig/binprec.pdf", height=5)
allprecs <- pblapply(seq(-11,23,0.2), binprec, xlim=log10(c(0.5,80)), 
                     ylim=lim0(6), col="deepskyblue1", freq=FALSE)
dev.off()

# 2.1. SSD computation ---------------------------------------------------------

source("Code_aid.R"); aid$load("PREC"); library(extremeStat)

ransample <- function(simn, trunc=0) # random sample generator
  {
  set.seed(simn) # reproducible 'random' numbers
  out <- lapply(aid$n, function(nn) sample(PREC,nn) )
  if(trunc>0) out <- lapply(out, function(x) x[x>=quantileMean(x,trunc)])
  out
  }
  
# custom weights come from section 2.3, determined from 400 runs (*_firstrun folders)
# truncation of 80% comes from section 2.4.
qn <- function(simn)
  {
  cat("\n------- sim run ", simn, " started ", as.character(Sys.time()), " -------\n", 
      file=paste0("simlogs/",simn,".txt"), append=TRUE)
  berryFunctions::tryStack({
  # Object and file name (with simulation run number):
  obname <- paste0("QN",simn)
  fname <- paste0("sim/QN",simn,".Rdata")
  if(file.exists(fname)) return()
  # Quantile estimation function:
  Qest <- function(x) 
    extremeStat::distLquantile(x, probs=aid$probs, truncate=0.8, 
    quiet=TRUE, order=FALSE, weightc=NA)#dweight)
  # random samples
  rs <- ransample(simn)
  # Hardcore computation returning a 3D array:
  QN <- vapply(rs, Qest, FUN.VALUE=array(0, dim=c(38,5)) )
  # Dimnames
  dimnames(QN)[3] <- list(paste(aid$n))
  names(dimnames(QN)) <- c("distr","prob", "n")
  # Saving to disc:
  assign(obname, QN, envir=environment())
  save(list=obname, file=fname)
  }, if(simn>15) warn=0, file=paste0("simlogs/",simn,".txt")  )# tryStack end
  cat("\n------- sim run ", simn, " finish  ", as.character(Sys.time()), " -------\n", 
      file=paste0("simlogs/",simn,".txt"), append=TRUE)
  }

# long computing time (2-2.5 minutes per simulation run):
# 400 in 2 hours on 7 cores
library(parallel) # for parallel lapply execution
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("aid", "PREC", "qn", "ransample"))#, "dweight"))
dummy <- pblapply(X=1:500, cl=cl, FUN=qn)
stopCluster(cl)
rm(cl, dummy)
rm(qn, ransample)



# 2.2. SSD aggregation -------------------------------------------------------

# Read in simulation results (1.4 GB for 2000 simulations!)
simEnv <- new.env()
dummy <- pblapply(dir("sim", full=TRUE), load, envir=simEnv) # 6 secs / 1 min for 500 sims
simQ <- as.list(mget(gtools::mixedsort(ls(envir=simEnv)), envir=simEnv))
simQ <- l2array(simQ) # this takes a minute, 400MB for 745 sims
save(simQ, file="dataprods/simQ.Rdata")
rm(simEnv, dummy)

# aggregate (2 minutes for 500 runs):
load("dataprods/simQ.Rdata")
simQA <- pbapply(simQ, MARGIN=1:3, quantileMean, 
                 probs=c(seq(0,1,0.1),0.25,0.75), na.rm=TRUE) 
                #probs=c(0.3,0.5,0.7), na.rm=TRUE)
save(simQA, file="dataprods/simQA.Rdata")
load("dataprods/simQA.Rdata")


# checks:
str(simQ)
sort(table(rownames(which(simQ[1:35, "99.9%", ,]>500, arr.ind=TRUE))))


dim(simQA)
str(simQA)
dimnames(simQA)
simQA["0%", ,"99.9%", as.character(40:45)]

which(simQ<0, arr.ind=TRUE) # some laplace
simQ[which(simQ<0)] # nothing really bad
which(simQ[,1:4,,]<0.5, arr.ind=TRUE)
which(simQ["kap", "99.9%", ,]>200, arr.ind=TRUE) # 2
str(which(simQ[1:35, , ,]>900, arr.ind=TRUE)) # 12k (45k>200) at 400 sims
toolarge <- which(simQ[1:35,"99%", ,]>500, arr.ind=TRUE)
head(toolarge)
table(rownames(toolarge))

kap <- as.vector(simQ["kap","99.9%",,])
val <- logSpaced(min=120, max=250, n=100, plot=F)
numlarger <- sapply(val, function(x) sum(kap>x, na.rm=T))
plot(val, numlarger, log="y", type="o", las=1)

val <- logSpaced(min=50, max=10000, n=30, plot=F)
numlarger <- pbsapply(val, function(x) sum(simQ[20:35, "99.9%",,]>x, na.rm=T))
plot(val, numlarger, log="yx", type="o", axes=F, main="Q99.9% GPD estimates larger than val")
logAxis(1); logAxis(2)



# 2.3. Distribution weights ----------------------------------------------------

load("dataprods/simQA.Rdata"); load("dataprods/PREC.Rdata"); source("Code_aid.R")

# _a. bias -----
# 136k Rainfall hours between 10 and 12 degrees for all stations (PREC)
Qtrue <- quantileMean(PREC, 0.999)
drmse <- sapply(dimnames(simQA)[[2]][1:35], function(d) rmse(simQA["50%",d,"99.9%",], rep(Qtrue,512)))
w_bias <- 3-drmse[1:17]
w_bias[w_bias<0] <- 0
w_bias <- w_bias/sum(w_bias)

pdf("fig/simQn.pdf", height=5)
par(mar=c(3.5,3.5,2,0.5), mgp=c(2.1,0.7,0), las=1)
dn <- c(names(sort(drmse)),"")#[1:2]
dcol <- rep("grey80", length(dn))
names(dcol) <- dn
dcol[grepl("GPD_", dn)] <- addAlpha("red")
dummy <- pblapply(dn, function(d)
  {
  plot(1, type="n", xlim=c(24,1900), ylim=c(5,20), log="x", xaxs="i", main=d, xaxt="n",
       xlab="sample size", ylab="Random sample 99.9% quantile estimate")
  logAxis(1)
  for(dd in dn[dn!=""]) lines(aid$n, simQA["50%",dd,"99.9%",], col=dcol[dd])
  abline(h=quantileMean(PREC, 0.999), lty=3)
  if(dn=="") return()
  ciBand(yl=simQA["30%",d,"99.9%",], 
         ym=simQA["50%",d,"99.9%",],
         yu=simQA["70%",d,"99.9%",], x=aid$n, colm="blue", add=TRUE)
  title(main=round(drmse[d],2), adj=0)
})
rm(dummy)
dev.off()

# ToDo: paper Fig 5 + 6_potsdamQn
# ToDo: recreate paper Fig 6 SSD GPD 



# _b. goodness of fit -----

library(extremeStat)
dlf <- distLfit(PREC, truncate=0.8, weightc=dweight)
dlf$gof

load("dataprods/PT.Rdata")
PRECstats <- lapply(PT, function(x) x[x$temp5>10 & x$temp5<12, "prec"] )
dlfstats <- pblapply(PRECstats, distLfit, truncate=0.8, progbars=FALSE, 
                     time=FALSE, plot=FALSE) # 2 mins with plot, 1 min without

dlfgofs <- lapply(dlfstats, "[[", "gof")
dlfgofs <- do.call(rbind, dlfgofs)
dlfgofs$dist <- sapply(strsplit(rownames(dlfgofs), ".", fix=T), "[", 2)
dn <- unique(dlfgofs$dist)
dn <- sapply(dn, function(d) mean(dlfgofs[dlfgofs$dist==d, "RMSE"]) )
dn <- names(sort(dn))

pdf("fig/distribution_gofs.pdf", height=5)
for(d in dn){
lh <- logHist(dlfgofs$RMSE, breaks=60, las=1, col=addAlpha("darkorange"), border=NA, main=d)
logHist(dlfgofs[dlfgofs$dist==d, "RMSE"], breaks=lh$breaks, col=1, logargs=list(xaxt="n"), add=TRUE)
}
rm(d)
dev.off()


# _c. error rate -----

load("dataprods/simQ.Rdata"); source("Code_aid.R")
simQ[simQ<0.5] <- NA
simQ[1:35,,,][simQ[1:35,,,]>200] <- NA

simNA <- apply(simQ, 1:3, function(x) mean(is.na(x)))
w_error <- apply(simQ, 1, function(x) mean(is.na(x)))
dn <- names(sort(w_error))


NAs <- apply(simQ[1:35, "99%", ,], MARGIN=1, function(x) mean(is.na(x))*100)
data.frame(NA_percent=sort(NAs))


pdf("fig/distribution_errorrates.pdf", height=5)
dummy <- pblapply(dn, function(d)
{
plot(1, xlim=c(25,2000), ylim=lim0(0.7), type="n", log="x", xaxt="n", 
     xlab="Sample size", ylab="Proportion of NAs across simulations", las=1)
logAxis(1)
title(main=d)
for(dd in dn) for(p in names(aid$probcols)) lines(aid$n, simNA[dd,p,], col=8)
for(p in names(aid$probcols)) lines(aid$n, simNA[d,p,], col=aid$probcols[p])
})
rm(dummy)
dev.off()


# _d. Visualization -----
save(dweight, file="dataprods/dweight.Rdata")

pdf("fig/distribution_weights.pdf", height=5)
distLgofPlot(dlf, ranks=FALSE, lwd=1.5)
barplot(sort(drmse), horiz=TRUE, las=1, main="Deviation from true quantile (RMSE along sample sizes)")
barplot(sort(w_bias), horiz=TRUE, las=1, main="Relative weights")
dev.off()




# 2.4. Truncation dependency ---------------------------------------------------
# ToDo: recreate paper Fig 3 + 4 (new)


# 2.5. tempdep Wakeby distribution ---------------------------------------------
# ToDo: recreate paper Fig 7



# 2.6. Sample size bias lin vs log ---------------------------------------------

load("dataprods/PREC.Rdata"); source("Code_aid.R")
# lower n-dependency if log10(randomsample) is used for fitting:
log_n <- c(25:100,seq(150,500, by=50))
log_qn <- function(...) berryFunctions::tryStack( vapply(log_n, function(nn){
    ransam <- sample(PREC,nn)
    lin <- extremeStat::distLquantile(ransam, sel=c("wak","gpa"), probs=aid$probs, 
    truncate=0.8, addinfo=FALSE, weighted=FALSE, gpd=FALSE, order=FALSE, quiet=TRUE)
    log <- extremeStat::distLquantile(log10(ransam), sel=c("wak","gpa"), probs=aid$probs, 
    truncate=0.8, addinfo=FALSE, weighted=FALSE, gpd=FALSE, order=FALSE, quiet=TRUE)
    rownames(log) <- paste0("log_",rownames(log))
    rbind(lin,log) }, FUN.VALUE=array(0, dim=c(6,4)) ), silent=TRUE)

# 80? runs per minute on 4 cores:
library(parallel)
cl <- makeCluster( detectCores()-0 )
clusterExport(cl, c("PREC", "log_qn", "aid", "log_n"))
log_QN <- pblapply(X=1:800, cl=cl, FUN=log_qn)
save(log_n, log_QN, file="dataprods/log_QN.Rdata")
stopCluster(cl) ; rm(cl)


load("dataprods/log_QN.Rdata") ; load("dataprods/PREC.Rdata")
log_QNA <- l2array( log_QN[sapply(log_QN, class)=="array"] )
log_QNAg <- pbapply(log_QNA, MARGIN=1:3, quantileMean, probs=c(seq(0,1,0.1),0.25,0.75), na.rm=TRUE)

log_plot <- function(d, fun=I, ...) 
  {
  ciBand(yl=fun(log_QNAg["30%",d,"99.9%",]), xlab="Sample size", ylim=c(0.7,1.8),
         ym=fun(log_QNAg["50%",d,"99.9%",]), yaxt="n", 
         yu=fun(log_QNAg["70%",d,"99.9%",]), x=log_n, ...)
  logAxis(2)
  abline(h=c(quantileMean(PREC,0.999), quantileMean(log10(PREC),0.999)), lty=3)
  legend("topright", d, bty="n", inset=0.04)
  }

pdf("fig/loglin_samplesizedependency.pdf", height=5)
par(mfrow=c(1,2), mar=c(3,3,2,0.4), mgp=c(1.9,0.6,0))
log_plot("wak", fun=log10,          ylab="log( distLquantile(sample) )", main="linear")
log_plot("log_wak",                 ylab="distLquantile( log(sample) )", main="logarithmic")
log_plot("gpa", fun=log10,          ylab="log( distLquantile(sample) )", main="linear")
log_plot("log_gpa",                 ylab="distLquantile( log(sample) )", main="logarithmic")
log_plot("quantileMean", fun=log10, ylab="log( distLquantile(sample) )", main="linear")
log_plot("log_quantileMean",        ylab="distLquantile( log(sample) )", main="logarithmic")
dev.off()



# 3. Hourly Prec-Temp relationship ---------------------------------------------

# 3.1. Raw data visualisation --------------------------------------------------

load("dataprods/meta.Rdata"); load("dataprods/PT.Rdata"); source("Code_aid.R")
range(sapply(PT, function(x)   max(x$prec, na.rm=T))) # 21.5 80.8
range(sapply(PT, function(x) range(x$temp, na.rm=T))) # -16.4  34.6
range(sapply(PT, function(x) range(x$temp5,na.rm=T))) # -21.3  29.8
hist(sapply(PT, nrow), breaks=30, col="salmon")

library("mapdata") ;  map <- maps::map('worldHires','Germany') ;  dev.off()

PT5 <- pblapply(1:142, function(i){
  x <- PT[[i]]
  data.frame(x=x$temp5[x$prec>2], y=x$prec[x$prec>2])
  })

# Plot PT-graph with 5 hour preceding dewpoint temp: 
pdf("fig/RawData_T5.pdf", height=5)
dummy <- pblapply(1:142, function(i){  # order(meta$ele), ...
  aid$stationplot(i, meta, map, ylim=c(2.5,70) )
  points(PT5[[i]], pch=16, col=addAlpha(1))
  text(-3,50, 50)
  })
dev.off()

PT5df <- do.call(rbind, PT5)
PT5df$y <- log10(PT5df$y)

pdf("fig/RawData_T5_all.pdf", height=5)
dummy <- pblapply(20:200, function(n){
gplots::hist2d(PT5df, col=seqPal(n), nbins=n, yaxt="n", xlim=c(-5.1,24),
xlab="Dew point temperature (mean of preceding 5 hours)  [ \U{00B0}C]",
ylab="Precipitation  [mm/h]", main=n)#"All stations")
logAxis(2)
})
dev.off()


# toDo: recreate paper Fig 2 (PT emp) or use fig/PTQstats p104 with different ylim

rm(dummy, map)


# 3.2. PT-quantiles computation ------------------------------------------------

load("dataprods/PT.Rdata"); source("Code_aid.R") # mid, probs
load("dataprods/cweights.Rdata")

# Change aid$mid to start at 4.8 for plotting, change 201 to 203

# long computing time (1 minute per station)
library(parallel) # for parallel lapply execution
cl <- makeCluster( detectCores()-0 )
clusterExport(cl, c("PT","aid", "cweights"))
PTQ <- pblapply(X=1:142, cl=cl, FUN=function(i)
  {
  x <- PT[[i]]   #  x <- PT[[3]]; t=19.5
  # Quantile estimates per temperature bin
  binQ <- lapply(aid$mid, function(t) {
    seldat <- x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)]
    extremeStat::distLquantile(seldat[!is.na(seldat)], probs=aid$probs, truncate=0.8, addinfo=TRUE,
                  weightc=NA, order=FALSE, ssquiet=TRUE, time=FALSE, progbars=FALSE)
    })
  # Transform into array for faster subsetting:
  binQ2 <- array(unlist(binQ), dim=c(38, 4, 201),  # c(nrow(binQ[[1]]), ncol(binQ[[1]]), length(binQ))
       dimnames=list(distr=rownames(binQ[[1]]), prob=colnames(binQ[[1]]), temp=aid$mid))
  return(binQ2)
  })
save(PTQ, file="dataprods/PTQ.Rdata")     # 30 mins
stopCluster(cl)
rm(cl)

length(PTQ) # 142 stations
str(PTQ[[1]]) # each at 201 temperature bins
aid$mid[100] # for 14.9 °C:
PTQ[[1]][,,100]

head(PT[[113]])
n113 <- sapply(aid$mid, function(t) sum(PT[[113]]$temp5>(t-1) & PT[[113]]$temp5<=(t+1)) )
all(n113 == PTQ[[113]]["n_full",1,] )
rm(n113)


# 3.3. PT-quantiles visualization ----------------------------------------------
load("dataprods/PTQ.Rdata"); source("Code_aid.R")

PTQlines <- function(
prob="",  
dn="",
cut=150,
col=addAlpha("black"),
...
)
{
stats <- sapply(PTQ, function(x){
         lines(aid$mid, x[dn, prob, ], col=col, ...)
         x[dn, prob, ]})
# average of stations (bins with >50 values, rainQuantile>150 ignored):
stats <- replace(stats, stats>cut, NA)       ### toDo: make sure this is mentioned in paper!
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
statav
}

pdf("fig/PTQ.pdf", height=5)
# _a. empirical and parametric quantiles in one plot ----------------------------
par(mar=c(4,4,2,0.5), mgp=c(2.5,0.6,0))
dn <- c("quantileMean", "weighted2")              ### later: use weightedc!
dc <- addAlpha(c("red", "blue"))
for(prob in c("90%", "99%", "99.9%", "99.99%"))
{
aid$PTplot(prob=prob)
dummy <- sapply(PTQ, function(x){
         for(d in 1:2) lines(aid$mid, x[dn[d], prob, ], col=dc[d]) })
}
# _b. Qemp/par side by side -----------------------------------------------------
par(mfrow=c(1,2), mar=c(2,2,1.5,0.5), oma=c(1.5,1.5,1.5,0) )
for(prob in c("90%", "99%", "99.9%", "99.99%"))
{
for(d in 1:2)
{
aid$PTplot(prob=prob, outer=TRUE, line=0)
title(main=dn[d])
statav <- PTQlines(prob=prob, dn=dn[d], col=dc[d])
lines(aid$mid, statav, lwd=2)  
}}
#
# _c. PTQ for each distribution function ----------------------------------------
par(mfrow=c(1,1), mar=c(4,4,2,0.5), oma=c(0,0,0,0) )
dummy <- pblapply(dimnames(PTQ[[1]])$distr, function(dn){ 
ylim <- c(5,130)
cut <- 150
if(dn=="n_full")    {ylim <- c(5,2000); cut <- 1e5}
if(dn=="n")         {ylim <- c(1, 400); cut <- 1e5}
if(dn=="threshold") {ylim <- c(0.5,25); cut <- 1e5}
aid$PTplot(prob="99.9%", ylim=ylim, main=dn, cc=!dn %in% c("n_full","n","threshold"))
statav <- PTQlines(prob="99.9%", dn=dn, cut=cut)
lines(aid$mid, statav, lwd=2, col="red")
if(dn=="n_full") abline(h=25, lwd=2, col="red")
if(dn=="n") abline(h=5, lwd=2, col="red")
})
rm(dummy, prob, d, dn, dc, statav)
dev.off()

# ToDo: recreate paper Fig 8


# 3.4. PTQ per station ---------------------------------------------------------

load("dataprods/PTQ.Rdata"); load("dataprods/PT.Rdata"); load("dataprods/meta.Rdata"); source("Code_aid.R")
library("mapdata") ;  map <- maps::map('worldHires','Germany') ;  dev.off()


dn <- c("quantileMean","weighted2","gpa", "wak")
dc <- c("brown1", "deepskyblue4", "darkolivegreen4", "peru")
pdf("fig/PTQ_stats.pdf", height=5)
dummy <- pblapply(1:142, function(i){
  aid$stationplot(i, meta, map, xlim=c(4.8,21), ylim=c(5,130) )
  legend("topleft", rev(dn), col=rev(dc), lwd=2, inset=c(0.064,0), bg="white")
  points(PT[[i]][PT[[i]]$prec>5,c("temp5","prec")], pch=16, col=addAlpha(1))
  for(d in 1:4) lines(aid$mid, PTQ[[i]][dn[d], "99.9%", ], col=dc[d], lwd=2)
  })
rm(dn, dc, dummy)
dev.off()

