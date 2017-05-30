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
# In Rstudio, turn on the outline menu (topright button in script window, CTRL + SHIFT + O)

# 0. Packages

# 1. Data 
# 1.1. Select DWD stations
# 1.2. Meta data weather stations
# 1.3. Download DWD station data
# 1.4. Dewpoint temperature cumputation
# 1.5. Dew point temperature of previous 5 hours
# 1.6. Raw data visualisation

# 2. Sample size dependency
# 2.1. SSD computation
# 2.2. SSD checks
# 2.3. SSD visualisation


# 3. Hourly Prec-Temp relationship
# 3.1. Raw data visualisation
# 3.2. PT-quantiles computation
# 3.3. PT-quantiles visualization
# 3.4. PTQ per station


# 4. tempdep distribution



# 0. Packages ------------------------------------------------------------------
if(FALSE){ # You need to download and install the packages only once
packs <- c("berryFunctions", "extremeStat", "pbapply", "maps", "gplots", "gtools",
         "mapdata", "OSMscale", "RCurl", "rdwd", "hexbin")
isinst <- sapply(packs, requireNamespace, quietly=TRUE)
install.packages(packs[!isinst])
rm(packs, isinst)
# If OSMscale does not install automatically, check out:
browseURL("https://github.com/brry/OSMscale#intro")
}
library(berryFunctions); library(pbapply) # potentially needed in every section



# ......... ----
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
source("Code_aid.R")

# CHECK: dew point temperature:
pdf("fig/Potsdam_Temp-Hum.pdf", height=5)
x <- PT_all[[104]] # Potsdam 104, ID 3987
hist(x$REL_FEUCHTE, col="orange", breaks=40 )
library(hexbin) # gplot.hexbin
plot(hexbin(x$LUFTTEMPERATUR, x$REL_FEUCHTE, xbins=80), colramp=seqPal,
             xlab="hourly air temperature  [\U{00B0}C]", ylab="Relative Humidity  [%]",
             main="Potsdam 1995-2015")
# attr(methods(class=class(hexbin(x$temp, x$hum))), "info")
x$dewtemp <- aid$dewtemp(x$LUFTTEMPERATUR, x$REL_FEUCHTE) 
p <- plot(hexbin(x$LUFTTEMPERATUR, x$dewtemp, xbins=80), colramp=seqPal) 
pushHexport(p$plot.vp) ; grid::grid.abline(0,1) ; rm(p)
# german wikipedia formula:
x$dewtemp2 <- aid$dewtemp2(x$LUFTTEMPERATUR, x$REL_FEUCHTE) 
hist(x$dewtemp-x$dewtemp2, breaks=40, col="purple") 
plot(x$LUFTTEMPERATUR, x$NIEDERSCHLAGSHOEHE, log="y", pch=16, xlim=c(-20,30), col=addAlpha(1)) 
plot(x$dewtemp,        x$NIEDERSCHLAGSHOEHE, log="y", pch=16, xlim=c(-20,30), col=addAlpha(1)) 
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


# 1.6. Raw data visualisation --------------------------------------------------

source("Code_aid.R"); aid$load("meta", "PT")
library("mapdata") ;  map <- maps::map('worldHires','Germany') ;  dev.off()

range(sapply(PT, function(x)   max(x$prec, na.rm=T))) # 21.5 80.8
range(sapply(PT, function(x) range(x$temp, na.rm=T))) # -16.4  34.6
range(sapply(PT, function(x) range(x$temp5,na.rm=T))) # -21.3  29.8
hist(sapply(PT, nrow), breaks=30, col="salmon")

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

rm(dummy, map, PT5, PT5df)



# ......... ----
# 2. SSD: Sample size dependency -----------------------------------------------
# 2.1. large dataset ----
load("dataprods/PT.Rdata")
PREC <- unname(unlist(lapply(PT, "[", "prec")))
logHist(PREC, breaks=50)
save(PREC, file="dataprods/PREC.Rdata")

# Distribution fit quality 
source("Code_aid.R"); aid$load("PREC"); library(extremeStat)

dlf <- distLfit(log10(PREC), truncate=0.9) # hangs at 82%+88% CDFS -> 3.3 min
dlq <- distLquantile(dlf=dlf, probs=aid$probs, truncate=0.9, gpd=FALSE, list=TRUE) 
dlq <- distLquantile(log10(PREC), truncate=0.9, sel="gpa", gpd=F,weight=F, prob=0.999)




pdf("fig/dists_fullsample.pdf")
p <- c(0.999,0.9999)
plotLquantile(distLquantile(dlf=dlf,probs=p,truncate=0.9,gpd=F,list=T), 
              log=TRUE, legargs=list(bg="white"), nbest=5, breaks=50)
points(quantile(log10(PREC),p), c(0.3,0.3), pch=16)
text(quantile(log10(PREC), 0.999), 0.6, "99.9%", adj=0.8)
text(quantile(log10(PREC), 0.9999), 0.6, "99.99%")
dev.off()

dlf <- distLfit(log10(PREC), truncate=0.95) # 1 min
plotLfit(dlf, log=T)
plotLquantile(distLquantile(dlf=dlf,probs=p,truncate=0.95,gpd=F,list=T), 
              log=TRUE, legargs=list(bg="white"), nbest=5, breaks=50)
points(quantile(log10(PREC),p), c(0.3,0.3), pch=16)

rm(p,dlq,dlf)


# 2.2. SSD computation ---------------------------------------------------------

source("Code_aid.R"); aid$load("PREC"); library(extremeStat)
if(packageVersion("extremeStat")<"1.3.2") stop("extremeStat too old. Please install from github.com/brry/extremeStat")

ransample <- function(simn) # random sample generator
  {
  set.seed(simn) # reproducible 'random' numbers
  out <- lapply(aid$n, function(nn) sample(log10(PREC),nn) )
  out
  }
  
# truncation of 90% decided in supplement
qn <- function(simn)
  berryFunctions::tryStack({
  # Object and file name (with simulation run number):
  obname <- paste0("QN",simn)
  fname <- paste0("sim_ssd/QN",simn,".Rdata")
  if(file.exists(fname)) return()
  # Quantile estimation function:
  Qest <- function(x) 
    extremeStat::distLquantile(x, probs=aid$probs, truncate=0.9, quiet=TRUE, 
                               weighted=FALSE, sel="gpa", order=FALSE)
  # random samples
  rs <- ransample(simn)
  # Hardcore computation returning a 3D array:
  QN <- vapply(rs, Qest, FUN.VALUE=array(0, dim=c(19,5)) )
  # Dimnames
  dimnames(QN)[3] <- list(paste(aid$n))
  names(dimnames(QN)) <- c("distr","prob", "n")
  # Saving to disc:
  assign(x=obname, value=QN, envir=environment())
  save(list=obname, file=fname)
  }, file=paste0("sim_log/",simn,".txt")  )

if(!file.exists("sim_ssd")) dir.create("sim_ssd")
if(!file.exists("sim_log")) dir.create("sim_log")
# long computing time (1.5 minutes per simulation run): 500 in 1:36 hours on 7 cores (4h on 3)
library(parallel) # for parallel lapply execution
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("aid", "PREC", "qn", "ransample"))
dummy <- pblapply(X=501:1000, cl=cl, FUN=qn)
stopCluster(cl)
rm(cl, dummy)
rm(qn, ransample)



# 2.3. SSD array + checks ------------------------------------------------------

# Read in simulation results + convert to array
# (1.4 GB for 2000 simulations!)
simEnv <- new.env()
dummy <- pblapply(dir("sim_ssd", full.names=TRUE)[1:10], load, envir=simEnv)
dummy <- pblapply(dir("sim_ssd", full.names=TRUE), load, envir=simEnv) # 6 secs / 1 min for 500 sims
simQL <- as.list(mget(gtools::mixedsort(ls(envir=simEnv)), envir=simEnv))
simQ <- l2array(simQL)
names(dimnames(simQ))[4] <- "simnr"
rm(simEnv, dummy, simQL)

# SSD simulation aggregation:
simQA <- pbapply(simQ, MARGIN=1:3, quantileMean, probs=c(0.3,0.5,0.7), na.rm=TRUE) # 1 min
save(simQA, file="dataprods/simQA.Rdata")

# checks:
str(simQ)
dim(simQ)
lapply(dimnames(simQ), head, 50)
apply(simQ[1:15, "99.9%", as.character(40:45),], 1:2, min, na.rm=TRUE)
min(simQ[1:15,,,], na.rm=T)

# Estimates that are obviously too large
val <- logSpaced(min=50, max=10000, n=30, plot=F)
numlarger <- pbsapply(val, function(x) sum(simQ[-(2:3), "99.9%",,]>log10(x), na.rm=T))
plot(val, numlarger, log="yx", type="o", axes=F, main="Q99.9% GPD estimates larger than val")
logAxis(1); logAxis(2)

toolarge99 <- which(simQ[-(2:3),"99%",,]>log10(500), arr.ind=TRUE)
head(toolarge99)
sort(table(rownames(toolarge99)))
rm(toolarge99, numlarger, val)



# 2.4. SSD visualisation -----
source("Code_aid.R"); aid$load("simQA", "PREC")

pdf("fig/fig2.pdf", height=3, width=3.5, pointsize=10)
par(mar=c(3,3,0.2,0.2), mgp=c(1.8,0.7,0), las=1, lend=1)
plot(1, type="n", xlim=c(50,730), ylim=log10(c(5,30)), xaxs="i", main="", yaxt="n",
       xlab="sample size n", ylab="")
logAxis(2)
for(d in c("empirical","gpa"))
  ciBand(yl=simQA["30%",d,"99.9%",], 
         ym=simQA["50%",d,"99.9%",],
         yu=simQA["70%",d,"99.9%",], x=aid$n, colm=if(d=="gpa") "blue" else "green3", add=TRUE)
abline(h=quantileMean(log10(PREC), probs=0.999), lty=3)
#abline(h=c(1.281135,1.25894), lty=3)
legend("bottomright", c("Parametric quantile", "Empirical quantile", 
                        "Central 40% of simulations", "Quantile of full sample"),
       lwd=c(2,2,11,1), lty=c(1,1,1,3), col=c("blue","green3",8,1), bg="white", cex=0.8)
text(80, log10(c(8.4, 22)), c("empirical","parametric"), col=c("green3","blue"), adj=0)
title(ylab="Random sample 99.9% quantile  [mm/h]       ")
dev.off()



# ......... ----
# 3. Hourly Prec-Temp relationship ---------------------------------------------

# 3.1. PT-quantiles computation ------------------------------------------------
source("Code_aid.R"); aid$load("PT"); library(extremeStat)

# long computing time (30 seconds per station, 10 min on 7 cores)
library(parallel) # for parallel lapply execution
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("PT","aid"))
PTQL <- pblapply(X=1:142, cl=cl, FUN=function(i){
  x <- PT[[i]]   #  x <- PT[[3]]; t=19.5
  # Quantile estimates per temperature bin
  binQ <- lapply(aid$mid, function(t) {
    seldat <- x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)]
    extremeStat::distLquantile(log10(seldat), probs=aid$probs, truncate=0.9,  
                               sel="gpa", order=FALSE, quiet=TRUE, weighted=FALSE)
    })
  binQ2 <- berryFunctions::l2array(binQ)
  names(dimnames(binQ2)) <- c("distr","prob","temp") ; dimnames(binQ2)[[3]] <- aid$mid
  return(binQ2)
  })
stopCluster(cl); rm(cl)
PTQ <- l2array(PTQL)
names(dimnames(PTQ))[4] <- "stat"; dimnames(PTQ)[[4]] <- names(PT)
save(PTQ, file="dataprods/PTQ.Rdata") 
rm(PTQL)

dim(PTQ) # 142 stations
str(PTQ) # each at 203 temperature bins
aid$mid[100] # for 14.7 °C:
PTQ[,,100,1]

head(PT[[113]])
n113 <- sapply(aid$mid, function(t) sum(PT[[113]]$temp5>(t-1) & PT[[113]]$temp5<=(t+1)) )
all(n113 == PTQ["n_full",1,,113], na.rm=TRUE)
rm(n113)


# 3.2. PT-quantiles visualisation ----------------------------------------------
source("Code_aid.R"); aid$load("PTQ")

PTQlines <- function(
prob="",  
dn="",
cut=150,
col=addAlpha("black"),
...
)
{
if(dn%in%c("n_full","n","threshold")) 
  {
  prob <- dimnames(PTQ)[[2]][1]
  for(i in 1:142) lines(aid$mid,    PTQ[dn,prob,,i], col=col, ...) 
  } else
  for(i in 1:142) lines(aid$mid, 10^PTQ[dn,prob,,i],  col=col, ...) 
stats <- PTQ[dn,prob,,]
# average of stations (bins with <50 values, log10(rainQuantile)>150 ignored):
stats <- replace(stats, stats>cut, NA)  
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
statav
}


pdf("fig/fig3.pdf", height=5)
par(mfrow=c(1,2), mar=c(2,2,0.5,0.5), oma=c(1.5,1.5,0,0) )
for(alpha in 0.15)# 1:6/20)
{
aid$PTplot(prob="99.9%", outer=TRUE, line=0, xlim=c(4.8,20.3), ylim=c(4,120), 
           main="", cc=FALSE)
statav_e <- PTQlines(prob="99.9%", dn="empirical", col=addAlpha("green3", alpha))
lines(aid$mid, 10^statav_e, lwd=3)  
aid$cc_lines(NA, mainargs=list(col=2, lty=2))
legend("topleft", "Empirical", bty="n")
#
aid$PTplot(prob="99.9%", outer=TRUE, line=5, xlim=c(4.8,20.3), ylim=c(4,120), 
           main="", cc=FALSE)
statav <- PTQlines(prob="99.9%", dn="gpa", col=addAlpha("blue", alpha))
lines(aid$mid, 10^statav_e, col="green3", lwd=3) 
lines(aid$mid, 10^statav, lwd=3) 
legend("topleft", "Parametric", bty="n")
aid$cc_lines(NA, mainargs=list(col=2, lty=2))
#title(main=alpha, line=-3, outer=T)
}
rm(alpha)
dev.off()


# 3.3. PTQ per station ---------------------------------------------------------

source("Code_aid.R"); aid$load("PTQ", "PT", "meta")
library("mapdata") ;  map <- maps::map('worldHires','Germany') ;  dev.off()

dn <- c("empirical","gpa")
dc <- c("green3", "blue")
pdf("fig/PTQ_stats.pdf", height=5)
dummy <- pblapply(order(meta$ele), function(i){
  aid$stationplot(i, meta, map, xlim=c(4.8,21), ylim=c(5,130) )
  legend("topleft", rev(dn), col=rev(dc), lwd=2, inset=c(0.064,0), bg="white")
  points(PT[[i]][PT[[i]]$prec>4,c("temp5","prec")], pch=16, col=addAlpha(1))
  for(d in 1:2) lines(aid$mid, 10^PTQ[dn[d], "99.9%", , i], col=dc[d], lwd=2)
  })
rm(d, dn, dc, dummy)
dev.off()




# 3.4. outlier station ----

plot(10^PTQ["gpa","99.9%","11",])
which.max(PTQ["gpa","99.9%","11",]) # IDberry 37, IDdwd 1346 Feldberg/Schwarzwald




# ......... ----
# 4. tempdep distribution ----------------------------------------------------
# 4.1. Parameter ----
source("Code_aid.R"); aid$load("PT"); library(lmomco); library(parallel) 

cl <- makeCluster( detectCores()-1 ) ; clusterExport(cl, "PT")

tdtemp <- seq(4,16,by=2) # 25 secs
tdpar <- pblapply(tdtemp, cl=cl, FUN=function(t) sapply(PT, function(x) # 7 x 20 secs
  {
  x <- x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)]
  lmomco::pargpa(lmom=lmomco::lmoms(log10(x)))$para
  }))
stopCluster(cl); rm(cl)


tdpar <- l2array(tdpar)
names(dimnames(tdpar)) <- c("par","stat","temp") ; dimnames(tdpar)[[3]] <- tdtemp
str(tdpar)
save(tdpar,tdtemp, file="dataprods/tdpar.Rdata") 


# 4.2. Computation ----
source("Code_aid.R"); aid$load("PT","tdpar"); library(lmomco); library(parallel)
xys <- lapply(1:3, function(param)
  {
  xy <- function(t,p) cbind(x=rep(tdtemp[t],142), y=tdpar[p,,t])
  xy <- do.call(rbind, lapply(seq_along(tdtemp), xy, p=param))
  xy <- as.data.frame(xy, check.names=FALSE)
  xy$station <- rownames(xy)
  rownames(xy) <- NULL
  xy
  })

# temp-dep simulation: temperature, sample size n, parameters, random numbers, quantiles
tdsimt <- 5:20
tdsimn <- pbsapply(tdsimt, function(t) sapply(PT, function(x) 
                                           sum(x$temp5>(t-1) & x$temp5<=(t+1))))
tdsimn <- colMeans(tdsimn)
tdsimp <- function(t) list(type="gpa", para=c(
   zeta=unname(predict(lm(y~x, data=xys[[1]]),newdata=data.frame(x=t))),
   beta=unname(predict(lm(y~x, data=xys[[2]]),newdata=data.frame(x=t))),
  delta=unname(predict(lm(y~x, data=xys[[3]]),newdata=data.frame(x=t)))
  ), source="pargpa")
tdsimr <- function() lapply(seq_along(tdsimt), function(i) lmomco::rlmomco(
  n=tdsimn[i], para=tdsimp(tdsimt[i]) ))

library(extremeStat)
tdsimq <- function(seed) 
  {
  set.seed(seed)
  out <- lapply(tdsimr(), extremeStat::distLquantile, truncate=0.9, gpd=FALSE, 
                order=FALSE, sel="gpa", weighted=FALSE, probs=aid$probs, quiet=TRUE)
  out <- berryFunctions::l2array(out)
  names(dimnames(out)) <- c("distr","prob","temp") ; dimnames(out)[[3]] <- tdsimt
  out
  }

cl <- makeCluster( detectCores()-1 ) # 1000 runs 2 min on 7 cores
clusterExport(cl, c(paste0("tdsim",c("t","n","p","r","q")), "xys","aid")) 
tdsimL <- pblapply(1:1000, cl=cl, FUN=tdsimq)
stopCluster(cl) ; rm(cl)

tdsim <- l2array(tdsimL)
save(tdsimt,tdsimn,tdsimp,tdsimr,tdsimq,tdsim,xys, file="dataprods/tdsim.Rdata") 


# 4.3. Visualisation ----
source("Code_aid.R"); aid$load("tdsim", "tdpar")
tdsimA <- apply(tdsim, 1:3, quantileMean, probs=c(0.3,0.5,0.7), na.rm=TRUE)
tdsimA <- tdsimA[,,,1:15] # remove temps with NA only
tdsimt <- tdsimt[1:15]

pdf("fig/fig4.pdf", height=4, pointsize=11) 
layout(matrix(c(1:3, rep(4,3)), ncol=2), widths=c(4,6))
par(mar=c(0,3,0,0.3), oma=c(3.5,0,0.1,0), mgp=c(2.1, 0.8,0), las=1, lend=1, cex=1)
# plot temperature dependent parameters:
leg <- function(i) {abline(lm(y~x, data=xys[[i]]), col="orange") 
     legend("topleft", legend=dimnames(tdpar)[[1]][i], inset=c(-0.0, -0.0), bty="n")}
xysplot <- function(i) plot(xys[[i]][,1],#+rnorm(994,sd=0.2), 
                            xys[[i]][,2], 
                            xlim=c(4,16), pch=16, col=addAlpha(1,0.2), ann=FALSE, axes=F)
xysplot(1); axis(2,at=seq(-0.5,-0.1,0.2)) ; leg(1); box()
xysplot(2); axis(2,at=seq(0.5,1.5,0.5)) ; leg(2); box()
xysplot(3); axis(2,at=seq(0.4,0.8,0.2)) ; leg(3); box()
axis(1)
rm(leg, xysplot)
title(xlab="Dewpoint temperature bin midpoint  [\U{00B0}C]", outer=TRUE)
##
par(mar=c(0,3.5,0,0.3))
plot(1, type="n", xlim=c(9,19.3), ylim=c(7,35), log="y", xlab="", yaxt="n",
     xaxs="i", ylab="GPD random sample 99.9% quantile  [mm/h]")
logAxis(2)
ciBand(yu=10^tdsimA["70%","empirical","99.9%",], nastars=FALSE,
       ym=10^tdsimA["50%","empirical","99.9%",],
       yl=10^tdsimA["30%","empirical","99.9%",], x=tdsimt, colm="green3", add=TRUE)
ciBand(yu=10^tdsimA["70%","gpa","99.9%",], nastars=FALSE,
       ym=10^tdsimA["50%","gpa","99.9%",],
       yl=10^tdsimA["30%","gpa","99.9%",], x=tdsimt, colm="blue", add=TRUE)
aid$cc_lines(NA, mainargs=list(col=2, lty=2))
tplot <- seq(5,21,1)
lines(tplot, sapply(tplot, function(t) 10^lmomco::quagpa(f=0.999, tdsimp(t))), lwd=2, col="orange")
legend("topleft", c("CC-scaling", "Real value from parameters in the left panel", 
                    "Parametric quantile", 
                    "Empirical quantile", "Central 40% of 1000 simulations"),
       lwd=c(1,2,2,2,11), lty=c(3,1,1,1), col=c(2,"orange","blue","green3",8), bg="white", cex=0.8)
text(20, c(1.25, 1.45), c("empirical","parametric") )
dev.off()



