# Analysis of temperature-dependency of rainfall intensity
# With 142 freely available hourly time series (>15 years) from DWD FTP-server
# Berry Boessenkool, 2016
browseURL("http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183")

# Each section should run independently in a clean R session
# Frequently needed auxiliary functions are in Code_aid.R
# some code is commented out to avoid accidental calling (e.g. if time consuming)
# computing times are noted in minutes (e.g. # 5 min)
# Historically, section 3 came before section 2, but this is organized by content

# Contents:

# 0. Packages

# 1. Data 
# 1.1. Select DWD stations
# 1.2. Meta data weather stations
# 1.3. Download DWD station data
# 1.4. Dewpoint temperature

# 2. Sample size dependency
# 2.1. SSD all distributions
# 2.2. Distribution weights
# 2.3. SSD GPD 
# 2.4. Truncation dependency 
# 2.5. tempdep Wakeby distribution

# 3. Hourly Prec-Temp relationship
# 3.1. Raw data visualisation
# 3.2. PT-quantiles computation
# 3.3. PT-quantiles visualization
# 3.4. PTQ per station


# 0. Packages ------------------------------------------------------------------
if(FALSE){ # You need to download and install the packages only once
packinst <- function(n) if(!requireNamespace(n, quietly=TRUE)) install.packages(n)
sapply(c("berryFunctions", "extremeStat", "pblapply", "maps", "gplots",
         "mapdata", "OSMscale", "RCurl"), packinst)
berryFunctions::instGit("brry/extremeStat")  # must be version >= 0.5.24 (2016-11-28)
berryFunctions::instGit("brry/rdwd") # not yet on CRAN must be >= 0.5.1  (2016-11-21)
berryFunctions::instGit("brry/OSMscale") #     must be version >= 0.3.11 (2016-11-22)
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
save(meta, file="meta.Rdata")

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


# Dew point temperature of previous 5 hours:
PT <- pblapply(PT_all, function(x) {                 # 10 seconds
  # Reduce filesize: column selection
  x <- x[, c("MESS_DATUM", "LUFTTEMPERATUR", "REL_FEUCHTE", "NIEDERSCHLAGSHOEHE")]
  colnames(x) <- c("date", "temp", "hum", "prec")
  # dew point temperature
  x$dewtemp <- aid$dewtemp(x$temp, x$hum) 
  # fill single missing NAs in dewtemp
  x$dtna <- approx(x$dewtemp, n=nrow(x))$y
  x$dtna[is.na(x$dewtemp[-1]) & is.na(x$dewtemp[-nrow(x)])] <- NA
  # dew point temperature of the preceding 5 hours (excluding the hour of rainfall):
  x$temp5 <- filter(x$dtna, c(0,rep(1/5,5)), sides=1)
  x$dtna <- NULL
  ### x$temp51 <- pbsapply(1:nrow(x), function(i) mean(x[pmax(1,i-5):(i-1),"dewtemp"], na.rm=TRUE)) # ca 20-30 secs per run!
  ### x$temp52 <- filter(x$dewtemp, c(0,rep(1/5,5)), sides=1) # ca 0.01 sec
  # Reduce filesize: keep only rainfall values >0.5 mm/h:
  x <- x[!is.na(x$prec) & !is.na(x$temp5) & x$prec>=0.5, ]
  x})
names(PT) <- names(PT_all)
rm(PT_all)
save(PT, file="PT.Rdata")



# 2. Sample size dependency ----------------------------------------------------

# 2.1. SSD all distributions ---------------------------------------------------

# we use custom weights from section 2.2
# ToDo: paper Fig 5 + 6_potsdamQn


# 2.2. Distribution weights ----------------------------------------------------

# we use a truncation of 80%, as will be examined in section 2.4.
load("PT.Rdata")
library(extremeStat)
# get weights for distribution functions:
dn <- c("exp", "gam", "gev", "glo", "gno", "gpa", "gum", "kap", "lap",
            "ln3", "nor", "pe3", "ray", "revgum", "rice", "wak", "wei")
RMSE <- pblapply(PT, function(x){   # computes 2 minutes
   sapply(c(5,10,15,20,25,30), function(t){
      d <- distLfit(x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)], truncate=0.8, plot=F, quiet=T)
      d$gof[dn, "RMSE"]})})
RMSE2 <- array(unlist(RMSE), dim=c(17,6,150))

w <- apply(RMSE2, 1:2, mean, na.rm=TRUE)
cweights <- max(rowMeans(w[,-6]))-rowMeans(w[,-6])
names(cweights) <- dn
cweights[c("revgum", "nor", "rice", "ray")] <- 0 #, "gam", "gum"
cweights <- cweights/sum(cweights)


pdf("Weights.pdf", height=5)
par(mfrow=c(1,1), mar=c(4,4,2,0.2), mgp=c(2.7,0.6,0), oma=c(0,0,0,0), las=1)
plot(1:17, type="n", xlim=lim0(0.15), yaxt="n", main="mean RMSE across stations",
     ylab="Distribution Function", xlab="RMSE (lower = better fit)")
cols <- seqPal(5, col=c("blue", "red"))
for(i in 1:5) lines(w[order(rowMeans(w[,-6])),i], 1:17, col=cols[i])
axis(2,1:17, dn[order(rowMeans(w[,-6]))])
legend("bottomright", legend=paste(c(5,10,15,20,25),"°C"), title="Temp bin midpoint",
       col=cols, lty=1)
cwplot <- cweights[order(rowMeans(w[,-6]))]
lines(replace(cwplot, cwplot==0, NA), 1:17, lwd=3, type="o", pch=16)
text(0.05, 13.8, "Weights")
#
plot(c(5,10,15,20,25,30), 1:6, type="n", ylim=lim0(0.1), las=1,
     main="min RMSE per temp", xlab="Temp bin midpoint", ylab="RMSE")
for(i in 1:150) points(c(5,10,15,20,25,30), apply(RMSE2, 2:3, min, na.rm=TRUE)[,i],
                       pch=16, col=addAlpha(1))
text(17.5, 0.01, "One dot per station")
dev.off()

save(cweights, file="cweights.Rdata")
rm(dn, RMSE, RMSE2, w, cols, cwplot, i)


# 2.3. SSD GPD -----------------------------------------------------------------
# ToDo: recreate paper Fig 6


# 2.4. Truncation dependency ---------------------------------------------------
# ToDo: recreate paper Fig 3 + 4 (new)


# 2.5. tempdep Wakeby distribution ---------------------------------------------
# ToDo: recreate paper Fig 7


# 3. Hourly Prec-Temp relationship ---------------------------------------------

# 3.1. Raw data visualisation --------------------------------------------------

load("meta.Rdata"); load("PT.Rdata"); source("Code_aid.R") # aid$cc_lines
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
dummy <- pblapply(1:142, function(i){  # pbapply(order(meta$ele), ...
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

load("PT.Rdata"); source("Code_aid.R") # mid, probs
load("cweights.Rdata")

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
save(PTQ, file="PTQ.Rdata")     # 30 mins
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
load("PTQ.Rdata"); source("Code_aid.R")

PTQplot <- function(
prob="",
main=prob,
xlab="Dewpoint temperature (mean of preceding 5 hours)  [ \U{00B0}C]", 
ylab=paste("Precipitation ", prob, "quantile  [mm/h]"),
outer=FALSE,
xlim=c(4.8,21),
ylim=c(2,130),
line=NA,
cc=TRUE,
...
)
{
plot(1, type="n", xlim=xlim, ylim=ylim, log="y", yaxt="n", xaxs="i", ann=FALSE)
title(main=main, xlab=xlab, ylab=ylab, outer=outer, line=line) 
logAxis(2)
if(cc) aid$cc_lines(NA)
}

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
stats <- replace(stats, stats>cut, NA)       ### make sure this is mentioned in paper!
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
statav
}

pdf("fig/PTQ.pdf", height=5)
# a. empirical and parametric quantiles in one plot ----------------------------
par(mar=c(4,4,2,0.5), mgp=c(2.5,0.6,0))
dn <- c("quantileMean", "weighted2")              ### later: use weightedc!
dc <- addAlpha(c("red", "blue"))
for(prob in c("90%", "99%", "99.9%", "99.99%"))
{
PTQplot(prob=prob)
dummy <- sapply(PTQ, function(x){
         for(d in 1:2) lines(aid$mid, x[dn[d], prob, ], col=dc[d]) })
}
# b. Qemp/par side by side -----------------------------------------------------
par(mfrow=c(1,2), mar=c(2,2,1.5,0.5), oma=c(1.5,1.5,1.5,0) )
for(prob in c("90%", "99%", "99.9%", "99.99%"))
{
for(d in 1:2)
{
PTQplot(prob=prob, outer=TRUE, line=0)
title(main=dn[d])
statav <- PTQlines(prob=prob, dn=dn[d], col=dc[d])
lines(aid$mid, statav, lwd=2)  
}}
#
# c. PTQ for each distribution function ----------------------------------------
par(mfrow=c(1,1), mar=c(4,4,2,0.5), oma=c(0,0,0,0) )
dummy <- pblapply(dimnames(PTQ[[1]])$distr, function(dn){ 
ylim <- c(5,130)
cut <- 150
if(dn=="n_full")    {ylim <- c(5,2000); cut <- 1e5}
if(dn=="n")         {ylim <- c(1, 400); cut <- 1e5}
if(dn=="threshold") {ylim <- c(0.5,25); cut <- 1e5}
PTQplot(prob="99.9%", ylim=ylim, main=dn, cc=!dn %in% c("n_full","n","threshold"))
statav <- PTQlines(prob="99.9%", dn=dn, cut=cut)
lines(aid$mid, statav, lwd=2, col="red")
if(dn=="n_full") abline(h=25, lwd=2, col="red")
if(dn=="n") abline(h=5, lwd=2, col="red")
})
rm(dummy, prob, d, dn, dc, statav)
dev.off()

# ToDo: recreate paper Fig 8


# 3.4. PTQ per station ---------------------------------------------------------

load("PTQ.Rdata"); load("PT.Rdata"); load("meta.Rdata"); source("Code_aid.R")
library("mapdata") ; map <- maps::map('worldHires','Germany') ; dev.off()

dn <- c("quantileMean","weighted2","gpa", "wak")
dc <- c("brown1", "deepskyblue4", "darkolivegreen4", "peru")
pdf("fig/PTQ_stats.pdf", height=5)
dummy <- pblapply(1:142, function(i){
  aid$stationplot(i, meta, map, xlim=c(5,21), ylim=c(5,130) )
  legend("topleft", rev(dn), col=rev(dc), lwd=2, inset=c(0.064,0), bg="white")
  points(PT[[i]][PT[[i]]$prec>5,c("temp5","prec")], pch=16, col=addAlpha(1))
  for(d in 1:4) lines(aid$mid, PTQ[[i]][dn[d], "99.9%", ], col=dc[d], lwd=2)
  })
rm(dn, dc, dummy)
dev.off()

