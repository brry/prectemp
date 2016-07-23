# Examplatory analysis of temperature-dependency of rainfall intensity
# With now freely available (but only 20 years) time series from DWD FTP-server
# Berry Boessenkool, July 2016
browseURL("http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183")

# Each section should run independently in a clean R session
# Auxiliary functions often required are in Code_aid.R

# 0. Packages

# 1. Download and process data from DWD server
# 1.1. Metadata and actually available files
# 1.2. Download files with longest records
# 1.3. Reduce filesize and add temp5

# 2. Hourly Prec-Temp relationship
# 2.1. Raw data visualisation
# 2.2. Outlier examination
# 2.3. Distribution weights
# 2.4. PT-quantiles computation
# 2.5. PT-quantiles visualisation


# 0. Packages ------------------------------------------------------------------
if(FALSE){ # You need to download and install the packages only once
install.packages("extremeStat")
install.packages("pblapply")
install.packages("RCurl")
install.packages("maps")
install.packages("mapdata")
source("http://raw.githubusercontent.com/brry/berryFunctions/master/R/instGit.R")
instGit("brry/berryFunctions") # must be version >= 1.10.25 (2016-07-20)
instGit("brry/extremeStat") # must be version >= 0.5.21 (2016-07-23)
instGit("brry/OSMscale")
}


# 1. Download and process data from DWD server ---------------------------------
library("berryFunctions")
# 1.1. Metadata and actually available files -----------------------------------
metatemp <- dataDWD("TU_Stundenwerte_Beschreibung_Stationen.txt",
                    base2="hourly/air_temperature/historical")
metaprec <- dataDWD("RR_Stundenwerte_Beschreibung_Stationen.txt", 
                    base2="hourly/precipitation/historical")
availabletemp <- dataDWD("", meta=2, base2="hourly/air_temperature/historical")
availableprec <- dataDWD("", meta=2, base2="hourly/precipitation/historical")

# select stations with file available for both prec and temp:
availabletemp2 <- suppressWarnings(as.integer(substr(availabletemp, 17,21))) # NA for metadata files
availableprec2 <- suppressWarnings(as.integer(substr(availableprec, 17,21)))
metatemp2 <- metatemp[metatemp$Stations_id %in% availabletemp2 &
                      metatemp$Stations_id %in% availableprec2 ,]
metaprec2 <- metaprec[metaprec$Stations_id %in% availabletemp2 &
                      metaprec$Stations_id %in% availableprec2 ,]
metatemp2 <- metatemp2[metatemp2$Stations_id %in% metaprec2$Stations_id,]
metaprec2 <- metaprec2[metaprec2$Stations_id %in% metatemp2$Stations_id,]

# sort by ID and combine
metatemp2 <- sortDF(metatemp2, "Stations_id")
metaprec2 <- sortDF(metaprec2, "Stations_id")
stopifnot(all(metatemp2$Stations_id                == metaprec2$Stations_id))
stopifnot(all(metatemp2$Stationshoehe              == metaprec2$Stationshoehe))
stopifnot(all(metatemp2$geoBreite                  == metaprec2$geoBreite))
stopifnot(all(metatemp2$geoLaenge                  == metaprec2$geoLaenge))
stopifnot(all(as.character(metatemp2$Stationsname) == as.character(metaprec2$Stationsname)))
stopifnot(all(as.character(metatemp2$Bundesland)   == as.character(metaprec2$Bundesland)))
colnames(metatemp2)[2:3] <- paste0("T_", colnames(metatemp2)[2:3])
colnames(metaprec2)[2:3] <- paste0("P_", colnames(metaprec2)[2:3])
metadata <- cbind(metatemp2[,c(1,4:8,2:3)], metaprec2[,2:3])

# Add filenames:
metadata$T_name <- availabletemp[match(metadata$Stations_id, availabletemp2)]
metadata$P_name <- availableprec[match(metadata$Stations_id, availableprec2)]

# sort stations for longest time series:
metadata$T_dur <- as.numeric(substr(metadata$T_name,32,39))/1e4 -
                  as.numeric(substr(metadata$T_name,23,30))/1e4
metadata$P_dur <- as.numeric(substr(metadata$P_name,32,39))/1e4 -
                  as.numeric(substr(metadata$P_name,23,30))/1e4
metadata$m_dur <- pmin(metadata$T_dur,metadata$P_dur)
metadata$T_dur_Beschr <- metadata$T_bis_datum/1e4 - metadata$T_von_datum/1e4
metadata$P_dur_Beschr <- metadata$P_bis_datum/1e4 - metadata$P_von_datum/1e4
metadata <- sortDF(metadata, "m_dur", decreasing=TRUE)


# Add UTM-coordinates (for distance computations)
coord <- OSMscale::projectPoints(metadata$geoBreite, metadata$geoLaenge)
metadata$utm32_x <- coord[,1]
metadata$utm32_y <- coord[,2]
rm(coord)

# Save metadata:
rownames(metadata) <- NULL
save(metadata, file="metadata.Rdata")
##View(metadata)
colPoints(geoLaenge, geoBreite, m_dur, data=metadata, add=F, asp=1.5)
colPointsHist(metadata$m_dur, mar=c(6,9,0,0), right=F, breaks=-1:20+0.5,
              nbins=22, bb=-1:20+0.5, bg="transparent")


# 1.2. Download files with longest records -------------------------------------
load("metadata.Rdata")
PTdown <- function(i)
  {
  Sys.sleep(runif(1,0,5)) # to avoid getting kicked off the FTP
  dataDWD(metadata[i,"T_name"], base2="hourly/air_temperature/historical",
                  read=FALSE, format=NULL, quiet=TRUE)
  Sys.sleep(runif(1,0,5))
  dataDWD(metadata[i,"P_name"], base2="hourly/precipitation/historical",
                  read=FALSE, format=NULL, quiet=TRUE)
  }

##dummy <- pblapply(1:150, PTdown) # 1 minute for 10 stations (mostly sleep time)
## # Commented out so I don't acccidentally rerun that.
## # Also, I did this in batches. Never tested how sensitive the FTP is...
rm(dummy, PTdown)

# read and process files:
# Note: missing entries in rainfall data are added with NA-rows
PTread <- function(i)
  {
  temp <- readDWD(paste0("DWDdata/", metadata[i,"T_name"]), format=NULL)
  prec <- readDWD(paste0("DWDdata/", metadata[i,"P_name"]), format=NULL)
  # all hours in time window:
  rng <- as.POSIXct(as.character(headtail(prec$MESS_DATUM)), format="%Y%m%d%H", tz="GMT")
  hours <- as.integer(format(seq.POSIXt(rng[1], rng[2], by="hour"), "%Y%m%d%H"))
  prec <- merge(prec, data.frame(MESS_DATUM=hours), all=TRUE)
  # merge prec and temp:
  PT <- merge(prec, temp, by="MESS_DATUM", all.x=TRUE, all.y=FALSE)
  PT
  }

PT1 <- pblapply(  1: 50, PTread) ;  save(PT1, file="PT1.Rdata") # 06m 05s
PT2 <- pblapply( 51:100, PTread) ;  save(PT2, file="PT2.Rdata") # 05m 55s
PT3 <- pblapply(101:150, PTread) ;  save(PT3, file="PT3.Rdata") # 05m 00s


# 1.3. Reduce filesize and add temp5 -------------------------------------------
load("PT1.Rdata") ; load("PT2.Rdata"); load("PT3.Rdata"); load("metadata.Rdata")
PT <- c(PT1,PT2,PT3)
metadata$missing <- NA # number of hours missing in precipitation data
metadata$missing[1:150] <- sapply(PT, function(x) sum(
                         is.na(x$NIEDERSCHLAGSHOEHE) | is.na(x$LUFTTEMPERATUR)))
save(metadata, file="metadata.Rdata")
which.max(metadata$missing) # 131
logHist(metadata$missing, breaks=30, logargs=list(base=1, exponent=6), col="yellow3")
hist(sapply(PT, nrow), breaks=30, col="tomato")

# this takes about 15 minutes for 150 stations:
PT <- pblapply(PT, function(x) {
      x <- x[, c("MESS_DATUM", "LUFTTEMPERATUR", "NIEDERSCHLAGSHOEHE")]
      colnames(x) <- c("date", "temp", "prec")
      # temperature of the preceding 5 hours (excluding the hour of rainfall):
      x$temp5 <- sapply(1:nrow(x), function(i) mean(x[pmax(1,i-5):(i-1),"temp"], na.rm=TRUE))
      # keep only rainfall values >0.5 mm/h:
      x <- x[!is.na(x$prec) & x$prec>=0.5, ]
      x})

save(PT, file="PT.Rdata")



# 2. Hourly Prec-Temp relationship ---------------------------------------------
# 2.1. Raw data visualisation --------------------------------------------------
load("metadata.Rdata"); load("PT.Rdata"); source("Code_aid.R") # cc_lines
range(sapply(PT, function(x)   max(x$prec, na.rm=T))) # 19.5 80.8
range(sapply(PT, function(x) range(x$temp, na.rm=T))) # -16.4  34.6
range(sapply(PT, function(x) range(x$temp5,na.rm=T))) # -16.9  38.8
hist(sapply(PT, nrow), breaks=30, col="salmon")

library("mapdata")
map <- maps::map('worldHires','Germany')
dev.off()

# Plot PT-graph with current and preceding temp (30secs each):
for(t5 in 1:2){
pdf(c("RawData.pdf","RawData_T5.pdf")[t5], height=5)
dummy <- pblapply(1:150, function(i){
  x <- PT[[i]]    # pbapply(order(metadata$Stationshoehe[1:150]), ...
  aid$stationplot(i, metadata, map, if(t5==1) xlab="Temperature  [°C]",
                  ylim=c(2.5,70) )
  points(x=x[x$prec>2,c("temp","temp5")[t5]], y=x$prec[x$prec>2], pch=16, col=addAlpha(1))
  lines(c(-16:10,10), c(aid$cc_outlier(-16:10),100), lty=3)
  text(-10,50, 50)
  })
dev.off()
}
rm(dummy, t5)

rm(map)


# 2.2. Outlier examination -----------------------------------------------------
load("PT1.Rdata") ; load("PT2.Rdata"); load("PT3.Rdata"); source("Code_aid.R")
PT_all <- c(PT1,PT2,PT3)  ; rm(PT1,PT2,PT3)
# Outlier definition: point above line in previous images
PT_outlier <- pblapply(seq_along(PT), function(i){
              x <- PT[[i]]
              index <- which(x$prec>aid$cc_outlier(x$temp5))
              if(length(index)<1) return(NA)
              y <- PT_all[[i]]
              y$outlier <- 0
              y$IDberry <- i
              y[y$MESS_DATUM %in% x$date[index], "outlier"] <- 1
              index2 <- c(sapply(which(y$outlier==1), "+", -5:5))
              index2 <- sort(unique(index2))
              index2 <- index2[index2>0]
              y_out <- y[index2 ,c(1,11,5,14:15, 2:4,6:10,12:13) ]
              y_out$MESS_DATUM <- sapply(y_out$MESS_DATUM, function(d)
                     paste(substring(d, c(1,5,7,9), c(4,6,8,10)), collapse="-"))
              substr(y_out$MESS_DATUM, 11,11) <- "_"
              insert <- which(diff(index2)!=1)
              if(length(insert)>1) y_out <- insertRows(y_out, insert)
              write.table(y_out, file=paste0("outlier/outlier_",i,".txt"),
                          quote=F, row.names=F) # sep="\t",
              y_out})

files <- dir("outlier", pattern=".txt", full=T)
names <- sapply(strsplit(files,"_"),"[",2)
names <- sapply(strsplit(names,".",fix=T),"[",1)
files <- files[order(as.integer(names)) ]

combineFiles(files, outFile="outlier/outlier all", sep=" ", names=F)
# that is copied into excel file and formatted
rm(files, names)

# No outliers are removed. These might be heavy snowsotmrs or data encoding errors.
# The PT-quantile computation will be restricted to >5 °C
# Only station 25 Muenchen has a weird outlier there.


# 2.3. distribution weights ----------------------------------------------------
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


# 2.4. PT-quantiles computation ------------------------------------------------
load("PT.Rdata"); load("cweights.Rdata"); source("Code_aid.R") # mid, probs

# long computing time (1 minute per station)
library(parallel) # for parallel lapply execution
cl <- makePSOCKcluster(detectCores(), outfile="PTQlog.txt")
clusterExport(cl, c("PT","aid", "cweights"))
dummy <- clusterEvalQ(cl, library(extremeStat))
begintime <- Sys.time(); begintime
PTQ <- parLapply(cl, 1:150, function(i)
  {
  x <- PT[[i]]
  # Quantile estimates per temperature bin
  binQ <- lapply(aid$mid, function(t) {
    seldat <- x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)]
    distLquantile(seldat[!is.na(seldat)], probs=aid$probs, truncate=0.8, addinfo=TRUE,
                  weightc=cweights, order=FALSE, ssquiet=TRUE, time=FALSE, progbars=FALSE)
    })
  message("--------------------------------\n", Sys.time(),
          "\nbinQ computed for station ID ", i, "\n--------------------------------")
  # Transform into array for faster subsetting:
  binQ2 <- array(unlist(binQ), dim=c(37, 4, 301),
       dimnames=list(distr=rownames(binQ[[1]]), prob=colnames(binQ[[1]]), temp=aid$mid))
  return(binQ2)
  })
stopCluster(cl)
save(PTQ, file="PTQ.Rdata")
Sys.time() - begintime # on 4 cores 31 minutes
rm(cl, begintime, dummy)

d <- readLines("PTQlog.txt")
d <- as.numeric(substr(d[grepl("binQ",d)], 30,100))
paste(round(length(d)/150*100,1), "% done")
d

which(!1:150 %in% d)
# should be empty, but 113 is missing

rm(d)
length(PTQ) # 150 stations
str(PTQ[[1]]) # each at 301 temperature bins
aid$mid[200] # for 24.9 °C:
PTQ[[1]][,,200]

head(PT[[113]])
n113 <- sapply(aid$mid, function(t) sum(PT[[113]]$temp5>(t-1) & PT[[113]]$temp5<=(t+1)) )
all(n113 == PTQ[[113]]["n_full",1,] )
rm(n113)


# 2.5. PT-quantiles visualisation ----------------------------------------------
load("PTQ.Rdata"); source("Code_aid.R")

pdf("PTQ.pdf", height=5)
#
message("Creating plot 1/3"); flush.console()
#
for(prob in c("90%", "99%", "99.9%", "99.99%"))
{
plot(1, type="n", xlim=c(5,30), ylim=c(2,90), log="y", yaxt="n", main=prob,
      xlab="Temperature mean of preceding 5 hours  [°C]", ylab="Precipitation  [mm/h]")
logAxis(2)
aid$cc_lines(NA)
dummy <- sapply(seq_along(PTQ), function(i){
     x <- PTQ[[i]]
     lines(aid$mid, x["weightedc",    prob, ], col=addAlpha("blue") )
     lines(aid$mid, x["quantileMean", prob, ], col=addAlpha("red") )
     })
}
#
message("Creating plot 2/3"); flush.console()
#
par(mfrow=c(1,2), mar=c(2,2,1.5,0.5), oma=c(1.5,1.5,1.5,0), mgp=c(2.5,0.6,0) )
for(prob in c("90%", "99%", "99.9%", "99.99%"))
{
for(type in 1:2)
{
dn <- c("quantileMean","weighted2")[type]
plot(1, type="n", xlim=c(8,28), ylim=c(5,130), log="y", yaxt="n", main=dn, ylab="", xlab="")
title(main=prob, xlab="Temperature mean of preceding 5 hours  [°C]",
      ylab="Precipitation  [mm/h]", outer=TRUE, line=0)
logAxis(2)
aid$cc_lines(NA)
stats <- sapply(PTQ, function(x){
            lines(aid$mid, x[dn, prob, ], col=addAlpha(c("red","blue")[type]) )
            x[dn, prob, ]})
stats <- replace(stats, stats>150, NA)  ###
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
lines(aid$mid, statav, lwd=2)
}}
#
message("Creating plot 3/3"); flush.console()
#
par(mfrow=c(1,1), mar=c(4,4,2,0.5), oma=c(0,0,0,0), mgp=c(2.5,0.6,0) )
dummy <- pblapply(dimnames(PTQ[[1]])$distr, function(dn){
ylim <- c(5,130)
cut <- 150
if(dn=="n_full")    {ylim <- c(5,2000); cut <- 1e5}
if(dn=="n")         {ylim <- c(1,400) ; cut <- 1e5}
if(dn=="threshold") {ylim <- c(0.5,25); cut <- 1e5}
plot(1, type="n", xlim=c(8,28), ylim=ylim, log="y", yaxt="n", main=dn,
      xlab="Temperature mean of preceding 5 hours  [°C]",
      ylab="Precipitation  99.9% quantile  [mm/h]")
logAxis(2)
if(!dn %in% c("n_full","n","threshold")) aid$cc_lines(NA)
stats <- sapply(PTQ, function(x){
            lines(aid$mid, x[dn, "99.9%", ], col=addAlpha(1) )
            x[dn, "99.9%", ]})
stats <- replace(stats, stats>cut, NA)
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
lines(aid$mid, statav, lwd=2, col="red")
if(dn=="n_full") abline(h=25, lwd=2, col="red")
if(dn=="n") abline(h=5, lwd=2, col="red")
})
rm(dummy, prob, dn, type, stats, statav)
dev.off()



load("PTQ.Rdata"); load("PT.Rdata"); load("metadata.Rdata"); source("Code_aid.R")
library("mapdata")
map <- maps::map('worldHires','Germany')
dev.off()

dn <- c("quantileMean","weighted2","gpa", "wak")
dc <- c("brown1", "deepskyblue4", "darkolivegreen4", "peru")
pdf("PTQ_stats.pdf", height=5)
dummy <- pblapply(1:150, function(i){
  aid$stationplot(i, metadata, map, xlim=c(8,28), ylim=c(5,130) )
  legend("bottomright", rev(dn), col=rev(dc), lwd=2, bg="white")
  points(PT[[i]][PT[[i]]$prec>5,c("temp5","prec")], pch=16, col=addAlpha(1))
  for(d in 1:4) lines(aid$mid, PTQ[[i]][dn[d], "99.9%", ], col=dc[d], lwd=2)
  })
rm(dn, dc, dummy)
dev.off()

