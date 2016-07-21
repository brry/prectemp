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
# 2.3. PT-quantiles computation
# 2.4. PT-quantiles aggregation
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
  plot(1, type="n", xlim=c(-15,35), ylim=c(0.5,70), log="y", yaxt="n",
      xlab=c("Temperature  [°C]","Temperature mean of preceding 5 hours  [°C]")[t5],
      ylab="Precipitation  [mm/h]", main=metadata$Stationsname[i])
  title(main=paste0(metadata$Stationshoehe[i]," meter asl\n", metadata$Stations_id[i],
                   " ID DWD\n", i, " ID berry\n",
                   round(metadata$missing[i]/(metadata$m_dur[i]*365*24)*100,1),
                   "% missing"), adj=1, cex.main=1, font.main=1)
  # station IDs in vicinity
  utmx <- metadata$utm32_x[1:150]
  utmy <- metadata$utm32_y[1:150]
  # d <- lapply(1:150, function(i) {d <- distance(utmx, utmy, utmx[i], utmy[i])
  #                                 order(d)[2:sum(d < 80*1000)]})
  # hist(sapply(d, length)) # for 80 km mostly between 4 and 8 (quartiles)
  dist <- distance(utmx, utmy, utmx[i], utmy[i])
  closeby <- order(dist)[2:sum(dist < 80*1000)]
  title(main=paste0("\n\n\n           within 80 km: ", toString(closeby)),
        adj=0, cex.main=1, font.main=1)
  logAxis(2)
  cc_lines(NA)
  smallPlot({
             plot(map, type="l", axes=F, ann=F)
             points(geoBreite~geoLaenge, data=metadata, col="gray95", pch=16, cex=0.6)
             #points(geoBreite~geoLaenge, data=metadata[1:150,], col="gray85", pch=16)
             lines(map)
             colPoints(geoLaenge, geoBreite, Stationshoehe, data=metadata[1:150,],
                       cex=0.6, col=seqPal(150, gb=T), legend=F)
             points(geoBreite~geoLaenge, data=metadata[i,], lwd=2, col="red")
             },
             x=c(0,17), y=c(70,100), mar=rep(0,4), bg="white")
  points(x=x[,c("temp","temp5")[t5]], y=x$prec, pch=16, col=addAlpha(1))
  lines(c(-16:10,10), c(cc_outlier(-16:10),100), lty=3)
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
              index <- which(x$prec>cc_outlier(x$temp5))
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


# 2.3. PT-quantiles computation ------------------------------------------------
load("PT.Rdata"); source("Code_aid.R") # mid, probs
library(extremeStat)

# long computing time (1 minute per station)
library(parallel) # for parallel lapply execution
cl <- makePSOCKcluster(detectCores())
clusterExport(cl, c("PT","mid","probs"))
clusterEvalQ(cl, library(extremeStat))
begintime <- Sys.time(); begintime
PTQ <- parLapply(cl, PT, function(x)
  {
  # Quantile estimates per temperature bin
  binQ <- lapply(mid, function(t) {
    seldat <- x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)]
    distLquantile(na.omit(seldat), probs=probs, truncate=0.8, addinfo=TRUE,
                  quiet=TRUE, time=FALSE, progbars=FALSE)
    })
  return(binQ)
  })
stopCluster(cl)
Sys.time() - begintime # on 8 cores __prognosed ca 16__ minutes
save(PTQ, file="PTQ.Rdata")
rm(cl, begintime)



# 2.4. PT-quantiles aggregation ------------------------------------------------
load("PT.Rdata"); load("PTQ.Rdata")
length(PTQ) # 150 stations
length(PTQ[[1]]) # each at 301 temperature bins
mid[150] # for 19.9 °C:
PTQ[[1]][[150]]


# 2.5. PT-quantiles visualisation ----------------------------------------------












