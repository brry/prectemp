# Examplatory analysis of temperature-dependency of rainfall intensity
# With now freely available (but only 20 years) time series from DWD FTP-server
# Berry Boessenkool, July 2016
browseURL("http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183")

# 0. Packages
# 1. Download and process data from DWD server
# 1.1. Metadata and actually available files
# 1.2. Download files with longest records
# 2. Hourly Prec-Temp relationship


# 0. Packages ------------------------------------------------------------------
if(FALSE){ # download and install the packages only once
install.packages("extremeStat")
install.packages("RCurl")
source("http://raw.githubusercontent.com/brry/berryFunctions/master/R/instGit.R")
instGit("brry/berryFunctions") # must be version >= 1.10.22 (2016-07-19)
}

library("berryFunctions")


# 1. Download and process data from DWD server ---------------------------------
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
rownames(metadata) <- NULL

# sort stations for longest time series:
metadata$T_dur <- as.numeric(substr(metadata$T_name,32,39))/1e4 -
                  as.numeric(substr(metadata$T_name,23,30))/1e4
metadata$P_dur <- as.numeric(substr(metadata$P_name,32,39))/1e4 -
                  as.numeric(substr(metadata$P_name,23,30))/1e4
metadata$m_dur <- pmin(metadata$T_dur,metadata$P_dur)
metadata$T_dur_Beschr <- metadata$T_bis_datum/1e4 - metadata$T_von_datum/1e4
metadata$P_dur_Beschr <- metadata$P_bis_datum/1e4 - metadata$P_von_datum/1e4
metadata <- sortDF(metadata, "m_dur", decreasing=TRUE)
metadata$missing <- NA  # number of hours missing after merging P and T

# Save metadata:
save(metadata, file="metadata.Rdata")
View(metadata)
colPoints(geoLaenge, geoBreite, m_dur, data=metadata, add=F, asp=1.5)
colPointsHist(metadata$m_dur, mar=c(6,9,0,0), right=F, breaks=-1:20+0.5,
              nbins=22, bb=-1:20+0.5, bg="transparent")


# 1.2. Download files with longest records -------------------------------------
load("metadata.Rdata")
PTdown <- function(i)
  {
  Sys.sleep(runif(1,0,5)) # to avoid getting kicked off the FTP
  temp <- dataDWD(metadata[i,"T_name"], base2="hourly/air_temperature/historical", format=NULL, quiet=TRUE)
  Sys.sleep(runif(1,0,5))
  prec <- dataDWD(metadata[i,"P_name"], base2="hourly/precipitation/historical", format=NULL, quiet=TRUE)
  PT <- merge(prec, temp, by="MESS_DATUM")
  PT
  }

##PT <- pblapply( 1:10, PTdown) # 2 minutes    # already done, thus commented out
##save(PT, file="PT.Rdata")
load("PT.Rdata")
##PT <- c(PT, pblapply(11:20, PTdown)) # 2 minutes
##PT <- c(PT, pblapply(21:30, PTdown)) # 2 minutes
PT <- c(PT, pblapply(31:40, PTdown))

# In case you miss saving, but download has been complete:
PTread <- function(i)
  {
  temp <- readDWD(paste0("DWDdata/", sub(".zip", "", metadata[i,"T_name"])), format=NULL)
  prec <- readDWD(paste0("DWDdata/", sub(".zip", "", metadata[i,"P_name"])), format=NULL)
  PT <- merge(prec, temp, by="MESS_DATUM")
  PT
  }
###PT <- c(PT, pblapply(21:30, PTread)) # 1 minute

metadata$missing[1:length(PT)] <- sapply(PT, function(x){
  rng <- as.POSIXct(as.character(headtail(x$MESS_DATUM)), format="%Y%m%d%H", tz="GMT")
  length(seq.POSIXt(rng[1], rng[2], by="hour")) - nrow(x)})

save(metadata, file="metadata.Rdata")
save(PT, file="PT.Rdata")



# 2. Hourly Prec-Temp relationship ---------------------------------------------
load("metadata.Rdata"); load("DWDdata/PT.Rdata")
source("Code_aid.R") # cc_lines
library("extremeStat")

range(sapply(PT, function(x) max(x$NIEDERSCHLAGSHOEHE, na.rm=T))) #

pdf("RawData.pdf", height=5)
dummy <- pblapply(seq_along(PT), function(i){
  x <- PT[[i]]
  plot(1, type="n", xlim=c(-10,30), ylim=c(0.5,70), log="y", yaxt="n",
      xlab="Temperature  [°C]", ylab="Precipitation  [mm/h]",
      main=metadata$Stationsname[i])
  logAxis(2)
  cc_lines(NA)
  smallPlot({plot(geoBreite~geoLaenge, data=metadata, asp=1.5, pch=16, col="gray", axes=F, ann=F);
             points(geoBreite~geoLaenge, data=metadata[i,], pch=16)},
             x=c(0,20), y=c(70,100), mar=rep(0,4), bg="white")
  points(NIEDERSCHLAGSHOEHE~LUFTTEMPERATUR, data=x[x$NIEDERSCHLAGSHOEHE>=0.5,],
        pch=16, col=addAlpha(1))
  })
rm(dummy)
dev.off()
















