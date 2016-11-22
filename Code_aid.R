
aid <- list(
# cc_lines ---------------------------------------------------------------------
# Clausius Clapeyron relationship for saturation vapor pressure
# Magnus-Tetens approximation using the August-Roche-Magnus formula
# shifted up and down to fit graphic range
cc_lines = function(
startpunkte,  # Y values of CC-rate lines at location 'anfang'
anfang,       # leftmost x value. DEFAULT: left end of graph
maincc=TRUE,  # add line of CC-relationship?
mainargs=NULL,# List of arguments passed to lines for main CC-line
ratecc=TRUE,  # add line(s) of CC-rate parallel to main CC line?
lty=3,        # Line type of rate lines
...)
  {
  if(missing(anfang)) anfang <- round(par("usr")[1])
  temp <- seq(anfang,42, len=100)
  # http://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#Meteorology_and_climatology
  cc <- function(temp) 6.1094*exp(17.625*temp/(temp+243.04))
  # Change rate (6 to 7% per ?C)
  changerate <- diff(cc(anfang:43))/cc(anfang:42) + 1
  if(missing(startpunkte)) startpunkte <- if(par("ylog")) 10^(seq( -1.95, 2, len=13)) else
                                                         1.1^(seq(-15,40, len=14))
  for(startp in startpunkte)
  if(ratecc) lines(anfang:43, cumprod(c(startp, changerate)), lty=lty, ...)#col=8)
  if(maincc) do.call(lines, owa(list(x=temp, y=cc(temp)), mainargs))#, lwd=2)
  # return values of ccrate for last startpunkte
  invisible(data.frame(x=anfang:43, y=cumprod(c(startp, changerate))))
  }
,

# dewtemp ----------------------------------------------------------------------
# dew point temperature estimate from air temperature and relative humidity
dewtemp=function(airtemp, relhum)
  {
  #browseURL("https://en.wikipedia.org/wiki/Dew_point#Calculating_the_dew_point")
  # Buck, A. L. (1981), "New equations for computing vapor pressure and enhancement factor", 
  # J. Appl. Meteorol. 20: 1527-1532
  a <- 6.1121
  b <- 17.368
  c <- 238.88
  gamfun <- function(at, rh) log(rh/100)+ b*at/(c+at) 
  c*gamfun(airtemp, relhum) / (b-gamfun(airtemp, relhum))
  }
,

dewtemp2=function(airtemp, relhum)
  {
  #browseURL("https://de.wikipedia.org/wiki/Taupunkt#Abh.C3.A4ngigkeit_der_Taupunkttemperatur_von_relativer_Luftfeuchtigkeit_und_Lufttemperatur")
  K1 <- 6.112 ;   K2 <- 17.62 ;   K3 <- 243.12
  K3*(K2*airtemp/(K3+airtemp)+log(relhum/100))/(K2*K3/(K3+airtemp)-log(relhum/100))
  }

,

# mid --------------------------------------------------------------------------
# midpoints of temperature bins
mid = seq(5, 35, by=0.1)
,

# probs, probcols --------------------------------------------------------------
# Colors and values for Quantiles
probcols = c("orange","forestgreen","darkblue","red")
,
probs = c(0.9, 0.99, 0.999, 0.9999)
,

# stationplot ------------------------------------------------------------------

stationplot = function(   # set up empty plot for station i with map (for section 2.1.)
i, # station ID
meta,
map,
xlab="Dew point temperature (mean of preceding 5 hours)  [ \U00B0 C]",
xlim=c(-5,25),
ylim=c(0.5,70),
onlymap=FALSE
)
  {
  plot(1, type="n", xlim=xlim, ylim=ylim, log="y", yaxt="n",
      xlab=xlab, ylab="Precipitation  [mm/h]", main=meta$name[i])
  if(!onlymap) 
  {
  title(main=paste0(meta$ele[i]," meter asl\n", meta$id[i]," ID DWD\n",i," ID berry "),
        adj=1, cex.main=1, font.main=1)
  # station IDs in vicinity
  dist <- OSMscale::earthDist(lat, long, data=meta, i=i)
  closeby <- order(dist)[2:sum(dist < 80)]
  title(main=paste0("\n\n\n           within 80 km: ", toString(closeby)),
        adj=0, cex.main=1, font.main=1)
  }
  logAxis(2)
  aid$cc_lines(NA)
  smallPlot({
       plot(map, type="l", axes=F, ann=F)
       colPoints(long, lat, ele, data=meta, cex=0.6, col=seqPal(100, gb=T), legend=F)
       points(lat~long, data=meta[i,], lwd=2, col="red")
       },
       x1=0, x2=0.17, y1=0.70, y2=1, mar=rep(0,4), bg="white")
  } # end station plot

# xxx --------------------------------------------------------------

)
