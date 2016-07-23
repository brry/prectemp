# Clausius Clapeyron relationship for saturation vapor pressure
# Magnus-Tetens approximation using the August-Roche-Magnus formula
# shifted up and down to fit graphic range

aid <- list(
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
  # Change rate (6 to 7% per °C)
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

cc_outlier = function(temp)
  {
  temp[temp>10] <- NA
  vals <- aid$cc_lines(3, anf=-16, ratecc=F, maincc=F)[1:28,]
  approx(vals$x, vals$y, xout=temp)$y
  }
,

# midpoints of temperature bins
mid = seq(5, 35, by=0.1)
,

# Colors and values for Quantiles
probcols = c("orange","forestgreen","darkblue","red")
,
probs = c(0.9, 0.99, 0.999, 0.9999)
,

stationplot = function(
i, # station ID
metadata,
map,
xlab="Temperature mean of preceding 5 hours  [°C]",
xlim=c(-15,35),
ylim=c(0.5,70)
)
  {
  plot(1, type="n", xlim=xlim, ylim=ylim, log="y", yaxt="n",
      xlab=xlab, ylab="Precipitation  [mm/h]", main=metadata$Stationsname[i])
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
  aid$cc_lines(NA)
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
  } # end station plot

)
