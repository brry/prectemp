# Clausius Clapeyron relationship for saturation vapor pressure
# Magnus-Tetens approximation using the August-Roche-Magnus formula
# shifted up and down to fit graphic range
cc_lines <- function(startpunkte, anfang, maincc=TRUE, mainargs=NULL, lty=3, ...)
  {
  if(missing(anfang)) anfang <- round(par("usr")[1])
  temp <- seq(anfang,32, len=100)
  # http://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#Meteorology_and_climatology
  cc <- function(temp) 6.1094*exp(17.625*temp/(temp+243.04))
  # Change rate (6 to 7% per °C)
  changerate <- diff(cc(anfang:33))/cc(anfang:32) + 1
  if(missing(startpunkte)) startpunkte <- if(par("ylog")) 10^(seq( -1.95, 2, len=13)) else
                                                         1.1^(seq(-15,40, len=14))
  for(startp in startpunkte)
     lines(anfang:33, cumprod(c(startp, changerate)), lty=lty, ...)#col=8)
  if(maincc) do.call(lines, owa(list(x=temp, y=cc(temp)), mainargs))#, lwd=2)
  }



