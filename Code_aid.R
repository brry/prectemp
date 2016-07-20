# Clausius Clapeyron relationship for saturation vapor pressure
# Magnus-Tetens approximation using the August-Roche-Magnus formula
# shifted up and down to fit graphic range
cc_lines <- function(
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



cc_outlier <- function(temp)
  {
  temp[temp>10] <- NA
  vals <- cc_lines(3, anf=-16, ratecc=F, maincc=F)[1:28,]
  approx(vals$x, vals$y, xout=temp)$y
  }


