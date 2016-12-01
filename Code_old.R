# Analysis of temperature-dependency of rainfall intensity
# Old stuff, may not run without changes anymore
# Berry Boessenkool, 2016


# Outlier examination -----------------------------------------------------
load("PT.Rdata") ; source("Code_aid.R")
cc_outlier = function(temp)
  {
  temp[temp>10] <- NA
  vals <- aid$cc_lines(3, anf=-16, ratecc=F, maincc=F)[1:28,]
  approx(vals$x, vals$y, xout=temp)$y
  }
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

# No outliers are removed. These might be heavy snowstorms or data encoding errors.
# The PT-quantile computation will be restricted to >5 ?C
# Only station 25 Muenchen has a weird outlier there.


# old distribution weights ----------------------------------------------------
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



n <- seq(200, 300, len=4)
qn <- function(n) vapply(n, function(nn)
  {
  #sinkfile <- file("log.txt", open="wt")  
  #sink(sinkfile, append=TRUE)
  #sink(sinkfile, append=TRUE, type="message")
  
  d <- capture.output( 
    message("nn=", nn),
    distLquantile(sample(PREC, 30), probs=aid$probs, truncate=0.8, 
                     addinfo=TRUE, weightc=NA, ssquiet=TRUE, time=FALSE, 
                     progbars=FALSE, order=FALSE),  file="log.txt", append=TRUE, type="message")
  #sink(type="message")
  #sink() 
  d
  }, FUN.VALUE=array(0, dim=c(38,4)))


