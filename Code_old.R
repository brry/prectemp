# Analysis of temperature-dependency of rainfall intensity
# Old stuff, may not run without changes anymore
# Berry Boessenkool, 2016



# gpa for random gamma values ----


dlf <- distLfit(log10(PREC), sel="gam", truncate=0.9)
dlf$parameter$gam #       alpha=11.2,  beta=0.05
GAMPAR <- list(type="gam", para=c(alpha=11.2, beta=0.05))
GAM <- 10^lmomco::rlmomco(1e4, GAMPAR)
10^lmomco::qlmomco(0.999, GAMPAR) # 16.6

dlf <- distLfit(log10(GAM), truncate=0.9)
plotLfit(dlf)
dlq <- distLquantile(dlf=dlf, probs=0.999, truncate=0.9, emp=F, list=T)
plotLquantile(dlq)

dlq80 <- distLquantile(log10(GAM), truncate=0.9, probs=0.999, weighted=F, gpd=F)
dlq90 <- distLquantile(log10(GAM), truncate=0.9, probs=0.999, weighted=F, gpd=F)

plot(10^dlq80[1:19,1], ylim=c(10,24)) ; points(19.5, 10^dlq80["gpa",1], pch=16)
plot(10^dlq90[1:19,1], ylim=c(10,24)) ; points(19.5, 10^dlq90["gpa",1], pch=16)

gamq <- pbsapply(aid$n, function(n)
  distLquantile(lmomco::rlmomco(n, GAMPAR), probs=0.999,
                truncate=0.9, sel="gpa", weight=F, gpd=F, quiet=T)[1:2,1])

plot(aid$n, 10^gamq[2,], type="l", col="green3", log="y", yaxt="n")
logAxis(2)
lines(aid$n, 10^gamq[1,], col=addAlpha("blue"))
abline(h=10^lmomco::qlmomco(0.999, GAMPAR), lty=3)



# 2.1. Select population -------------------------------------------------------

# select temperature bin for population for SSD experiments
# - range where there is no CC drop yet
# - distribution should be representative:
pdf("fig/binprec.pdf", height=5)
allprecs <- pblapply(seq(-11,23,0.2), function(temp)   # 15 secs
  {
  bin <- temp+c(-1,1)
  PREC <- unlist(lapply(PT, function(x) x[x$temp5>bin[1] & x$temp5<bin[2], "prec"] ))
  logHist(PREC, breaks=80, main=paste0("Histogramm of all hourly rainfall ",
          "records at 142 stations in ", formatC(round(temp,1), format="f", digits=1),
          "\U{00B0}C bin"), xlab="Precipitation  [mm/h]", 
          xlim=log10(c(0.5,80)), ylim=lim0(6), col="deepskyblue1", freq=FALSE)
  title(sub=format(length(PREC), big.mark="'"), adj=1)
  d <- density(log10(PREC))
  lines(d$x, d$y*4)
  })
dev.off()

# All 136k rainfall records between 10 and 12 degrees event dewpoint temperature 
PREC <- unlist(lapply(PT, function(x) x[x$temp5>10 & x$temp5<12, "prec"] ))
PREC <- unname(PREC)
hist(PREC, breaks=50, col="deepskyblue1")
logHist(PREC)
logHist(PREC, breaks=80)




# add custom weighted quantile estimates:
library(extremeStat)
old <- simQ[,,100,1]
new <- q_weighted(old, weightc=weights)
stopifnot(all(old==new, na.rm=TRUE))





which(simQ["kap", "99.9%", ,]>200, arr.ind=TRUE) # 2
kap <- as.vector(simQ["kap","99.9%",,])
val <- logSpaced(min=120, max=250, n=100, plot=F)
numlarger <- sapply(val, function(x) sum(kap>x, na.rm=T))
plot(val, numlarger, log="y", type="o", las=1)




# bias:
d_bias <- sapply(dimnames(simQA)[[2]][1:35], function(d) 
  rmse(simQA["50%",d,"99.9%",], rep(quantileMean(PREC, 0.999),512)))
w_bias <- d_bias["kap"]-d_bias[1:17]
w_bias[w_bias<0] <- 0
w_bias <- w_bias/sum(w_bias)

# goodness of fit:
d_gof <- apply(simQ[1:35,"RMSE",,], 1, mean, na.rm=TRUE)
w_gof <- d_gof[1:17]
w_gof <- w_gof/sum(w_gof)

# error rate:
simQ[1:35,,,][simQ[1:35,,,]>500] <- NA
d_error <- apply(simQ[1:35,1:4,,], 1, function(x) mean(is.na(x)))
w_error <- 0.02-d_error[1:17]
w_error[w_error<0] <- 0
w_error <- w_error/sum(w_error)

# total weights (for all dists, incl. GPD_):
d_weights <- data.frame(d_bias, d_gof, d_error)
d_weights <- d_weights[rownames(d_weights) !="weightedc",]
d_weights <- d_weights[rownames(d_weights) !="GPD_BAY_extRemes",]
d_weights <- apply(d_weights, 2, function(x) x/sum(x,na.rm=TRUE)*100)
d_weights <- as.data.frame(d_weights)
d_weights$mean <- rowMeans(d_weights, na.rm=TRUE)
d_weights <- sortDF(d_weights, "mean", decreasing=FALSE)



# _b. goodness of fit -----

#dlf <- extremeStat::distLfit(PREC, truncate=0.8)
#dlf$gof

load("dataprods/PT.Rdata")
PRECstats <- lapply(PT, function(x) x[x$temp5>10 & x$temp5<12, "prec"] )
dlfstats <- pblapply(PRECstats, distLfit, truncate=0.8, progbars=FALSE, 
                     time=FALSE, plot=FALSE) # 2 mins with plot, 1 min without
dlfgofs <- lapply(dlfstats, "[[", "gof")
dlfgofs <- do.call(rbind, dlfgofs)
dlfgofs$dist <- sapply(strsplit(rownames(dlfgofs), ".", fix=T), "[", 2)
dn <- unique(dlfgofs$dist)
dn <- sapply(dn, function(d) mean(dlfgofs[dlfgofs$dist==d, "RMSE"]) )
dn <- names(sort(dn))

pdf("fig/distribution_gofs_stations.pdf", height=5)
for(d in dn){
lh <- logHist(dlfgofs$RMSE, breaks=60, las=1, col=addAlpha("darkorange"), border=NA, main=d)
logHist(dlfgofs[dlfgofs$dist==d, "RMSE"], breaks=lh$breaks, col=1, logargs=list(xaxt="n"), add=TRUE)
}
rm(d)
dev.off()



load("sim/QN1.Rdata")
QN1[c("kap","ln3"),,1:20]
rm(QN1)



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


