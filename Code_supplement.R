# Analysis of temperature-dependency of rainfall intensity
# supplementary paper
# Berry Boessenkool, 2017
browseURL("http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183")



# _ SSD GPD implementations -----
source("Code_aid.R"); aid$load("simQA", "PREC", "weights")

dn5 <- dimnames(simQA)[[2]][c(33,23,24,26,25)]
dn6 <- dimnames(simQA)[[2]][27:32] ; dn6 <- sort(dn6)
col5 <- RColorBrewer::brewer.pal(5, "Set2") ; names(col5) <- dn5
col6 <- RColorBrewer::brewer.pal(6, "Set2") ; names(col6) <- dn6

pdf("fig/sup_SSD.pdf", height=3, width=3.5, pointsize=10)
#for(smooth in c(1,3,5,7,9,11,13,15)){
smooth <- 9; {
par(mfrow=c(1,2), mar=c(2,0,0.2,0.4), oma=c(1,3,0,0), mgp=c(1.8,0.7,0), las=1)
## panel 1
plot(1, type="n", xlim=c(25,830), ylim=log10(c(6,21)), xaxs="i", main="", yaxt="n",
       xlab="", ylab="")
logAxis(2)
for(d in dn6) lines(aid$n, movAv(simQA["50%",d,"99.9%",], smooth), col=col6[d])
abline(h=quantileMean(log10(PREC), probs=0.999), lty=3)
text(50, log10(19), "MLE", adj=0)
legend("bottomright", replace(dn6, 6, "GPD_MLE_Renext"), lwd=2, col=col6, bg="white", cex=0.65)
## panel 2
#par(mar=c(3,0,0.2,0.4))
plot(1, type="n", xlim=c(25,830), ylim=log10(c(6,21)), xaxs="i", main="", yaxt="n",
       xlab="", ylab="")
logAxis(2, labels=FALSE)
for(d in dn5) lines(aid$n, movAv(simQA["50%",d,"99.9%",], smooth), col=col5[d])
abline(h=quantileMean(log10(PREC), probs=0.999), lty=3)
text(50, log10(19), "LM / PWM", adj=0)
legend("bottomright", dn5, lwd=2, col=col5, bg="white", cex=0.65)
#
dens <- density(log10(PREC))
#lines(25+dens$y*1000, dens$x)
title(xlab="sample size n", outer=TRUE, mgp=c(-0.2,-1,0), xpd=NA)
title(ylab="Random sample 99.9% quantile  [mm/h]", outer=TRUE, mgp=c(1.8,1,0), xpd=NA)
 } # end for loop smooth
dev.off()

rm(dn5,dn6,col5,col6, smooth, dens)



# 2.4. Truncation dependency ---------------------------------------------------
source("Code_aid.R"); aid$load("PREC", "weights"); rm(BGE,weights_all,weights_dn)

trunc <- seq(0,0.95, len=50)
library(parallel) ; library(extremeStat)
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("aid", "PREC", "weights", "trunc")) # 20 min on 3 cores
trqL <- pblapply(trunc, cl=cl, FUN=function(tr) 
  berryFunctions::tryStack(
  extremeStat::distLquantile(log10(PREC), truncate=tr, quiet=TRUE, order=FALSE, 
                             probs=aid$probs, emp=FALSE, weighted=TRUE, weightc=weights),
  file=paste0("simlogs/trunc",formatC(round(tr,4),format="f", digits=4),".txt") ))
stopCluster(cl) ; rm(cl)
trq <- l2array(trqL)
names(dimnames(trq)) <- c("distr","prob","trunc")
dimnames(trq)[[3]] <- formatC(round(trunc,4),format="f", digits=4)
save(trunc,trq, file="dataprods/trunc.Rdata") 
rm(trqL)


# GPD et al:
aid$load("PT")
trunc_stats <- c(0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98) 
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("aid", "PT", "weights", "trunc_stats")) # 30 min on 3 cores
trqL_stats <- pblapply(seq_along(PT), cl=cl, FUN=function(stat)
  berryFunctions::tryStack(
  {
  tqL <- lapply(trunc_stats, FUN=function(tr) extremeStat::distLquantile(log10(PT[[stat]]$prec), 
                 truncate=tr, quiet=TRUE, order=FALSE, probs=aid$probs, weightc=weights))
  tq <- berryFunctions::l2array(tqL)
  names(dimnames(tq)) <- c("distr","prob","trunc")
  dimnames(tq)[[3]] <- trunc_stats
  tq
  }, 
  file=paste0("simlogs/truncstat",stat,".txt") ))
stopCluster(cl) ; rm(cl)
trq_stats <- l2array(trqL_stats)
names(dimnames(trq_stats))[4] <- "stat"
dimnames(trq_stats)[[4]] <- names(PT)
save(trunc_stats,trq_stats, file="dataprods/trunc_stats.Rdata") 
rm(trqL_stats)


# _ Visualisation -----

source("Code_aid.R"); aid$load("PREC","trunc")

pdf("fig/fig6.pdf", height=3.5, width=3.5, pointsize=11)
par(mar=c(3,2.8,0.2,0.4), mgp=c(1.8,0.5,0))
plot(1, type="n", xlab="Truncation proportion", xlim=0:1, xaxs="i", xaxt="n",
     ylab="99.9% quantile estimate  [mm/h]", ylim=c(0.8, 1.9), yaxt="n")
logAxis(2)
axis(1, 0:5*0.2, c("0",1:4*0.2,"1"))
dn <- c("weightedc","wei","gpa")
dnlegend <- c("weighted","Weibull","GPD")
cols <- RColorBrewer::brewer.pal(5, "Set1")[-(1:2)]; names(cols) <- dn
for(d in dimnames(trq)[[1]]) lines(trunc, trq[d,"99.9%",], col=8)
for(d in dn) lines(trunc, trq[d,"99.9%",], col=cols[d], lwd=2)
abline(v=0.8)
abline(h=quantileMean(log10(PREC), probs=0.999), lty=3)
legend("topright", c(dnlegend, "other"),
       col=c(cols,8), lty=1, lwd=c(rep(2,6),1), bg="white", cex=0.8)
text(0.01, log10(15), "empirical quantile (full sample)", cex=0.8, adj=0)
dev.off()
rm(d,dn,dnlegend,cols)


# GPD per station ----
source("Code_aid.R"); aid$load("trunc_stats", "weights"); rm(BGE,weights_all,weights_dn)

dn <- c(names(weights)[1:12], dimnames(trq_stats)[[1]][c(22:34)])
dnlegend <- dn; names(dnlegend) <- dn
dnlegend <- sub("_","-", dnlegend)
dnlegend <- sub("_","\n",dnlegend)
dnlegend <- sub("-","_", dnlegend)

pdf("fig/fig7.pdf", height=5, pointsize=11)
par(mfrow=c(5,5), oma=c(2,2.9,0.2,0.4), mar=c(0,0,0,0), xpd=F, mgp=c(1.8,0.5,0), cex=1)
for(d in dn)
{
plot(1, type="n", xlab="", ylim=c(0.8, 1.9), yaxt="n", ylab="", 
     xlim=0:1, xaxs="i", xaxt="n")
if(d %in% dn[c( 1,21)]) logAxis(2) else logAxis(2, labels=FALSE)
if(d %in% dn[c(21,25)]) axis(1, 0:5*0.2, c("0",1:4*0.2,"1"))
abline(v=0.8, col=8)
for(i in 1:142) lines(trunc_stats, trq_stats[d,"99.9%",,i], col=addAlpha("blue"))
textField(0.5, log10(50), dnlegend[d], cex=0.8)
}
title(xlab="Truncation proportion", outer=TRUE, line=0.7, cex=1)
title(ylab="Quantile estimate  [mm/h]", outer=TRUE, line=1.5, cex=1)
dev.off()

rm(dn,d, dnlegend, i)





# PTQ aggregate across distributions ----

source("Code_aid.R"); aid$load("PTQ")

dn <- dimnames(PTQ)[[1]][1:33]
dc <- rep("grey", length(dn)); names(dc) <- dn
dc[grepl("GPD_",dn)] <- addAlpha("red")
dc[dn=="quantileMean"] <- "blue"
statavs <- pblapply(dn, function(d)
{
  stats <- PTQ[d,"99.9%",,]
  stats <- replace(stats, stats>10, NA)
  statav <- rowMeans(stats, na.rm=TRUE)
  statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
  statav
})
names(statavs) <- dn

pdf("fig/statavs.pdf", height=5)
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(1.5,1.5,0,0) )
aid$PTplot(prob="99.9%", outer=TRUE, line=0, ylim=c(4,120), main="", cc=FALSE)
for(d in dn) lines(as.numeric(names(statavs[[d]])), 10^statavs[[d]], col=dc[d])
aid$cc_lines(NA)
legend("bottomright", c("17 lmomco dists", "11 GPD fits", "empirical", "CC scaling"),
       col=c("grey","red","blue","black"), lty=1, bg="white")
dev.off()
