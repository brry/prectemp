# Analysis of temperature-dependency of rainfall intensity
# supplementary paper
# Berry Boessenkool, 2017
browseURL("http://www.nat-hazards-earth-syst-sci-discuss.net/nhess-2016-183")



# 1. station map ----

# CHECK: coordinates #   View(meta) # Map of stations
load("dataprods/meta.Rdata") ; library(OSMscale)
map <- pointsMap(lat,long, data=meta, zoom=6, type="maptoolkit-topo")
meta <- cbind(meta, 
              projectPoints(lat,long, data=meta, from=pll(), to=pmap(map)) )

pdf("fig/sup1_stationsmap.pdf", width=5)
par(mar=rep(0,4))
plot(map, removeMargin=FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col=addAlpha("white",0.5), border=NA)
#colPoints(x, y, m_dur, data=meta, zlab="Time series length  [Years]", 
#          Range=c(14.9,20.04), legargs=list(mar=c(0.8,0.5,1,0.5), lines=F),
#          density=FALSE, x1=0.5, y1=0.93, y2=1,x2=1, cex=2)
colPoints(x, y, ele, data=meta, zlab="Station elevation  [m asl]", 
          legargs=list(mar=c(0.8,0.5,1,0.5), lines=F), #col=terrain.colors(110)[1:100],
          density=FALSE, x1=0.5, y1=0.93, y2=1,x2=1, cex=2)
points(y~x, data=meta, cex=2)
dev.off()



# 2. ssd GPD fitting methods ----

source("Code_aid.R"); aid$load("simQA", "PREC")

dn_mle <- c("GPD_MLE_evd", "GPD_MLE_evir", "GPD_MLE_extRemes", "GPD_MLE_fExtremes", 
            "GPD_MLE_ismev", "GPD_MLE_Renext_Renouv")
dn_mom <- c("GPD_LMO_extRemes", "GPD_GML_extRemes", "GPD_PWM_fExtremes", 
            "GPD_LMO_lmomco", "GPD_PWM_evir")
col_mle <- RColorBrewer::brewer.pal(6, "Set2") ; names(col_mle) <- dn_mle
col_mom <- RColorBrewer::brewer.pal(5, "Set2") ; names(col_mom) <- dn_mom

pdf("fig/sup2_SSD.pdf", height=3, width=3.5, pointsize=10)
#for(smooth in c(1,3,5,7,9,11,13,15)){
smooth <- 9; {
par(mfrow=c(1,2), mar=c(2,0,0.2,0.4), oma=c(1,3,0,0), mgp=c(1.8,0.7,0), las=1)
## panel 1
plot(1, type="n", xlim=c(50,830), ylim=log10(c(7,25)), xaxs="i", main="", yaxt="n",
       xlab="", ylab="")
logAxis(2)
for(d in dn_mle) lines(aid$n, movAv(simQA["50%",d,"99.9%",], smooth), col=col_mle[d])
abline(h=quantileMean(log10(PREC), probs=0.999), lty=3)
text(90, log10(24), "MLE", adj=0)
legend("bottomright", replace(dn_mle, 6, "GPD_MLE_Renext"), lwd=2, col=col_mom, 
       bg="white", cex=0.65)
## panel 2  #par(mar=c(3,0,0.2,0.4))
plot(1, type="n", xlim=c(50,830), ylim=log10(c(7,25)), xaxs="i", main="", yaxt="n",
       xlab="", ylab="")
logAxis(2, labels=FALSE)
for(d in dn_mom) lines(aid$n, movAv(simQA["50%",d,"99.9%",], smooth), col=col_mom[d])
abline(h=quantileMean(log10(PREC), probs=0.999), lty=3)
text(90, log10(24), "LM / PWM", adj=0)
legend("bottomright", dn_mom, lwd=2, col=col_mom, bg="white", cex=0.65)
#
#dens <- density(log10(PREC)); lines(25+dens$y*1000, dens$x)
title(xlab="sample size n", outer=TRUE, mgp=c(-0.2,-1,0), xpd=NA)
title(ylab="Random sample 99.9% quantile  [mm/h]", outer=TRUE, mgp=c(1.8,1,0), xpd=NA)
 } # end for loop smooth
dev.off()

rm(dn_mle,dn_mom,col_mle,col_mom, smooth)



# 3. Truncation dependency ----
# 3.1. truncdep computation ----
source("Code_aid.R"); aid$load("PT","PREC"); library(parallel) ; library(extremeStat)

# truncation dependency full dataset
truncs <- c(0.7,0.75,0.8, 0.85,seq(0.86,0.995, 0.005))
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("PREC", "truncs"))
gpds <- pbsapply(truncs, cl=cl, FUN=function(t)  # 1 minute for PREC_11, 3 for PREC_all
  extremeStat::distLquantile(log10(PREC), probs=0.999, truncate=t, sel="gpa", emp=F, quiet=T)[1,1])
stopCluster(cl) ; rm(cl)

pdf("fig/gpa_truncdep.pdf", height=5)
par(mar=c(3.5,3.8,1,0.5), mgp=c(2.1, 0.7,0))
plot(truncs, 10^gpds, type="n", axes=T, las=1, ylab="", xlab="Truncation proportion")
title(ylab="GPD quantile estimate of complete pooled dataset", mgp=c(2.8,1,0))
abline(h=10^quantile(log10(PREC), 0.999, type=8), col=8)
lines(truncs, 10^gpds, lwd=2, type="o")
dev.off()

truncs <- c(0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98) 
cl <- makeCluster( detectCores()-1 )
clusterExport(cl, c("aid", "PT", "truncs")) # 15 min on 7 cores
trqL_stats <- pblapply(seq_along(PT), cl=cl, FUN=function(stat)
  berryFunctions::tryStack(
  {
  tq <- lapply(truncs, FUN=function(tr) extremeStat::distLquantile(
               log10(PT[[stat]]$prec), truncate=tr, quiet=TRUE, probs=aid$probs, 
               sel=c("wei","pe3","wak","gno","gum"), order=FALSE, 
               emp=FALSE, weighted=FALSE, gpd=TRUE))
  tq <- berryFunctions::l2array(tq)
  names(dimnames(tq)) <- c("distr","prob","trunc")
  dimnames(tq)[[3]] <- truncs
  tq
  }, 
  file=paste0("sim_logs/truncstat",stat,".txt") ))
stopCluster(cl) ; rm(cl)
trunq <- l2array(trqL_stats)
names(dimnames(trunq))[4] <- "stat"
dimnames(trunq)[[4]] <- names(PT)
save(truncs,trunq, file="dataprods/trunc.Rdata") 
rm(trqL_stats)


# 3.2. truncdep visualisation ----
# Visualisation truncdep per station (mainly GPD comparison)
source("Code_aid.R"); aid$load("trunc")

dn <- c("GPD_MLE_evd", "GPD_MLE_extRemes", "GPD_MLE_ismev", "GPD_MLE_Renext_Renouv",
        "GPD_MLE_evir", "GPD_MLE_fExtremes", "dummy", "GPD_GML_extRemes",
        "GPD_PWM_fExtremes",  "GPD_PWM_evir", "GPD_LMO_extRemes", "GPD_LMO_lmomco")
bg <- c(rep("seashell",6), rep("lavender",2), rep("mistyrose",4))
names(bg) <- dn
dnlegend <- dn; names(dnlegend) <- dn
# Replace second "_" with newline:
dnlegend <- sub("_","-", dnlegend)
dnlegend <- sub("_","\n",dnlegend)
dnlegend <- sub("-","_", dnlegend)


pdf("fig/sup3_truncdep.pdf", height=5, pointsize=10)
par(mfrow=c(3,4), oma=c(2,2.9,0.2,0.4), mar=c(0,0,0,0), xpd=F, mgp=c(1.8,0.5,0), cex=1)
for(d in dn)
{
#par(bg=unname(bg[d]), mar=rep(1,4))
plot(1, type="n", xlab="", ylim=c(0.8, 1.9), yaxt="n", ylab="", 
     xlim=0:1, xaxs="i", xaxt="n")
if(d=="dummy") next
u <- par("usr") 
rect(u[1], u[3], u[2], u[4], col=bg[d], border=NA) 
rm(u)
if(d %in% dn[c( 1,9)]) logAxis(2) else logAxis(2, labels=FALSE)
if(d %in% dn[c(9,12)]) axis(1, 0:5*0.2, c("0",1:4*0.2,"1"))
abline(v=0.9, col=8)
for(i in 1:142) lines(truncs, trunq[d,"99.9%",,i], col=addAlpha("blue"))
textField(0.5, log10(50), dnlegend[d], cex=0.8, fill=bg[d])
}
title(xlab="Truncation proportion", outer=TRUE, line=0.7, cex=1)
title(ylab="Quantile estimate  [mm/h]", outer=TRUE, line=1.5, cex=1)
dev.off()

rm(dn,d, dnlegend, i, bg)



# 4. PrecTemp GPD ----

source("Code_aid.R"); aid$load("PTQ")

dn <- c("GPD_GML_extRemes", "GPD_LMO_extRemes", "GPD_LMO_lmomco", 
        "GPD_PWM_fExtremes", "GPD_PWM_evir", 
        "GPD_MLE_fExtremes", "GPD_MLE_Renext_Renouv", "GPD_MLE_evir", "GPD_MLE_evd",
        "GPD_MLE_extRemes", "GPD_MLE_ismev", "empirical")

dc <- rep(addAlpha("blue"), length(dn))
names(dc) <- dn
dc["GPD_LMO_lmomco"] <- "blue"
dc[grepl("GPD_MLE",dn)] <- addAlpha("red")
dc["empirical"] <- "green3"
statavs <- pblapply(dn, function(d)
{
  stats <- PTQ[d,"99.9%",,]
  stats <- replace(stats, stats>10, NA)
  statav <- rowMeans(stats, na.rm=TRUE)
  statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
  statav
})
names(statavs) <- dn

pdf("fig/sup4_prectemp_GPD.pdf", height=5)
#for(dist in dn){
par(mfrow=c(1,1), mar=c(2,2,0.5,0.5), oma=c(1.5,1.5,0,0) )
aid$PTplot(prob="99.9%", outer=TRUE, line=0.2, xlim=c(4.9,19.8), ylim=c(6,60), 
           xaxs="i", main="", cc=FALSE)
for(d in dn) lines(as.numeric(names(statavs[[d]])), 10^statavs[[d]], col=dc[d], lwd=2)
aid$cc_lines(NA, anf=4)
#lines(as.numeric(names(statavs[[dist]])), 10^statavs[[dist]], col=2, lwd=2)
#title(main=dist, line=-2)
#}
legend("bottomright", c("LM/PWM", "MLE", "Empirical", "CC scaling"),
       col=c("blue","red","green3","black"), lwd=2, bg="white")
dev.off()




# ......... ----
# 5. Independency check ----
load("dataprods/PT.Rdata")

# 5.1. Time between events ----

# threshold exceedances (at 90% truncation)
PTt <- lapply(PT, function(x) x[x$prec>quantile(x$prec,0.9),])

diffs <- sapply(PTt, function(x){d <- as.numeric(diff(x$date))
                                 c(d1=mean(d==1), d10=mean(d<=10), d24=mean(d<=24))} )
diffs <- as.data.frame(t(diffs*100))

png("fig/timediff.png", width=5, height=7, units="in", res=500)
par(mfrow=c(3,1), mar=c(3,3,2,1), oma=c(0,0,4,0), mgp=c(2,0.7,0), las=1)
hist(diffs$d1 , breaks=50, main=round(median(diffs$d1) ), col="moccasin")
hist(diffs$d10, breaks=50, main=round(median(diffs$d10)), col="moccasin")
hist(diffs$d24, breaks=50, main=round(median(diffs$d24)), col="moccasin")
title(main="Time difference between threshold exceedances
      Median and histogram of percentage 
      of differences below 1, 10 and 24 hours", outer=TRUE)
dev.off()


# 5.2. Pq-T for independent events ----
source("Code_aid.R"); aid$load("PT"); library(extremeStat)

# https://stackoverflow.com/questions/45440342

eventize <- function(x, mindiff=1) 
  {
  diffs <- as.difftime(diff(x$date))
  units(diffs) <- "hours"
  diffs <- as.numeric(diffs) 
  diffs <- c(0, diffs)
  event <- cumsum(diffs>mindiff) # distinct event if more than n hours apart
  xmax <- unlist(tapply(x$prec, event, FUN=function(v){
    out <- rep(0, length(v))
    out[which.max(v)] <- 1 # select first maximum value if there are ties
    out
    }))
  x[xmax==1,]
  }

PTe01 <- pbsapply(PT, eventize,             simplify=FALSE)
PTe05 <- pbsapply(PT, eventize, mindiff= 5, simplify=FALSE)
PTe12 <- pbsapply(PT, eventize, mindiff=12, simplify=FALSE)
rm(eventize)

PTQ <- function(i, obj)
  {
  x <- obj[[i]] # x <- PTe01[[3]]; t=15.5
  # Quantile estimates per temperature bin
  binQ <- lapply(aid$mid, function(t) {
    seldat <- x$prec[ x$temp5>(t-1) & x$temp5<=(t+1)]
    extremeStat::distLquantile(log10(seldat), probs=aid$probs, truncate=0.9,  
                               sel="gpa", gpd=FALSE, order=FALSE, quiet=TRUE, weighted=FALSE)
    })
  binQ2 <- berryFunctions::l2array(binQ)
  names(dimnames(binQ2)) <- c("distr","prob","temp") ; dimnames(binQ2)[[3]] <- aid$mid
  return(binQ2)
}
PTQe01 <- pblapply(X=1:142, FUN=PTQ, obj=PTe01) # 4:33 min
PTQe05 <- pblapply(X=1:142, FUN=PTQ, obj=PTe05) # 3:49 min
PTQe12 <- pblapply(X=1:142, FUN=PTQ, obj=PTe12) # 3:29 min
rm(PTQ)
l2a <- function(x) {x <- l2array(x); 
                    names(dimnames(x))[4] <- "stat"
                    dimnames(x)[[4]] <- names(PT)
                    x}
PTQe01 <- l2a(PTQe01)
PTQe05 <- l2a(PTQe05)
PTQe12 <- l2a(PTQe12)
rm(l2a)
save(PTQe01,PTQe05,PTQe12, file="dataprods/PTQe.Rdata") 


# 5.3. Pq-T events Vis ----
source("Code_aid.R"); aid$load("PTQe", "PTQ")

PTQlines <- function(
PTQ,
prob="99.9%",  
dn="gpa",
cut=150,
col=addAlpha("blue", 0.15),
empav=NULL,
...
)
{
aid$PTplot(xlim=c(4.8,20.3), ylim=c(5,150), cc=FALSE)
for(i in 1:142) lines(aid$mid, 10^PTQ[dn,prob,,i], col=col, ...) 
stats <- PTQ[dn,prob,,]
stats <- replace(stats, stats>cut, NA)  
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
if(!is.null(empav)) lines(aid$mid, 10^empav, col="green3", lwd=3) 
lines(aid$mid, 10^statav, lwd=3)
aid$cc_lines(NA, mainargs=list(col=2, lty=2))
return(invisible(statav))
}


pdf("fig/sup5_PTQ_independent_events.pdf", height=5)
par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), oma=c(2,2,0,0))
statave <- PTQlines(PTQ, dn="empirical", col=addAlpha("green3", 0.15))
legend("topleft", "Empirical, all rainfall observations", bty="n")
PTQlines(PTQe01, empav=statave); legend("topleft", "GPD, event maxima (>=1 hour apart)", bty="n")
PTQlines(PTQe05, empav=statave); legend("topleft", "GPD (>=5 hours apart)", bty="n")
PTQlines(PTQe12, empav=statave); legend("topleft", "GPD (>=12 hours apart)", bty="n")
title(xlab="Dewpoint temperature (mean of preceding 5 hours)  [ \U{00B0}C]", outer=TRUE, line=0.5)
title(ylab="Event maximum precipitation 99.9% quantile  [mm/h]", outer=TRUE, line=0.5)
dev.off()




# 6. PTQ all quantiles ----

source("Code_aid.R"); aid$load("PTQ")

PTQlines <- function(
prob="",  
dn="gpa",
cut=150,
col=addAlpha("blue", 0.15),
...
)
{
aid$PTplot(xlim=c(4.8,20.3), ylim=c(5,150), cc=FALSE)
for(i in 1:142) lines(aid$mid, 10^PTQ[dn,prob,,i], col=col, ...) 
stats <- PTQ[dn,prob,,]
stats <- replace(stats, stats>cut, NA)  
statav <- rowMeans(stats, na.rm=TRUE)
statav[apply(stats,1,function(x) sum(!is.na(x))<50)] <- NA
#
estats <- PTQ["empirical",prob,,]
estats <- replace(estats, estats>cut, NA)  
estatav <- rowMeans(estats, na.rm=TRUE)
estatav[apply(estats,1,function(x) sum(!is.na(x))<50)] <- NA
lines(aid$mid, 10^estatav, lwd=3, col="green3")
#
lines(aid$mid, 10^statav, lwd=3)
aid$cc_lines(NA, mainargs=list(col=2, lty=2))
legend("topleft", prob, bty="n")
return(invisible(statav))
}


pdf("fig/sup6_PT_allprobs.pdf", height=5)
par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), oma=c(2,2,0,0))
for(prob in names(aid$probcols) ) PTQlines(prob=prob)
title(xlab="Dewpoint temperature (mean of preceding 5 hours)  [ \U{00B0}C]", outer=TRUE, line=0.5)
title(ylab="Precipitation quantile  [mm/h]", outer=TRUE, line=0.5)
rm(prob)
dev.off()


pdf("fig/PT_99.pdf", height=5)
par(mar=c(4,4,0.5,0.5))
aid$PTplot(prob="99%", xlim=c(7,21), ylim=c(7,60), main="", line=2.5)
for(i in 1:142) lines(aid$mid, 10^PTQ["empirical","99%",,i], col=addAlpha("green3", 0.3)) 
for(i in 1:142) lines(aid$mid, 10^PTQ["gpa",      "99%",,i], col=addAlpha("blue"  , 0.3))
rm(i)
dev.off()


# 7. n per bin ----

source("Code_aid.R"); aid$load("PTQ")
statav <- rowMeans(PTQ["n_full","90%",,], na.rm=TRUE)

pdf("fig/sup7_nperbin.pdf", height=4)
par(mfrow=c(1,2), mar=c(3.4,4.5,0.5,0.5))
# linear axis:
plot(1, xlim=c(4.8,20.5), xaxs="i", ylim=lim0(2100), type="n", las=1, ann=F)
title(xlab="Dewpoint temperature (mean of preceding 5 hours)  [ \U{00B0}C]", 
      line=-1.3, outer=TRUE)
title(ylab="sample size per temperature bin", line=3.3)
for(i in 1:142) lines(aid$mid, PTQ["n_full","90%",,i], col=addAlpha(1) ) ; rm(i)
lines(aid$mid, statav, lwd=3, col="red")
# log scale:
aid$PTplot(xlim=c(4.8,20.5), ylim=c(50,2000), main="", cc=FALSE, xlab="",
           ylab="n per bin, log scale", line=3.3)
for(i in 1:142) lines(aid$mid, PTQ["n_full","90%",,i], col=addAlpha(1) ) ; rm(i)
lines(aid$mid, statav, lwd=3, col="red")
dev.off()

