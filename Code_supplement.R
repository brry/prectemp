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


