### R script for hierbaps ###

library(ape)
library(BactDating)
library(coda)

path <- "/home/user/Documents/GBS_typeIV"

tGBSVtemp <- loadCFML(paste(path,"FOLDER",paste("CFML","FOLDER",sep=""),sep="/"))

## export tip.label for dating order
write.csv(tGBSVtemp$tip.label,file=paste(path,"FOLDER",paste("CFML","FOLDER","-dates.csv",sep=""),sep="/"))

## use csv file from above to generate list of dates in tip.label order for input as "d"
dGBSVtemp=c(read.csv(file=paste(path,"FOLDER",paste("FOLDER","-dates.csv",sep=""),sep="/")))

## to reroot tree to best root
rootedGBSVtemp=initRoot(tGBSVtemp,dGBSVtemp$d)

## assess strength of temporal signal
rGBSVtemp=roottotip(rootedGBSVtemp,dGBSVtemp$d)

## assuming strong clock-like behavior main analysis; note this analysis uses lenght of 10000 for MCMC; will need to verify in subsequent steps
resGBSVtemp=bactdate(tGBSVtemp,dGBSVtemp$d,showProgress = T)

## check for MCMC convergence - effective sample size of the parameters of which should be at least 100 for each parameter
mcmcGBSVdefaulttemp=as.mcmc(resGBSVtemp)
effectiveSize(mcmcGBSVdefaulttemp)
## also check traces
plot(resGBSVtemp,'trace')

GBSVres100ktemp=bactdate(tGBSVtemp,dGBSVtemp$d,nbIts=100000,showProgress = T)
plot(GBSVres100ktemp,'trace')

modelcompare(resGBSVtemp,GBSVres100ktemp)

GBSVres1Mtemp=bactdate(tGBSVtemp,dGBSVtemp$d,nbIts=1000000,showProgress = T)
plot(GBSVres1Mtemp,'trace')

modelcompare(GBSVres100ktemp,GBSVres1Mtemp)
