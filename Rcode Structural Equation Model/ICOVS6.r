# ICOV model 

# load data
folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code Structural Equation model/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,datapath))


set.seed(1234)
icovNSP<-stan(
  file = paste0(folderpath,stanpath,"ICOVNonspatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15)
)
set.seed(1234)
icovSP<-stan(
  file = paste0(folderpath,stanpath,"ICOVSpatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15)
)

set.seed(1234)
sitenichel$SynIndex2 = sitenichel$SynIndex^2
icovSP2<-stan(
  file = paste0(folderpath,stanpath,"ICOVSpatial2.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15)
)

compare(icovNSP,icovSP,icovSP2) 
compare(icovNSP,icovSP) 
traceplot(icovSP,pars="beta_p",inc_warmup=T)

icovsamp <- extract.samples(icovSP)
icovsamp2 <- extract.samples(icovSP2)

siteniches2$logmeancomm = log(siteniches2$meancommsize)

# compare two different forms
coefplotter(icovSP,coefnames =c("Intercept","Average species stability","Synchrony","Log mean community abundance","Species richness"),type="fixed")
coefplotter(icovSP2,coefnames =c("Intercept","Average species stability","Synchrony","Log mean community abundance","Synchrony2"),type="fixed")

coefplotter(icovSP,coefnames =c("Intercept","Average species stability","Synchrony","Log mean community abundance","Species richness"),type="random")
coefplotter(icovSP2,coefnames =c("Intercept","Average species stability","Synchrony","Log mean community abundance","Synchrony2"),type="random")

spatialvarplot(icovsamp)


lineequationmab1 <- function(x,i)  icovsamp$beta_p[,1,i] + icovsamp$beta_p[,2,i] * x 
plotter(icovsamp,'SpeICOV','ICOV','speciesICOV','ICOV',lineeq = lineequationmab1,xlab="Average population stability",ylab = "Community stability")
margeq1 <- function(x,y,i) y - (mean(icovsamp$beta_p[,1,i]) + mean(icovsamp$beta_p[,3,i]) * x[,1] + mean(icovsamp$beta_p[,4,i]) * x[,2]+ mean(icovsamp$beta_p[,5,i]) * x[,3])
plotter(icovsamp,'SpeICOV','ICOV','speciesICOV','ICOV',lineeq = lineequationmab1,margeq =margeq1,xlab="Average population stability",ylab = "Community stability",m1="SynIndex",m2="LgMeanCom",m3="SppRich",yl=-2.3,yu=2.8,linetype = "sds")

lineequationmab2 <- function(x,i)  icovsamp$beta_p[,1,i] + icovsamp$beta_p[,3,i] * x  
plotter(icovsamp,'SynIndex','ICOV','synindex','ICOV',lineeq = lineequationmab2,xlab="Synchrony",ylab =  "Community stability")
margeq2 <- function(x,y,i) y - (mean(icovsamp$beta_p[,1,i]) + mean(icovsamp$beta_p[,2,i]) * x[,1] + mean(icovsamp$beta_p[,4,i]) * x[,2]+ mean(icovsamp$beta_p[,5,i]) * x[,3])
plotter(icovsamp,'SynIndex','ICOV','synindex','ICOV',lineeq = lineequationmab2,margeq =margeq2,xlab="Synchrony",ylab =  "Community stability",m1="SpeICOV",m2="LgMeanCom",m3="SppRich",yl=-2.8,yu=2.8,roundtox = 2)

lineequationmab3 <- function(x,i)  icovsamp$beta_p[,1,i] + icovsamp$beta_p[,4,i] * x
plotter(icovsamp,'LgMeanCom','ICOV','logmeancomm','ICOV',lineeq = lineequationmab3,xlab="Log mean community abundance",ylab = "Community stability")
margeq3 <- function(x,y,i) y - (mean(icovsamp$beta_p[,1,i]) + mean(icovsamp$beta_p[,2,i]) * x[,1] + mean(icovsamp$beta_p[,3,i]) * x[,2] + mean(icovsamp$beta_p[,5,i])* x[,3])
plotter(icovsamp,'LgMeanCom','ICOV','logmeancomm','ICOV',lineeq = lineequationmab3,margeq =margeq3,xlab="Log mean community abundance",ylab = "Community stability",m1="SpeICOV",m2="SynIndex",m3="SppRich",yl=-3,yu=3,linetype = "dss")


lineequationmab4 <- function(x,i)  icovsamp$beta_p[,1,i] + icovsamp$beta_p[,5,i] * x
plotter(icovsamp,'SppRich','ICOV','spprich','ICOV',lineeq = lineequationmab3,xlab="Species richness",ylab = "Community stability")
margeq4 <- function(x,y,i) y - (mean(icovsamp$beta_p[,1,i]) + mean(icovsamp$beta_p[,2,i]) * x[,1] + mean(icovsamp$beta_p[,3,i]) * x[,2] + mean(icovsamp$beta_p[,4,i])* x[,3])
plotter(icovsamp,'SppRich','ICOV','spprich','ICOV',lineeq = lineequationmab4,margeq =margeq4,xlab="Species richness",ylab = "Community stability",m1="SpeICOV",m2="SynIndex",m3="LgMeanCom",yl=-3,yu=3,linetype = "ddd",roundtox = 0)


# figures for keeping
setwd("")
tiff("icovfixed.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(icovSP,coefnames =c("Intercept","Average population \n stability","Synchrony","Log mean \n community abundance","Species \n richness"),type="fixed")
dev.off()

setwd("")
tiff("icovrandom.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(icovSP,coefnames =c("Intercept","Average population \n stability","Synchrony","Log mean \n community abundance","Species \n richness"),type="random")
dev.off()

setwd("s")
tiff("icovspatialvar.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
spatialvarplot(icovsamp)
dev.off()

setwd("")
tiff("icovSicovMmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(icovsamp,'SpeICOV','ICOV','speciesICOV','ICOV',lineeq = lineequationmab1,margeq =margeq1,xlab="Average population stability",ylab = "Community stability",m1="SynIndex",m2="LgMeanCom",m3="SppRich",yl=-2.3,yu=2.8,linetype = "sds")
dev.off()

setwd("")
tiff("icovSynMmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(icovsamp,'SynIndex','ICOV','synindex','ICOV',lineeq = lineequationmab2,margeq =margeq2,xlab="Synchrony",ylab =  "Community stability",m1="SpeICOV",m2="LgMeanCom",m3="SppRich",yl=-2.8,yu=2.8,roundtox = 2)
dev.off()

setwd("")
tiff("icovLgMCAmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(icovsamp,'LgMeanCom','ICOV','logmeancomm','ICOV',lineeq = lineequationmab3,margeq =margeq3,xlab="Log mean community abundance",ylab = "Community stability",m1="SpeICOV",m2="SynIndex",m3="SppRich",yl=-3,yu=3,linetype = "dss")
dev.off()

setwd("")
tiff("icovSpRichmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(icovsamp,'SppRich','ICOV','spprich','ICOV',lineeq = lineequationmab4,margeq =margeq4,xlab="Species richness",ylab = "Community stability",m1="SpeICOV",m2="SynIndex",m3="LgMeanCom",yl=-3,yu=3,linetype = "ddd",roundtox = 0)
dev.off()



# details for a table should they be required 
summary(icovSP,pars="beta_p")$summary


# variance partition plots 
## write a function to do it fast across functions 
labelnames = c("Average \n population \n stability","Synchrony","Log mean \n community \n abundance","Species \n richness", "Site","Residuals")
varnames = c("SpeICOV","SynIndex","LgMeanCom","SppRich")
targetname = c("ICOV")
findplots = FindVarPart(varnames=varnames,targetname = targetname,labelnames = labelnames,postsamp = icovsamp,ymax=1.29,box=T)


setwd("")
tiff("ICOvVarpar.tif", res=600, compression = "lzw", height=5, width=17, units="in")
plot_grid(findplots[[1]],findplots[[2]] + theme(axis.title.y = element_blank()),findplots[[3]] + theme(axis.title.y = element_blank()) ,ncol=3,labels=c("a","b","c"))
dev.off()

findplots[[4]] #sp
findplots[[5]] #fi
findplots[[6]] #uk
