# synchrony model 

# load data
folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
listpath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityButterflyData.csv"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code Structural Equation model/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,listpath))
CommunityData = read.csv(paste0(folderpath,datapath))


# model structure syn ~ int + jaccard + spprich + spatial intercept
# run the two models 
set.seed(1234)
synNSP=rstan::stan(
  file = paste0(folderpath,stanpath,"synindexNonspatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)
set.seed(1234)
synSP=rstan::stan(
  file = paste0(folderpath,stanpath,"synindexSpatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)
set.seed(1234)
# below compares with and without jaccard included
synSP2=rstan::stan(
  file = paste0(folderpath,stanpath,"synindexSpatial2.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)

compare(synSP,synNSP,synSP2)
traceplot(synSP,pars="beta_p",inc_warmup=T)

synsamp = extract.samples(synSP)

coefplotter(synSP,coefnames =c("Intercept","Jaccard similarity","Species richness"),type="fixed")
coefplotter(synSP,coefnames =c("Intercept","Jaccard similarity","Species richness"),type="random")
spatialvarplot(synsamp)


lineequationmab1 = function(x,i)  synsamp$beta_p[,1,i] + synsamp$beta_p[,2,i] * x
plotter(synsamp,'Jaccard','SynIndex','jaccard','synindex',lineeq = lineequationmab1,xlab="Jaccard similarity",ylab = "Synchrony",linetype = "dds",roundtoy = 2,roundtox = 2)
margeq1 = function(x,y,i) y - (mean(synsamp$beta_p[,1,i]) + mean(synsamp$beta_p[,3,i]) * x[,1])
plotter(synsamp,'Jaccard','SynIndex','jaccard','synindex',lineeq = lineequationmab1,margeq =margeq1 ,xlab="Jaccard similarity",ylab = "Synchrony",m1="SppRich",linetype = "dds",yu=4,yl=-3)


lineequationmab2 = function(x,i)  synsamp$beta_p[,1,i] + synsamp$beta_p[,3,i] * x
plotter(synsamp,'SppRich','SynIndex','spprich','synindex',lineeq = lineequationmab2,xlab="Species richness",ylab = "Synchrony",roundtox=0,linetype="sss")
margeq2 = function(x,y,i) y - (mean(synsamp$beta_p[,1,i]) + mean(synsamp$beta_p[,2,i]) * x[,1])
plotter(synsamp,'SppRich','SynIndex','spprich','synindex',lineeq = lineequationmab2,margeq =margeq2 ,xlab="Species richness",ylab = "Synchrony",m1="Jaccard",roundtox=0,yu=4,yl=-2)

# check other model for differnences - not really any
synsamp2 = extract.samples(synSP2)
lineequationmab3 = function(x,i)  synsamp2$beta_p[,1,i] + synsamp2$beta_p[,2,i] * x
plotter(synsamp2,'SppRich','SynIndex','spprich','synindex',lineeq = lineequationmab3,xlab="Species richness",ylab = "Synchrony",roundtox=0)


# figures for keeping
setwd("")
tiff("synfixed.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(synSP,coefnames =c("Intercept","Jaccard \n similarity","Species \n richness"),type="fixed")
dev.off()

setwd("")
tiff("synrandom.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(synSP,coefnames =c("Intercept","Jaccard \n similarity","Species \n richness"),type="random")
dev.off()

setwd("")
tiff("synspatialvar.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
spatialvarplot(synsamp)
dev.off()

setwd("")
tiff("synjaccardmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(synsamp,'Jaccard','SynIndex','jaccard','synindex',lineeq = lineequationmab1,margeq =margeq1 ,xlab="Jaccard similarity",ylab = "Synchrony",m1="SppRich",linetype = "dds",yu=4,yl=-3,roundtox=2)
dev.off()

setwd("")
tiff("synspprichmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(synsamp,'SppRich','SynIndex','spprich','synindex',lineeq = lineequationmab2,margeq =margeq2 ,xlab="Species richness",ylab = "Synchrony",m1="Jaccard",roundtox=0,yu=4,yl=-2)
dev.off()
#I'm not going to do the non-marginals because they are misleading - and the spatial field doesn't do much

# details for a table should they be required 
summary(synSP,pars="beta_p")$summary


## write a function to do it fast across functions 
labelnames = c("Jaccard \n similarity","Species \n richness","Site","Residuals")
varnames = c("Jaccard","SppRich")
targetname = c("SynIndex")
findplots = FindVarPart(varnames=varnames,targetname = targetname,labelnames = labelnames,postsamp = synsamp,ymax=1.4,box=T)


setwd("")
tiff("synVarpar.tif", res=600, compression = "lzw", height=5, width=15, units="in")
plot_grid(findplots[[1]],findplots[[2]] + theme(axis.title.y = element_blank()),findplots[[3]] + theme(axis.title.y = element_blank()) ,ncol=3,labels=c("a","b","c"))
dev.off()

findplots[[4]] #sp
findplots[[5]] #fi
findplots[[6]] #uk





