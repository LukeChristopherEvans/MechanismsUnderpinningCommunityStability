# Mean Abundance model 

# load data
folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
listpath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityButterflyData.csv"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code Structural Equation model/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,listpath))
CommunityData = read.csv(paste0(folderpath,datapath))


set.seed(1234)
meanabNSP<-stan(
  file = paste0(folderpath,stanpath,"logmeanabNonspatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)


set.seed(1234)
meanabSP<-stan(
  file = paste0(folderpath,stanpath,"logmeanabSpatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)

compare(meanabNSP,meanabSP)
traceplot(meanabSP,pars="beta_p",inc_warmup=T)

logmeanpoststamp= extract.samples(meanabSP)

spatialvarplot(logmeanpoststamp)
coefplotter(meanabSP,coefnames =c("Intercept","Niche Mismatch","Mean niche volume","Species \n richness"),type="fixed")
coefplotter(meanabSP,coefnames =c("Intercept","Niche Mismatch","Mean niche volume","Species \n richness"),type="random")


# add for convience
siteniches2$logmeanab = log(siteniches2$meanab)

# interesting line is niche volume
lineequationmab1 <-function(x,i)  logmeanpoststamp$beta_p[,1,i] + logmeanpoststamp$beta_p[,3,i] * x
plotter(logmeanpoststamp,'MeanVol','LgMeanAb','meanvol','logmeanab',lineeq = lineequationmab1,xlab="Mean niche volume",ylab = "Log mean abundance",roundtox = 0,linetype = "dds")

# marginal versions
marginalmab1 <- function(x,y,i) y - ( mean(logmeanpoststamp$beta_p[,1,i]) + mean(logmeanpoststamp$beta_p[,2,i]) * x[,1] )
plotter(logmeanpoststamp,'MeanVol','LgMeanAb','meanvol','logmeanab',lineeq = lineequationmab1,xlab="Mean niche volume",margeq = marginalmab1,ylab = "Log mean abundance",m1="Mismatch",linetype = "dds",yu=3.1,yl=-3.1)

# interesting line is niche volume
lineequationmab2 <-function(x,i)  logmeanpoststamp$beta_p[,1,i] + logmeanpoststamp$beta_p[,2,i] * x
plotter(logmeanpoststamp,'Mismatch','LgMeanAb','meanvol','mismatch',lineeq = lineequationmab2,xlab="Niche mismatch",ylab = "Log mean abundance",roundtox = 0,linetype = "ddd")

# marginal versions
marginalmab2 <- function(x,y,i) y - ( mean(logmeanpoststamp$beta_p[,1,i]) + mean(logmeanpoststamp$beta_p[,3,i]) * x[,1] )
plotter(logmeanpoststamp,'Mismatch','LgMeanAb','meanvol','mismatch',lineeq = lineequationmab2,xlab="Niche mismatch",margeq = marginalmab2,ylab = "Log mean abundance",m1='MeanVol',linetype = "ddd",yu=3.1,yl=-3.1)




# figures for keeping
setwd("")
tiff("lgmeanfixed.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(meanabSP,coefnames =c("Intercept","Niche Mismatch","Mean niche volume","Species \n richness"),type="fixed")
dev.off()

setwd("")
tiff("lgmeanrandom.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(meanabSP,coefnames =c("Intercept","Niche Mismatch","Mean niche volume","Species \n richness"),type="random")
dev.off()

setwd("")
tiff("lgmeanspatialvar.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
spatialvarplot(logmeanpoststamp)
dev.off()

setwd("")
tiff("lgmeanM1meanvolMarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(logmeanpoststamp,'MeanVol','LgMeanAb','meanvol','logmeanab',lineeq = lineequationmab1,xlab="Mean niche volume",margeq = marginalmab1,ylab = "Log mean abundance",m1="Mismatch",linetype = "dds",yu=3.1,yl=-3.1)
dev.off()

setwd("")
tiff("lgmeanM1meanvol_spatial.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(logmeanpoststamp,'MeanVol','LgMeanAb','meanvol','logmeanab',drawcor=T, lineeq = lineequationmab1,xlab="Mean niche volume",ylab = "Log mean abundance",roundtox = 0,linetype = "dds")
dev.off()

setwd("")
tiff("lgmeannichemismatch.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(logmeanpoststamp,'Mismatch','LgMeanAb','meanvol','mismatch',lineeq = lineequationmab2,xlab="Niche mismatch",margeq = marginalmab2,ylab = "Log mean abundance",m1='MeanVol',linetype = "ddd",yu=3.1,yl=-3.1)
dev.off()


summary(meanabSP,pars="beta_p")$summary

# variance partition plots 

## write a function to do it fast across functions 
labelnames = c("Niche mismatch","Mean \n niche volume", "Site","Residuals")
varnames = c("Mismatch","MeanVol")
targetname = c("LgMeanAb")
findplots = FindVarPart(varnames=varnames,targetname = targetname,labelnames = labelnames,postsamp = logmeanpoststamp,ymax=1.52)


setwd("")
tiff("LgMeanVarpar.tif", res=600, compression = "lzw", height=5, width=15, units="in")
plot_grid(findplots[[1]],findplots[[2]] + theme(axis.title.y = element_blank()),findplots[[3]] + theme(axis.title.y = element_blank()) ,ncol=3,labels=c("a","b","c"))
dev.off()

findplots[[4]] #sp
findplots[[5]] #fi
findplots[[6]] #uk


