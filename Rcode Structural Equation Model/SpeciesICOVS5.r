# species ICOV model 

# load data
folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
listpath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityButterflyData.csv"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code Structural Equation model/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,listpath))
CommunityData = read.csv(paste0(folderpath,datapath))


# model structure sICOV ~ int + mismatch + meanvol + logmeanab + spatial intercept
# run the two models 
set.seed(1234)
sicovSP<-stan(
  file = paste0(folderpath,stanpath,"SpeciesIcovSpatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)
set.seed(1234)
sicovNSP<-stan(
  file = paste0(folderpath,stanpath,"SpeciesIcovNonspatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.9999)
)

compare(sicovSP,sicovNSP)
traceplot(sicovSP,pars="beta_p",inc_warmup=T)
sisamp <- extract.samples(sicovSP)

# add logmean
siteniches2$logmeanab = log(siteniches2$meanab)

coefplotter(sicovSP,coefnames =c("Intercept","Niche mismatch","Mean niche volume","Log mean abundance"),type="fixed")
coefplotter(sicovSP,coefnames =c("Intercept","Niche mismatch","Mean niche volume","Log mean abundance"),type="random")
spatialvarplot(sisamp)



lineequationmab1 <- function(x,i)  sisamp$beta_p[,1,i] + sisamp$beta_p[,2,i] * x
plotter(sisamp,'Mismatch','SpeICOV','mismatch','speciesICOV',lineeq = lineequationmab1,xlab="Niche mismatch",ylab = "Average species stability")
margeq1 <- function(x,y,i) y - (mean(sisamp$beta_p[,1,i]) + mean(sisamp$beta_p[,3,i]) * x[,1] + mean(sisamp$beta_p[,4,i]) * x[,2])
plotter(sisamp,'Mismatch','SpeICOV','mismatch','speciesICOV',lineeq = lineequationmab1,margeq =margeq1,xlab="Niche mismatch",ylab = "Average species stability",m1="MeanVol",m2="LgMeanAb",yu=5,yl=-2)

lineequationmab2 <- function(x,i)  sisamp$beta_p[,1,i] + sisamp$beta_p[,4,i] * x
plotter(sisamp,'LgMeanAb','SpeICOV','logmeanab','speciesICOV',lineeq = lineequationmab2,xlab="Log mean abundance",ylab = "Average species stability")
margeq2 <- function(x,y,i) y - (mean(sisamp$beta_p[,1,i]) + mean(sisamp$beta_p[,2,i]) * x[,1] + mean(sisamp$beta_p[,3,i]) * x[,2])
plotter(sisamp,'LgMeanAb','SpeICOV','logmeanab','speciesICOV',lineeq = lineequationmab2,margeq =margeq2,xlab="Log mean abundance",ylab = "Average species stability",m1="Mismatch",m2="MeanVol",yu=6,yl=-2)



# figures for keeping
setwd("")
tiff("specicovfixed.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(sicovSP,coefnames =c("Intercept","Niche \n mismatch","Mean niche \n volume","Log \n mean abundance"),type="fixed")
dev.off()

setwd("")
tiff("specicovrandom.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(sicovSP,coefnames =c("Intercept","Niche mismatch","Mean niche \n volume","Log mean \n abundance"),type="random")
dev.off()

setwd("")
tiff("specicovspatialvar.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
spatialvarplot(sisamp)
dev.off()

setwd("")
tiff("specicovMisMmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(sisamp,'Mismatch','SpeICOV','mismatch','speciesICOV',lineeq = lineequationmab1,margeq =margeq1,xlab="Niche mismatch",ylab = "Average species stability",m1="MeanVol",m2="LgMeanAb",yu=5,yl=-2)
dev.off()

setwd("")
tiff("specicovLgMAbmarg.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(sisamp,'LgMeanAb','SpeICOV','logmeanab','speciesICOV',lineeq = lineequationmab2,margeq =margeq2,xlab="Log mean abundance",ylab = "Average species stability",m1="Mismatch",m2="MeanVol",yu=6,yl=-2)
dev.off()
#I'm not going to do the non-marginals because they are misleading - and the spatial field doesn't do much

# details for a table should they be required 
summary(sicovSP,pars="beta_p")$summary

# variance partition plots 

## write a function to do it fast across functions 
labelnames = c("Niche \n mismatch","Mean \n niche \n volume","Log mean \n abundance" ,"Site","Residuals")
varnames = c("Mismatch","MeanVol","LgMeanAb")
targetname = c("SpeICOV")
findplots = FindVarPart(varnames=varnames,targetname = targetname,labelnames = labelnames,postsamp = sisamp,ymax=1.73)


setwd("/")
tiff("specicovVarpar.tif", res=600, compression = "lzw", height=5, width=15, units="in")
plot_grid(findplots[[1]],findplots[[2]] + theme(axis.title.y = element_blank()),findplots[[3]] + theme(axis.title.y = element_blank()) ,ncol=3,labels=c("a","b","c"))
dev.off()

findplots[[4]] #sp
findplots[[5]] #fi
findplots[[6]] #uk

