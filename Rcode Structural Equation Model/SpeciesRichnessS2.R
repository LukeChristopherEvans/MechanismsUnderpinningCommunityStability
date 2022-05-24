# Species richness model 

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
spprichNSP<-stan(
  file = paste0(folderpath,stanpath,"spprichstanNonspatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)
set.seed(1234)
spprichSP<-stan(
  file = paste0(folderpath,stanpath,"spprichstanSpatial.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99,max_treedepth = 15)
)
# 

compare(spprichNSP,spprichSP)
traceplot(spprichSP,pars="beta_p",inc_warmup=T)

posteriorsamp =  extract.samples(spprichSP)

coefplotter(spprichSP,coefnames =c("Intercept","Overall Niche distance"),type="fixed")
coefplotter(spprichSP,coefnames =c("Intercept","Overall Niche distance"),type="random")

# tests 
lineequation <- function(x,i)  posteriorsamp$beta_p[,1,i] + posteriorsamp$beta_p[,2,i] * x
plotter(posteriorsamp,'OverallNiDist','SppRich','overallnichedist','spprich',lineeq = lineequation,drawcor = F,roundtoy=0,roundtox = 0,xlab="Niche distance",ylab = "Species richness",linetype = "sdd")

# figures for keeping
setwd("")
tiff("spprichfixed.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(spprichSP,coefnames =c("Intercept","Overall Niche distance"),type="fixed")
dev.off()

setwd("")
tiff("spprichrandom.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
coefplotter(spprichSP,coefnames =c("Intercept","Overall Niche distance"),type="random")
dev.off()

setwd("")
tiff("spprichspatialvar.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
spatialvarplot(posteriorsamp)
dev.off()

setwd("")
tiff("spprichm1.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(posteriorsamp,'OverallNiDist','SppRich','overallnichedist','spprich',lineeq = lineequation,drawcor = F,roundtoy=0,xlab="Niche distance",ylab = "Species richness",linetype = "sdd")
dev.off()

setwd("")
tiff("spprichm1_space.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
plotter(posteriorsamp,'OverallNiDist','SppRich','overallnichedist','spprich',lineeq = lineequation,drawcor = T,roundtoy=0,xlab="Niche distance",ylab = "Species richness",linetype = "sdd")
dev.off()

# details for a table should they be required 
summary(spprichSP,pars="beta_p")$summary


# variance partition plots 

## write a function to do it fast across functions 
labelnames = c("Niche distance", "Site","Residuals")
varnames = c("OverallNiDist")
targetname = c("SppRich")
findplots = FindVarPart(varnames=varnames,targetname = targetname,labelnames = labelnames,postsamp = posteriorsamp,ymax=2)


setwd("")
tiff("SpprichVarpar.tif", res=600, compression = "lzw", height=5, width=15, units="in")
plot_grid(findplots[[1]],findplots[[2]] + theme(axis.title.y = element_blank()),findplots[[3]] + theme(axis.title.y = element_blank()) ,ncol=3,labels=c("a","b","c"))
dev.off()

findplots[[4]] #sp
findplots[[5]] #fi
findplots[[6]] #uk


