folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code Correlations/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,datapath))

# small list for the correlation
smalllist = list(
  N =  CommunityList$nData,
  x =  cbind(CommunityLis$Mismatch,CommunityLis$OverallNiDist)
)

set.seed(1234)
d1=stan(
    file = paste0(folderpath,stanpath,"correlationtests.stan"),
  data=smalllist,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)

stan_trace(d1, pars=c("rho", "mu", "sigma")) 

setwd("")
tiff("cor2MiOn.tif", res=600, compression = "lzw", height=6, width=6.5, units="in")
stan_dens(d1,"rho")+ xlim(-1,1) +xlab("Niche mismatch vs Niche distance")
dev.off()