
folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code DSeparation/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,datapath))

# start a list 
setres = vector(mode='list',length=19)

set.seed(1234)
d1= stan(
  file = paste0(folderpath,stanpath,"m1dsep.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.99)
)

bci1 = extract.samples(d1)

# bayesian p 
#https://www.tandfonline.com/doi/full/10.1080/00031305.2020.1717621
bp1 = 2*(1- max(c(sum(bci1$ADSEP[,1] > 0) / 4000, sum(bci1$ADSEP[,1] < 0) /4000)) )

nextindex = which.min(!sapply(setres,is.null))

setres[[nextindex]] =c(bp1)


dataoutpath = "MechanismsUnderpinningCommunityStability/Data"
saveRDS(setres,paste0(folderpath,dataoutpath, "setres.RData" ))
