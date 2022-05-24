
folderpath = "" 
functionpath = "MechanismsUnderpinningCommunityStability/Rcode Structural Equation Model/PlottingFunctions.R"
datapath = "MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"
stanpath = "MechanismsUnderpinningCommunityStability/Stan code DSeparation/"

source(paste0(folderpath,functionpath)) # will load required packages
CommunityList = read_rds(paste0(folderpath,datapath))

d16=stan(
 file = paste0(folderpath,stanpath,"m16dsep.stan"),
  data=CommunityList,
  iter=2000,
  chains=4,
  cores=4,
  control = list(adapt_delta = 0.999)
)
bci1 = extract.samples(d16)

# bayesian p 
#https://www.tandfonline.com/doi/full/10.1080/00031305.2020.1717621
bp1 = 2*(1- max(c(sum(bci1$ADSEP[,1] >0) / 4000, sum(bci1$ADSEP[,1] <0) /4000)))

#to save - commented out here
#dataoutpath = "MechanismsUnderpinningCommunityStability/Data"
#FischersC = readRDS(paste0(folderpath,dataoutpath, "FischersC.RData" ))
#nextindex = which.min(!sapply(FischersC,is.null))
#FischersC[[nextindex]] =c(bp1)
#saveRDS(FischersC,paste0(folderpath,dataoutpath, "FischersC.RData" ))