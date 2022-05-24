library(tidyverse)
library(geodist)

# load data
folderpath = "" # fill in location of folder, the 'here' package apparently does this more easily but I couldn't get it to work
filepath = "MechanismsUnderpinningCommunityStability/Data/CommunityButterflyData.csv"

CommunityData = read.csv(paste0(folderpath,filepath))

# split by country
CommunityDataSplit = split(CommunityData,CommunityData$site)

SpainMat = geodist(CommunityDataSplit$CAT[,c(3,4)],measure = "haversine") /1000
FinMat = geodist(CommunityDataSplit$FIN[,c(3,4)],measure = "haversine") /1000
UKMat = geodist(CommunityDataSplit$UK[,c(3,4)],measure = "haversine") /1000

maxval =  max(UKMat) # scaled by the largest distance
SpainMat = SpainMat / maxval 
FinMat =  FinMat / maxval 
UKMat= UKMat / maxval 


makelist <- function(dt){
  
  dt$CountryIndex <- as.integer(as.factor(dt$site)) 
  
  dt.list1<- list(
    
    # integers
    nData=nrow(dt),
    nCountry=length(unique(dt$CountryIndex)),
    
    Intercept = rep(1,nrow(dt)),
    
    # variables
    SppRich = scale(dt$spprich)[,1],
    SynIndex = scale(dt$synindex)[,1],
    SpeICOV = scale(dt$speciesICOV)[,1],
    LgMeanAb =scale(log(dt$meanab))[,1],
    LgMeanCom =scale(log(dt$meancommsize))[,1],
    ICOV = scale(dt$ICOV)[,1],
    
    #niche variables 
    Mismatch =scale(dt$mismatch)[,1],
    MeanVol =scale(dt$meanvol)[,1],
    OverallNiDist =scale(dt$overallnichedist)[,1],
    Jaccard= scale(dt$jaccard)[,1],
    
    # random
    CountryID = dt$CountryIndex, # 1=spain, 2=fin,3=uk
    SiteID = dt$nsite, # site given country
    
    UKMat = UKMat,
    FinMat =FinMat,
    SpainMat=SpainMat,
    
    nUK = dim(UKMat)[1],
    nFin = dim(FinMat)[1],
    nSpan = dim(SpainMat)[1]
    
  )
  
}



# make a list for stan
CommunityList=makelist(CommunityData)

# save for quick access 
#saveRDS(CommunityList,paste0(folderpath,"/MechanismsUnderpinningCommunityStability/Data/CommunityList.rds"))
