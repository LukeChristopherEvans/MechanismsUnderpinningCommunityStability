#plotting functions 
library(tidyverse)
library(ggsci)
library(cowplot)
library(RColorBrewer)
library(scales)
library(rethinking)
library(mgcv)

# color scheme
#cols<-pal_lancet("lanonc")(3)
#display.brewer.all(colorblindFriendly = T)
cols<-brewer.pal(3,"Set2")

# used in plotter
drawxy <- function(f,x,i) {
  mu <-sapply(x,f,i)
  mu_u <- apply(mu,2,mean)
  mupi<- apply(mu,2,PI)
  rbind(mu_u,mupi)
}

givemarg <- function(f,x,y,i) {
  marginals <-f(x,y,i)
  marginals
}


# plot distance function
distanceplot <-  function(posteriorsamp, PriToPost=F,VarName ="x") {
  
  plot(NULL,type="l",lwd=2,cex.lab=1.3,xlab=("Distance (km)"),ylab=paste(VarName,"Standarised covariance"),ylim=c(0,1),xlim=c(0,1000),
       bty = "n",bty = "n", xaxt = "n", yaxt = "n")
  box("plot", bty = "l", lwd = 2)
  axis(side = 1, lwd = 0, lwd.ticks = 2,cex.axis=1.3)
  axis(side = 2, lwd = 0, lwd.ticks = 2, las = 2,cex.axis=1.3)
  
  
  for(i in 1:sitenichel$nCountry) {
    
    distrange <- seq(from=0,to=1,length.out=100)
    pmcov<- sapply(distrange, function(x) posteriorsamp$etasq[,i]*exp(-posteriorsamp$rhosq[,i]*x^2) )
    pmcov_mu<-apply(pmcov,2,mean)
    
    # Destandardise
    tdistrange <- distrange * maxval
    
    lines(tdistrange,pmcov_mu,lwd=3,col=cols[i])
    
  }
  
  if (PriToPost==T) {
    plot(NULL,type="l",lwd=2,cex.lab=1.3,xlab=("Distance (km)"),ylab="Standarised covariance",ylim=c(0,3),xlim=c(0,1000),
         bty = "n",bty = "n", xaxt = "n", yaxt = "n")
    box("plot", bty = "l", lwd = 2)
    axis(side = 1, lwd = 0, lwd.ticks = 2,cex.axis=1.3)
    axis(side = 2, lwd = 0, lwd.ticks = 2, las = 2,cex.axis=1.3)  
    
    samp = rexp(50,1)
    distrange <- seq(from=0,to=1,length.out=100)
    pmcov<- sapply(distrange, function(x) samp*exp(-samp*x^2) )
    # Destandardise
    tdistrange <- distrange * maxval
    
    for(j in 1:nrow(pmcov)) {
      lines(tdistrange,pmcov[j,],lwd=0.4,col="grey")
    }
    
    distrange <- seq(from=0,to=1,length.out=100)
    for(i in 1:sitenichel$nCountry) {
      for (j in 1:50) {
        lines(tdistrange, posteriorsamp$etasq[j,i]*exp(-posteriorsamp$rhosq[j,i]*distrange^2),lwd=0.4,col=cols[i])
      }
    }
    
  }
}


# Spatial variance mapped
spatialvarplot <-  function(posteriorsamp=posteriorsamp) {
  
  library(maps)
  library(mapdata)
  par(mar=rep(0,4))
  Spm<-map('worldHires','Spain',xlim=c(-10,5),ylim=c(35,44))
  par(mar=rep(0,4))
  Ukm<-map('worldHires','UK',	xlim=c(-11,3.2), ylim=c(49,60.9))
  par(mar=rep(0,4))
  Fm<-map('worldHires','Finland')
  maplist<-list(Spm,Fm,Ukm)
  
  par(mfrow=c(1,3))
  for(i in 1:sitenichel$nCountry) {
    #i <- 3
    dimi<-ifelse(i==3,sitenichel$nUK,ifelse(i==2,sitenichel$nFin,sitenichel$nSpan))
    mat<-if(i==3) {
      sitenichel$UKMat
    } else if(i==2) {
      sitenichel$FinMat
    } else {
      sitenichel$SpainMat}
    
    K <-matrix(0,nrow = dimi,ncol=dimi)
    
    for (j in 1:dimi) { 
      for (k in 1:dimi){
        K[j,k] <- median(posteriorsamp$etasq[,i]) * exp(-median(posteriorsamp$rhosq[,i])*mat[j,k]^2)
      }
    } 
    
    diag(K) <-median(posteriorsamp$etasq[,i]+0.1)
    
    # correlation matrix 
    rhos <-round(cov2cor(K),2)
    
    # long lat
    lonlat<-if(i==3) {
      dplyr::filter(siteniches2,site=="UK")[,c(4,5)]
    } else if(i==2) {
      dplyr::filter(siteniches2,site=="FIN")[,c(4,5)]
    } else {
      dplyr::filter(siteniches2,site=="CAT")[,c(4,5)]}
    
    
    par(mar=rep(0,4))
    map(maplist[[i]])
    points(lonlat$transect_lon,lonlat$transect_lat,pch=19,cex=1.1,col=cols[i])
    
    for(j in 1:dimi){
      for(k in 1:dimi){
       lines(c(lonlat$transect_lon[j],lonlat$transect_lon[k]),c(lonlat$transect_lat[j],lonlat$transect_lat[k]),lwd=0.3,col=col.alpha(cols[i], rhos[j,k]^2 *0.5) )
      }
    }
  }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1)) # put margins back
  par(mfrow=c(1,1))
}

# posterior plot function
plotter <- function(posteriorsamp,stanx,stany,dtx,dty,lineeq=T,margeq=F, drawcor=F ,xlab="x",ylab="y",m1=NA,m2=NA,m3=NA,m4=NA,roundtox=1,roundtoy=1,linetype=c("sss"),yl=0,yu=1) {
  
  stanxvals<-sitenichel[[stanx]]
  stanyvals<-sitenichel[[stany]]
  dtxvals <- siteniches2[[dtx]]
  dtyvals <- siteniches2[[dty]]
  
  #marginal values 
  m1val<-sitenichel[[m1]]
  m2val<-sitenichel[[m2]]
  m3val<-sitenichel[[m3]]
  m4val<-sitenichel[[m4]]
  
  # full x
  mxf<-max(stanxvals)
  mif<-min(stanxvals)
  xfrange <- seq(from=mif,to=mxf,length.out=10)
  
  # full y 
  mxy<-max(stanyvals)
  miy<-min(stanyvals)
  yrange <- seq(from=miy,to=mxy,length.out=10)# replaced with zero here
  
  
  # axis conversion
  truerange <- (xfrange ) *sd(dtxvals) + mean(dtxvals)
  truey <- yrange *sd(dtyvals) + mean(dtyvals)
  
  ylab<-if (typeof(margeq) != "logical" ){paste("Residual",ylab)
  } else {ylab}
  if (typeof(margeq) != "logical" ){
    plot(NULL,xlab=xlab,ylab=ylab, axes=F,xlim=c(mif,mxf),ylim=c(yl,yu),
         bty = "n",bty = "n", xaxt = "n", yaxt = "n",cex.lab=1.3) 
  } else {
    plot(NULL,xlab=xlab,ylab=ylab, axes=F,xlim=c(mif,mxf),ylim=c(miy,mxy),
         bty = "n",bty = "n", xaxt = "n", yaxt = "n",cex.lab=1.3) 
  }
  
  box("plot", bty = "l", lwd = 2)
  axis(1, at=xfrange,labels=round(truerange,roundtox),lwd = 0, lwd.ticks = 2,cex.axis=1.3)
  
  
  for(i in 1:sitenichel$nCountry) {
    if (typeof(lineeq) != "logical" ) {
    xhere<- stanxvals[sitenichel$CountryID ==i]
    yhere<- stanyvals[sitenichel$CountryID ==i]
    # segment x
    mx<-max(xhere)
    mi<-min(xhere)
    xrange <- seq(from=mi,to=mx,length.out=50)
    
    linedata<-drawxy(lineeq,xrange,i)
    if(substr(linetype,i,i) == "s") {
    lines(xrange,linedata[1,],type="l",col=col.alpha(cols[i], 0.7),lwd=5)
    } else {
      lines(xrange,linedata[1,],type="l",lty="dashed",col=col.alpha(cols[i], 0.7),lwd=5)
    }
    lines(xrange,linedata[2,],lty="dotted",col=col.alpha(cols[i], 0.7),lwd=2)
    lines(xrange,linedata[3,],lty="dotted",col=col.alpha(cols[i], 0.7),lwd=2)
    
    }
    # spatial correlation
    if (drawcor ==T) {
      dimi<-ifelse(i==3,sitenichel$nUK,ifelse(i==2,sitenichel$nFin,sitenichel$nSpan))
      mat<-if(i==3) {
        sitenichel$UKMat
      } else if(i==2) {
        sitenichel$FinMat
      } else {
        sitenichel$SpainMat}
      
      K <-matrix(0,nrow = dimi,ncol=dimi)
      
      for (j in 1:dimi) { 
        for (k in 1:dimi){
          K[j,k] <- median(posteriorsamp$etasq[,i]) * exp(-median(posteriorsamp$rhosq[,i]*mat[j,k]))
        }
      } 
      
      diag(K) <-median(posteriorsamp$etasq[,i]+0.1)
      
      # correlation matrix 
      rhos <-round(cov2cor(K),2)
      
      for(j in 1:dimi){
        for(k in 1:dimi){
          lines(c(xhere[j],xhere[k]),c(yhere[j],yhere[k]),lwd=0.2,col=col.alpha(cols[i], rhos[j,k]^2 *0.5 ) )
        }
      }
    }
    
    
    yrangemarg<-c()
    if(typeof(margeq) != "logical" ) {
      
      x1here<- m1val[sitenichel$CountryID ==i]
      x2here<- m2val[sitenichel$CountryID ==i]
      x3here<- m3val[sitenichel$CountryID ==i]
      x4here<- m4val[sitenichel$CountryID ==i]
      xshere<-cbind(x1here,x2here,x3here,x4here)
      yhere<- stanyvals[sitenichel$CountryID ==i]
      ymarg<-givemarg(margeq,xshere,yhere,i)
      points(xhere,ymarg,col=col.alpha(cols[i], 0.7),pch=19,cex=1.5)
      yrangemarg<-c(yrangemarg,ymarg)
      
      
    } else {
      
      points(stanxvals,stanyvals,col=alpha(cols[sitenichel$CountryID],0.8/3),pch=19,cex=1.5) #they will go on top of one another - but doesn't matter
      
    }
    
  }
  
  # add marginal axis
  if(typeof(margeq) == "logical" ) {
    axis(2, at=yrange,labels=round(truey,roundtoy),lwd = 0, lwd.ticks = 2, las = 2,cex.axis=1.3)
   }  else {
    mxym<-max(yrangemarg)
    miym<-min(yrangemarg)
    join = c(mxym,miym)
    sv = join[which.max(abs(join))]
    
    yrangemargin <- c(seq(from=sv,to=0,length.out=3), 0 + diff(seq(from=sv,0,length.out=3))[1], 0 + 2 * diff(seq(from=sv,0,length.out=3))[1])  # replaced with zero here
    truey2 <- yrangemargin *sd(dtyvals)
    axis(2, at=c(yrangemargin,0),labels=round(c(truey2,0.0),roundtoy),lwd = 0, lwd.ticks = 2, las = 2,cex.axis=1.3)
    
  }
}

# coefficient plotter
coefplotter<- function(model=m, coefnames, type="random") {
  
  countries <- c("Spain","Finland","UK")
  ncoefs<-length(coefnames)
  newnames<-c()
  
  if (type=="random") {
    for(i in 1:ncoefs) {
      for(j in countries) {
        newnames <-c(newnames, paste(coefnames[i],j))
      }
    }
    
    for(i in 1:ncoefs) {
      for(j in 1:3){
        namehere <- paste0("beta_p","[",i,",",j,"]")  
        names(model)[names(model) ==namehere] = newnames[1]
        newnames<-newnames[-1]
      }
    }
    
    plot(model, pars = c('beta_p')) +
      theme(legend.position = "none",axis.text=element_text(size=20),
            axis.title=element_text(size=20)) +
      xlab("Effect size (scaled)") +
      ylab("Coefficient ") +
      geom_vline(xintercept = 0,lty="dashed")
    
  } else {
    
    for(i in 1:ncoefs) {
      newnames <-c(newnames, coefnames[i])
    }
    
    for(i in 1:ncoefs) {
      namehere <- paste0("hyperbeta","[",i,"]")  
      names(model)[names(model) ==namehere] = newnames[1]
      newnames<-newnames[-1]
    }
    
    plot(model, pars = c('hyperbeta')) +
      theme(legend.position = "none",axis.text=element_text(size=20),
            axis.title=element_text(size=20)) +
      xlab("Effect size (scaled)") +
      ylab("Coefficient ") +
      geom_vline(xintercept = 0,lty="dashed")
    
  }
}

# function used for dsep tests
bcis = function(obj){
  e1=extract(obj,pars="ADSEP")
  res=c()
  for(i in 1:4) {
    s1=e1$ADSEP[,i]
    med1=median(s1)
    dist=abs(med1-0)
    bp=length(s1[s1 < med1 - dist | s1 > med1+dist]) / length(s1)
    res =c(res,bp)
  } 
  return(res)
}


# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
# Split violin plots
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



# Variance paritioning 
# code adapted from 
# https://www.biorxiv.org/content/10.1101/2021.10.17.464682v2.full.pdf

library(fastDummies)
library(Matrix)
library(tidyverse)
library(cowplot)

construcLinMat = function(params,Xcov,y,sitedummymat,siteoff){
  ACov <- as.matrix(Xcov %*% Diagonal(x = params[colnames(Xcov)])) 
  colnames(ACov) <- colnames(Xcov)
  AIid <- as.matrix(sitedummymat) %*% siteoff[params[length(params)],] # row id
  colnames(AIid) <- c('site')
  A <- cbind(ACov, AIid)
  ## Linear term for residuals (omiting the intercept above does not affect variance)
  A <- cbind(A, residuals = y - rowSums(A))
  return(as.matrix(A)) # the product to return 
  ####
}

# this par
lin_term_cov <- function(params, Xcov, y,sitedummymat,siteoff) {
  A <- construcLinMat(params, Xcov, y,sitedummymat,siteoff)
  K <- cov(A)
  dimnames(K) <- list(colnames(A), colnames(A))
  return(K)
}



colsextended<-brewer.pal(3,"Dark2")



FindVarPart = function(varnames,targetname,labelnames,postsamp,ymax=1.1,box=T){
  varnamesfull = c(varnames ,"CountryID","SiteID")
  dg =  cbind.data.frame(sitenichel[varnamesfull],sitenichel[targetname])
  # make country selection
  plotlist=vector(mode="list",length=6) # first 3 plots second three tables
  
  for(i in 1:3){
    
    tdg = filter(dg, CountryID==unique(dg$CountryID)[i])
    xvars = dplyr::select(tdg,all_of(varnames))
    
    # target
    y = dplyr::select(tdg,all_of(targetname))
    names(y) ="y"
    # params 
    params = postsamp$beta_p[,,i]
    colnames(params) = c("Intercept",varnames)
    
    # random effects for spain
    sitelabs = c("spK","finK","ukK")[i]
    siteoff=postsamp[[sitelabs]]
    
    # dummy matrix for site 
    sitedummymat =dummy_cols(data.frame(site= as.factor(tdg$SiteID) ))[,-1]
    
    # drop intercept 
    XCov <- as.matrix(xvars)
    params = cbind(params, ID=1:4000)
    
    # run the big function
    bigK = apply(params,1,lin_term_cov,Xcov=XCov,y=y,sitedummymat=sitedummymat,siteoff=siteoff)
    # sort
    ncols = length(varnames) + 2 # site and resids
    nvarnames = c(varnames,"site","residuals")
    dim(bigK) <- c(ncols, ncols, 4000)
    dimnames(bigK) <- list(rows = nvarnames, cols = nvarnames, iterations = NULL)
    
    # variance partition
    V_x <- apply(bigK,  'iterations', diag) / var(y$y)
    # V_x has samples as columns and rows as variables 
    dimnames(V_x) <- list(variance = nvarnames, iterations = NULL)
    
    # Diagonal variance partition
    P_x <- apply(bigK, 'iterations', function(x) {diag(x) / sum(diag(x))})
    # P_x has samples as columns and rows as variables 
    
    # correlation between var parition 
    cor_V_x <- cor(aperm(V_x, c('iterations', 'variance')))
    
    # plotting with gg
    dVx= pivot_longer(data.frame(t(V_x)),everything(), names_to="variable",values_to = "variance")
    dVx$type = "Vx"
    dPx= pivot_longer(data.frame(t(P_x)),everything(), names_to="variable",values_to = "variance")
    dPx$type = "Px"
    
    dVx =rbind.data.frame(dVx,dPx)
    
    dVx$variable=factor(dVx$variable,levels =c(nvarnames))
    dVx$type=factor(dVx$type,levels =c("Vx","Px"))
    
    if(box ==F){
      plotlist[[i]] = dVx %>%
        ggplot(aes(x=variable,y=variance,fill=type)) +geom_split_violin() +
        theme_cowplot() +
        xlab("")+
        ylab("Density")+
        ylim(0,ymax)+
        scale_x_discrete(labels = labelnames)+
        scale_fill_manual(values = c(cols[i],colsextended[i]) )+
      theme(legend.position = "none")
    } else {
      plotlist[[i]] = dVx %>%
        ggplot(aes(x=variable,y=variance,fill=type)) +geom_boxplot(outlier.alpha = 0.1) +
        theme_cowplot() +
        xlab("")+
        ylab("Density")+
        ylim(0,ymax)+
        scale_x_discrete(labels = labelnames)+
        scale_fill_manual(values = c(cols[i],colsextended[i]) )+
        theme(legend.position = "none")
    }
    
    # summary table
    t1=dVx %>% group_by(type,variable)%>%
      summarise( mean = mean(variance),
                 var = var(variance), 
                 min = min(variance),
                 q25 = quantile(variance, 0.25), 
                 median = median(variance), 
                 q75 = quantile(variance, 0.75), 
                 max = max(variance))
    plotlist[[i+3]] = t1
    
    
  }
  
  return(plotlist) 
}
