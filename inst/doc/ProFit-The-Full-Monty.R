## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github('asgr/ProFound')
#  install_github('ICRAR/ProFit')

## ------------------------------------------------------------------------
evalglobal=FALSE

## ------------------------------------------------------------------------
library(ProFit)
library(ProFound)
library(RColorBrewer)

## ------------------------------------------------------------------------
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magimage(image, hicut=1)

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  segim=profoundMakeSegim(image, magzero=30, pixscale=0.34, plot=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  head(segim$segstats)

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  segim_expand=profoundMakeSegimExpand(image, segim$segim, expand=6, expandsigma=3, skycut=-1, magzero=30, pixscale=0.34, rotstats=TRUE, plot=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  SBlim=profitFlux2SB(segim_expand$skyRMS, magzero=30, pixscale=0.34)
#  print(SBlim)

## ---- eval=evalglobal----------------------------------------------------
#  object_frac=length(which(segim_expand$objects==1))/length(image)
#  qnorm(1-object_frac)

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  temp=qqnorm(sort((image[segim_expand$objects==0]-segim_expand$sky)/segim_expand$skyRMS), plot.it = FALSE)
#  magplot(0,0,xlim=c(-5,5),ylim=c(-5,5), grid=TRUE, type='n', xlab='Normal Quantile', ylab='Sky Pixel Quantile')
#  lines(temp)
#  abline(0,1,col='red')

## ---- eval=evalglobal----------------------------------------------------
#  gain=profoundGainEst(image,objects=segim_expand$objects, sky=segim_expand$sky, skyRMS=segim_expand$skyRMS)
#  gain

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  sigma=profoundMakeSigma(image, objects=segim_expand$objects, gain=gain, sky=segim_expand$sky, skyRMS =segim_expand$skyRMS, plot=TRUE)

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  magplot(segim_expand$segstats$R50, segim_expand$segstats$con,log='x', ylim=c(0,1), col=hsv(magmap(segim_expand$segstats$axrat, flip=TRUE)$map), xlab='Major Axis / pix', ylab='Concentration')
#  magbar('topleft', title='Axrat', titleshift=1)
#  abline(h=c(0.4,0.65), lty=2)
#  abline(v=c(0.7,1.2), lty=2)
#  
#  magplot(segim_expand$segstats$SB_N50, segim_expand$segstats$con, ylim=c(0,1), col=hsv(magmap(segim_expand$segstats$axrat, flip=TRUE)$map), xlab='SB[50] / mag/asec^2', ylab='Concentration', grid=TRUE)
#  magbar('topleft', title='Axrat', titleshift=1)
#  abline(v=c(23.5,2), lty=2)
#  
#  magplot( segim_expand$segstats$asymm, segim_expand$segstats$con, ylim=c(0,1), col=hsv(magmap(segim_expand$segstats$axrat, flip=TRUE)$map), xlab='Asymmetry', ylab='Concentration', grid=TRUE)
#  magbar('topleft', title='Axrat', titleshift=1)
#  abline(v=c(0.4), lty=2)

## ---- eval=evalglobal----------------------------------------------------
#  starlist=segim_expand$segstats[segim_expand$segstats$axrat>0.9 & segim_expand$segstats$con>0.4 & segim_expand$segstats$con<0.65 & segim_expand$segstats$R50>0.7 & segim_expand$segstats$R50<1.2 & segim_expand$segstats$SB_N50<23.5 & segim_expand$segstats$asymm<0.2,]
#  print(starlist)

## ---- eval=evalglobal----------------------------------------------------
#  psf_image=magcutout(image, loc=starlist[1,c('xcen','ycen')], box=c(31,31))
#  psf_sigma=magcutout(sigma, loc=starlist[1,c('xcen','ycen')], box=c(31,31))
#  psf_segim=magcutout(segim_expand$segim, loc=starlist[1,c('xcen','ycen')], box=c(31,31))

## ---- fig.width=5, fig.height=5, dpi=40, eval=evalglobal-----------------
#  magimage(psf_image$image)
#  magimage(psf_sigma$image)
#  magimage(psf_segim$image)

## ---- eval=evalglobal----------------------------------------------------
#  psf_x=psf_image$loc[1]
#  psf_y=psf_image$loc[2]
#  psf_mag=starlist[1,'mag']
#  psf_fwhm=starlist[1,'R50']*2/0.339
#  psf_con=1/starlist[1,'con']
#  
#  psf_modellist=list(
#    moffat=list(
#      xcen=psf_x,
#      ycen=psf_y,
#      mag=psf_mag,
#      fwhm=psf_fwhm,
#      con=psf_con,
#      axrat=1,
#      box=0
#    )
#  )

## ---- eval=evalglobal----------------------------------------------------
#  psf_tofit=list(
#    moffat=list(
#      xcen=TRUE,
#      ycen=TRUE,
#      mag=TRUE,
#      fwhm=TRUE,
#      con=TRUE,
#      axrat=FALSE,
#      box=FALSE
#    )
#  )
#  
#  psf_tolog=list(
#    moffat=list(
#      xcen=FALSE,
#      ycen=FALSE,
#      mag=FALSE,
#      fwhm=TRUE,
#      con=TRUE,
#      axrat=FALSE,
#      box=FALSE
#    )
#  )
#  
#  psf_intervals=list(
#    moffat=list(
#      xcen=list(psf_x+c(-5,5)),
#      ycen=list(psf_y+c(-5,5)),
#      mag=list(psf_mag+c(-2,2)),
#      fwhm=list(c(1,10)),
#      con=list(c(1,10)),
#      axrat=list(c(0.1,1)),
#      box=list(c(-1,1))
#    )
#  )

## ---- eval=evalglobal----------------------------------------------------
#  psf_Data=profitSetupData(psf_image$image, sigma=psf_sigma$image, modellist=psf_modellist, tofit=psf_tofit, tolog=psf_tolog, magzero=30, algo.func='optim', intervals=psf_intervals)

## ---- eval=evalglobal, fig.width=8, fig.height=5, dpi=40-----------------
#  profitLikeModel(parm=psf_Data$init, Data=psf_Data, makeplots=TRUE, plotchisq=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  psf_fit=optim(psf_Data$init, profitLikeModel, method='BFGS', Data=psf_Data, control=list(fnscale=-1))

## ---- eval=evalglobal, fig.width=8, fig.height=5, dpi=40-----------------
#  profitLikeModel(parm=psf_fit$par, Data=psf_Data, makeplots=TRUE, plotchisq=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  psf_modellist_fit=profitRemakeModellist(parm=psf_fit$par, Data=psf_Data)$modellist
#  psf_modellist_fit$moffat$xcen=25/2
#  psf_modellist_fit$moffat$ycen=25/2
#  psf_model=profitMakeModel(modellist=psf_modellist_fit, dim=c(25,25))$z

## ---- fig.width=5, fig.height=5, dpi=40, eval=evalglobal-----------------
#  magimage(psf_model)

## ---- eval=evalglobal----------------------------------------------------
#  psf_x=starlist$xcen
#  psf_y=starlist$ycen
#  psf_mag=starlist$mag
#  psf_fwhm=starlist[1,'R50']*2/0.339
#  psf_con=1/starlist[1,'con']
#  
#  psf_modellist2=list(
#    pointsource = list(
#      xcen = psf_x,
#      ycen = psf_y,
#      mag = psf_mag
#    ),
#    psf=list(
#      moffat=list(
#        mag=0,
#        fwhm=psf_fwhm,
#        con=psf_con,
#        axrat=1,
#        box=0
#      )
#    )
#  )

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  segim_expand_psf=profoundMakeSegimExpand(image, segim$segim, expand=starlist$segID, expandsigma=3, skycut=-1, magzero=30, pixscale=0.34, plot=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  psf_region=matrix(0, dim(image)[1], dim(image)[2])
#  psf_region[segim_expand$segim %in% starlist$segID]=1

## ---- eval=evalglobal----------------------------------------------------
#  psf_tofit2=list(
#    pointsource = list(
#      xcen = rep(TRUE,5),
#      ycen = rep(TRUE,5),
#      mag = rep(TRUE,5)
#    ),
#    psf=list(
#      moffat=list(
#        mag=FALSE,
#        fwhm=TRUE,
#        con=TRUE,
#        axrat=FALSE,
#        box=FALSE
#      )
#    )
#  )
#  
#  psf_tolog2=list(
#    pointsource = list(
#      xcen = rep(FALSE,5),
#      ycen = rep(FALSE,5),
#      mag = rep(FALSE,5)
#    ),
#    psf=list(
#      moffat=list(
#        mag=FALSE,
#        fwhm=TRUE,
#        con=TRUE,
#        axrat=TRUE,
#        box=TRUE
#      )
#    )
#  )
#  
#  psf_intervals2=list(
#    pointsource = list(
#      xcen=list(starlist$xcen[1]+c(-5,5), starlist$xcen[2]+c(-5,5), starlist$xcen[3]+c(-5,5), starlist$xcen[4]+c(-5,5), starlist$xcen[5]+c(-5,5)),
#      ycen=list(starlist$ycen[1]+c(-5,5), starlist$ycen[2]+c(-5,5), starlist$ycen[3]+c(-5,5), starlist$ycen[4]+c(-5,5), starlist$ycen[5]+c(-5,5)),
#      mag=list(starlist$mag[1]+c(-2,2), starlist$mag[2]+c(-2,2), starlist$mag[3]+c(-2,2), starlist$mag[4]+c(-2,2), starlist$mag[5]+c(-2,2))
#    ),
#    psf=list(
#      moffat=list(
#        mag=list(c(-1, 1)),
#        fwhm=list(c(1,10)),
#        con=list(c(1,10)),
#        axrat=list(c(0.1,1)),
#        box=list(c(-1,1))
#      )
#    )
#  )

## ---- eval=evalglobal----------------------------------------------------
#  psf_Data2=profitSetupData(image, sigma=sigma, modellist=psf_modellist2, tofit=psf_tofit2, tolog=psf_tolog2, magzero=30, algo.func='optim', intervals=psf_intervals2, region=psf_region)

## ---- eval=evalglobal, fig.width=8, fig.height=5, dpi=40-----------------
#  profitLikeModel(parm=psf_Data2$init, Data=psf_Data2, makeplots=TRUE, plotchisq=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  psf_fit2=optim(psf_Data2$init, profitLikeModel, method='BFGS', Data=psf_Data2, control=list(fnscale=-1))

## ---- eval=evalglobal, fig.width=8, fig.height=5, dpi=40-----------------
#  profitLikeModel(parm=psf_fit2$par, Data=psf_Data2, makeplots=TRUE, plotchisq=TRUE)

## ---- eval=evalglobal, fig.width=5, fig.height=5, dpi=40-----------------
#  psf_modellist_fit2=profitRemakeModellist(parm=psf_fit2$par, Data=psf_Data2)$modellist
#  psf_model_full2=profitMakeModel(psf_modellist_fit2,dim=dim(image),magzero=30)$z
#  for(i in 1:5){
#    magimage(magcutout(image, loc=starlist[i,c('xcen','ycen')], box=c(51,51))$image, col=rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100)), magmap=FALSE, zlim=c(-2e2,2e2))
#    magimage(magcutout(image-psf_model_full2, loc=starlist[i,c('xcen','ycen')], box=c(51,51))$image, col=rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100)), magmap=FALSE, zlim=c(-2e2,2e2))
#  }

## ---- eval=evalglobal----------------------------------------------------
#  psf_modellist_fit2$psf$moffat$xcen=25/2
#  psf_modellist_fit2$psf$moffat$ycen=25/2
#  psf_model2=profitMakeModel(modellist=psf_modellist_fit2$psf, magzero=30, dim=c(25,25))$z

## ---- fig.width=5, fig.height=5, dpi=40, eval=evalglobal-----------------
#  magimage(psf_model2)

## ---- eval=evalglobal----------------------------------------------------
#  gal_image=magcutout(image, loc=segim_expand$segstats[6,c('xcen','ycen')], box=c(101,101))
#  gal_sigma=magcutout(sigma, loc=segim_expand$segstats[6,c('xcen','ycen')], box=c(101,101))
#  gal_segim=magcutout(segim_expand$segim, loc=segim_expand$segstats[6,c('xcen','ycen')], box=c(101,101))

## ---- fig.width=5, fig.height=5, dpi=40, eval=evalglobal-----------------
#  magimage(gal_image$image)
#  magimage(gal_sigma$image)
#  magimage(gal_segim$image)

## ---- eval=evalglobal----------------------------------------------------
#  gal_x=gal_image$loc[1]
#  gal_y=gal_image$loc[2]
#  gal_mag=segim_expand$segstats[6,'mag']
#  gal_re=segim_expand$segstats[6,'R50']/0.339
#  gal_ang=segim_expand$segstats[6,'ang']
#  gal_axrat=segim_expand$segstats[6,'axrat']
#  
#  gal_modellist=list(
#    sersic=list(
#      xcen=rep(gal_x,2),
#      ycen=rep(gal_y,2),
#      mag=rep(gal_mag,2)+c(2,0.7),
#      re=c(gal_re/4,gal_re),
#      nser=c(4,1),
#      ang=c(0,gal_ang),
#      axrat=c(1,gal_axrat),
#      box=rep(0,2)
#    )
#  )

## ---- eval=evalglobal----------------------------------------------------
#  gal_tofit=list(
#    sersic=list(
#      xcen= c(TRUE,NA), #We fit for xcen and tie the two togther
#      ycen= c(TRUE,NA), #We fit for ycen and tie the two togther
#      mag= c(TRUE,TRUE), #Fit for both
#      re= c(TRUE,TRUE), #Fit for both
#      nser= c(FALSE,FALSE), #Fit for neither
#      ang= c(FALSE,TRUE), #Fit for disk
#      axrat= c(FALSE,TRUE), #Fit for disk
#      box= c(FALSE,FALSE) #Fit for neither
#    )
#  )
#  
#  gal_tolog=list(
#    sersic=list(
#      xcen= c(FALSE,FALSE),
#      ycen= c(FALSE,FALSE),
#      mag= c(FALSE,FALSE),
#      re= c(TRUE,TRUE), #re is best fit in log space
#      nser= c(TRUE,TRUE), #nser is best fit in log space
#      ang= c(FALSE,FALSE),
#      axrat= c(TRUE,TRUE), #axrat is best fit in log space
#      box= c(FALSE,FALSE)
#    )
#  )
#  
#  gal_intervals=list(
#    sersic=list(
#      xcen=list(lim=gal_x+c(-5,5),lim=gal_x+c(-5,5)),
#      ycen=list(lim=gal_y+c(-5,5),gal_y+c(-5,5)),
#      mag=list(lim=c(10,30),lim=c(10,30)),
#      re=list(lim=c(1,5),lim=c(1,100)),
#      nser=list(lim=c(0.5,20),lim=c(0.5,20)),
#      ang=list(lim=c(-180,360),lim=c(-180,360)),
#      axrat=list(lim=c(0.001,1),lim=c(0.001,1)),
#      box=list(lim=c(-1,1),lim=c(-1,1))
#    )
#  )

## ---- eval=evalglobal----------------------------------------------------
#  gal_Data=profitSetupData(gal_image$image, sigma=gal_sigma$image, modellist=gal_modellist, tofit=gal_tofit, tolog=gal_tolog, magzero=30, algo.func='optim', intervals=gal_intervals, psf=psf_model, segim=gal_segim$image)

## ---- eval=evalglobal, fig.width=8, fig.height=5, dpi=40-----------------
#  profitLikeModel(parm=gal_Data$init, Data=gal_Data, makeplots=TRUE, plotchisq=TRUE)

## ---- eval=evalglobal----------------------------------------------------
#  gal_fit=optim(gal_Data$init, profitLikeModel, method='BFGS', Data=gal_Data, control=list(fnscale=-1))

## ---- eval=evalglobal----------------------------------------------------
#  gal_fit_modellist=profitRemakeModellist(gal_fit$par, Data=gal_Data)$modellist
#  print(gal_fit_modellist)

## ---- eval=evalglobal, fig.width=8, fig.height=5, dpi=40-----------------
#  profitLikeModel(parm=gal_fit$par, Data=gal_Data, makeplots=TRUE, plotchisq=TRUE)
#  profitEllipsePlot(Data=gal_Data, modellist=gal_fit_modellist, pixscale=0.34, SBlim=SBlim)

