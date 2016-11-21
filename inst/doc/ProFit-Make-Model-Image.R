## ------------------------------------------------------------------------
library(ProFit)

## ------------------------------------------------------------------------
ExampleImage0=profitMakeModel(list(sersic=list(xcen=50, ycen=50, re=4, nser=4, ang=30, axrat=0.5)))
str(ExampleImage0)

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  image(ExampleImage0)

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  magimage(ExampleImage0)

## ------------------------------------------------------------------------
model1 = list(
  sersic = list(
		xcen   = c(180, 60),
		ycen   = c(90, 10),
		mag = c(15, 13),
		re  = c(14, 5),
		nser  = c(3, 10),
		ang  = c(46, 80),
		axrat  = c(0.4, 0.6),
		box = c(0.5,-0.5)
	),
	pointsource = list(
		xcen   = c(34.3,10.1,150),
		ycen   = c(74.2,120.5,130.4),
		mag = c(10,13,16)
	),
	sky = list(
		bg = 3e-12
	)
)

## ------------------------------------------------------------------------
ExampleImage1=profitMakeModel(model=model1, dim=c(200,200))
str(ExampleImage1)

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  magimage(ExampleImage1)

## ------------------------------------------------------------------------
ExamplePSF = profitMakePointSource()
ExamplePSF = ExamplePSF/sum(ExamplePSF)
ExampleImagePSF1=profitMakeModel(model=model1, psf=ExamplePSF, dim=c(200,200))

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  magimage(ExampleImagePSF1)

## ---- fig.width=5, fig.height=5------------------------------------------
model1a=model1
model1a$psf = list(sersic=list(nser=0.5,mag=0,re=1,axrat=1,ang=0))
ExampleImagePSF1A=profitMakeModel(model=model1a, psf=ExamplePSF, dim=c(200,200))

## ------------------------------------------------------------------------
modelBD = list(
	sersic = list(
		xcen   = c(100, 100),
		ycen   = c(100, 100),
		mag = c(14, 12),
		re  = c(2, 15),
		nser  = c(4, 1),
		ang  = c(0, 60),
		axrat  = c(1, 0.3),
		box = c(0,0)
	)
)
ExampleImage = profitMakeModel(model=modelBD, psf=ExamplePSF, dim=c(200,200))

## ------------------------------------------------------------------------
magimage(ExampleImage)

## ------------------------------------------------------------------------
Finesample=3L
FinePSFmodel = model1a$psf
FinePSFmodel$sersic$re = FinePSFmodel$sersic$re*Finesample
ExamplePSFFine = profitMakePointSource(model=FinePSFmodel,image = matrix(0,27,27))
ExampleImageFine = profitMakeModel(model=modelBD, dim=c(200,200),finesample = Finesample)

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  magimage(ExampleImageFine$z-ExampleImage$z)

## ------------------------------------------------------------------------
model2 = list(
	sersic = list(
		xcen   = runif(20,0,200),
		ycen   = runif(20,0,200),
		mag = runif(20,15,20),
		re  = runif(20,1,100),
		nser  = runif(20,0.5,8),
		ang  = runif(20,0,180),
		axrat  = runif(20,0.3,1),
		box = runif(20,-0.3,0.3)
	),
	pointsource = list(
		xcen   = runif(10,0,200),
		ycen   = runif(10,0,200),
		mag = runif(10,15,20)
	),
	sky = list(
		bg = 3e-12
	)
)

## ---- fig.width=5, fig.height=5------------------------------------------
system.time(ExampleImagePSF2<-profitMakeModel(model=model2, psf=ExamplePSF, dim=c(200,200)))
magimage(ExampleImagePSF2)

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  model3 = list(
#  	sersic = list(
#  		xcen   = runif(20,0,1000),
#  		ycen   = runif(20,0,1000),
#  		mag = runif(20,15,20),
#  		re  = runif(20,1,100),
#  		nser  = runif(20,0.5,8),
#  		ang  = runif(20,0,180),
#  		axrat  = runif(20,0.3,1),
#  		box = runif(20,-0.3,0.3)
#  	),
#  	pointsource = list(
#  		xcen   = runif(10,0,1000),
#  		ycen   = runif(10,0,1000),
#  		mag = runif(10,15,20)
#  	),
#  	sky = list(
#  		bg = 3e-12
#  	)
#  )
#  system.time(ExampleImagePSF3<-profitMakeModel(model=model3, psf=ExamplePSF, dim=c(1000,1000)))
#  magimage(ExampleImagePSF3)

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  system.time(profitMakeModel(model=model3, psf=ExamplePSF, dim=c(1000,1000), remax=10))

## ---- fig.width=5, fig.height=5, eval=F----------------------------------
#  ExampleImageAdd=profitAddMats(ExampleImagePSF3$z, ExampleImagePSF2$z, c(300,400))
#  magimage(ExampleImageAdd)

