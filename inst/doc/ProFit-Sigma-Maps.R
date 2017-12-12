## ------------------------------------------------------------------------
library(ProFit)
library(ProFound)

## ------------------------------------------------------------------------
modellist=list(
  sersic=list(
    xcen=c(200, 200),
    ycen=c(200, 200),
    mag=c(17, 13),
    re=c(5, 10),
    nser=c(4, 1),
    ang=c(0, 135),
    axrat=c(1, 0.3),
    box=c(0,0)
  )
)

psf=profitMakeGaussianPSF(fwhm=3, dim=c(25,25))

image_elec=profitMakeModel(modellist=modellist, psf=psf, dim=c(400, 400), magzero=30)$z

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magimage(image_elec)
maghist(image_elec, log='x', majorn=c(10,5))

## ------------------------------------------------------------------------
RanPois=rpois(n=1e4, lambda=400)
RanPois_QE=rbinom(n=1e4, size=RanPois, prob=0.5)
RanPois_PureFilter=RanPois/2

sd(RanPois_QE)
sd(RanPois_PureFilter)

## ------------------------------------------------------------------------
sigma_elec=sqrt(image_elec)

## ------------------------------------------------------------------------
sky_elec=400
sky_elec_noise=sqrt(sky_elec)

image_elec=image_elec+sky_elec
sigma_elec=sqrt(sigma_elec^2+sky_elec_noise^2)

## ------------------------------------------------------------------------
dark_elec=9
dark_elec_noise=sqrt(dark_elec)

bias_elec=900
read_elec_noise=10
  
image_elec=image_elec+dark_elec+bias_elec
sigma_elec=sqrt(sigma_elec^2+dark_elec_noise^2+read_elec_noise^2)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
maghist(sqrt(image_elec)/sigma_elec)

## ------------------------------------------------------------------------
image_elec_sub=image_elec - sky_elec - dark_elec - bias_elec

## ------------------------------------------------------------------------
magzero=30
gain_elec_ADU=1/10^(-0.4*magzero)
image_ADU_sub=image_elec_sub/gain_elec_ADU

## ------------------------------------------------------------------------
image_ADU_sub_noise=image_ADU_sub+rnorm(length(image_ADU_sub), mean=0, sd=sigma_elec)/gain_elec_ADU

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magimage(image_ADU_sub_noise)

## ------------------------------------------------------------------------
sky_ADU=sky_elec/gain_elec_ADU
sigma_est_1=sqrt((image_ADU_sub_noise+sky_ADU)*gain_elec_ADU + read_elec_noise^2 + dark_elec_noise^2)/gain_elec_ADU

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
maghist((sigma_est_1/sigma_elec)*gain_elec_ADU)

## ------------------------------------------------------------------------
segim=image_elec>=quantile(image_elec,0.5)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magimage(segim)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
  sky_pix_ADU=image_ADU_sub_noise[segim==0]
maghist(sky_pix_ADU)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magplot(density(log10(sky_pix_ADU[sky_pix_ADU>0]), bw=0.01), col='red', type='l')
lines(density(log10(abs(sky_pix_ADU[sky_pix_ADU<0])), bw=0.01), col='blue')

## ------------------------------------------------------------------------
sky_pix_ADU_noise=sd(sky_pix_ADU)
sigma_est_2=sqrt(image_ADU_sub_noise/gain_elec_ADU+sky_pix_ADU_noise^2)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
maghist(sigma_est_2/sigma_elec*gain_elec_ADU)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
sigma_est_3=profoundMakeSigma(image=image_ADU_sub_noise, sky=0, skyRMS=sky_elec_noise, gain=gain_elec_ADU, readRMS=read_elec_noise, darkRMS=dark_elec_noise, sky_units='elec', read_units='elec', dark_units='elec')
magimage(sigma_est_3)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
sigma_est_4=profoundMakeSigma(image=image_ADU_sub_noise, sky=0, skyRMS=sky_pix_ADU_noise, gain=gain_elec_ADU)
magimage(sigma_est_4)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magplot(sigma_est_2, sigma_est_3, pch='.', log='xy', grid=TRUE)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magplot(sigma_est_3, sigma_est_4, pch='.', log='xy', grid=TRUE)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
sigma_est_5=profoundMakeSigma(image=image_ADU_sub_noise, objects=segim>0, sky=0, skyRMS=sky_pix_ADU_noise, gain=gain_elec_ADU)
magimage(sigma_est_5)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magplot(sigma_est_2, sigma_est_5,pch='.',log='xy', grid=TRUE)
abline(h=sky_pix_ADU_noise,col='red')

## ------------------------------------------------------------------------
sd((image_ADU_sub_noise-image_ADU_sub)/sigma_elec)*gain_elec_ADU

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magimage((image_ADU_sub_noise-image_ADU_sub)/sigma_est_1)
magimage((image_ADU_sub_noise-image_ADU_sub)/sigma_est_2)

magplot(density((image_ADU_sub_noise-image_ADU_sub)/sigma_est_1))
curve(dnorm, col='red', add=T)
magplot(density((image_ADU_sub_noise-image_ADU_sub)/sigma_est_2))
curve(dnorm, col='red', add=T)

sd((image_ADU_sub_noise-image_ADU_sub)/sigma_est_1)
sd((image_ADU_sub_noise-image_ADU_sub)/sigma_est_2)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
magimage((image_ADU_sub_noise*1.01-image_ADU_sub)/sigma_est_1)

## ---- fig.width=5, fig.height=5, dpi=40----------------------------------
maghist(sqrt(rpois(1e3,1e4)), breaks=20)
maghist(sqrt(rpois(1e3,1e2)), breaks=20)
maghist(sqrt(rpois(1e3,1e1)), breaks=20)

