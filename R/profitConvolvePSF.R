profitConvolvePSF=function(image, psf, calcregion, docalcregion=FALSE, 
  options=list(method="Bruteconv"), sky=0, plot=FALSE, ...){
  if(sky!=0){
    image=image-sky
  }
  if(missing(calcregion)){
    if(docalcregion){
      calcregion=matrix(1,dim(image)[1],dim(image)[2])
    }else{
      calcregion=matrix(1,1,1)
    }
  }
  
  if(all(dim(calcregion)==dim(image))==FALSE & docalcregion) {
    stop(paste("calcregion dimensions are ",dim(calcregion)[1],":",dim(calcregion)[2]," and they must be ",dim(image)[1],":",dim(image)[2],"!",sep=""))
  }
  
  if(dim(psf)[1]%%2==0 | dim(psf)[1]%%2==0){
    xrange=floor(-dim(psf)[1]/2):ceiling(dim(psf)[1]/2)
    yrange=floor(-dim(psf)[2]/2):ceiling(dim(psf)[2]/2)
    regrid=expand.grid(xrange,yrange)
    psf=matrix(profitInterp2d(regrid[,1],regrid[,2],psf)[,3],length(xrange),length(yrange))
  }
  psf=psf/sum(psf)
  isbc1 = options$method == "Bruteconv"
  isfftr = options$method == "FFTconv"
  isfftw = options$method == "FFTWconv"
  if(isfftw & !requireNamespace("fftw", quietly = TRUE)){
    stop('The fftw package is needed for the FFTWconv option to work. Please install it from CRAN.', call. = FALSE)
  }
  isfft = isfftr || isfftw
  if(isbc1)
  {
    output=profitBruteConv(image,psf,calcregion,docalcregion)
  } else if(isfft)
  {
    if(isfftr) {
      psffft = options$fft$psf$r
      if(is.null(psffft))
      {
        psfpad = matrix(0,options$fft$paddim[1],options$fft$paddim[2])
        psfpad[options$fft$psf$x,options$fft$psf$y] = psf
        psfpad
        psffft = fft(psfpad)
      }
    } else if(isfftw) {
      psffft = options$fft$psf$w
      if(is.null(psffft))
      {
        psfpad = matrix(0,options$fft$paddim[1],options$fft$paddim[2])
        psfpad[options$fft$psf$x,options$fft$psf$y] = psf
        if(is.null(options$fft$fftwplan)) psffft = fftw::FFT(psfpad)
        else psffft = fftw::FFT(psfpad,fftwplan=options$fft$fftwplan)
      }
    }
    imagepad = matrix(0,options$fft$paddim[1],options$fft$paddim[2])
    imagepad[options$fft$padimagex,options$fft$padimagey] = image
    if(isfftw) {
      if(is.null(options$fft$fftwplan)) imagepad = fftw::FFT(imagepad)
      else imagepad = fftw::FFT(imagepad, plan=options$fft$fftwplan)
    } else if(isfftr) {
      imagepad = fft(imagepad, inverse = FALSE)
    }
    imagepad = imagepad * psffft
    if(isfftw) {
      if(is.null(options$fft$fftwplan)) imagepad = fftw::IFFT(imagepad)
      else imagepad = fftw::IFFT(imagepad, plan=options$fft$fftwplan)
      dim(imagepad) = options$fft$paddim
    } else if(isfftr) {
      imagepad = fft(imagepad, inverse = TRUE)/options$fft$paddim[1]/options$fft$paddim[2]
    }
    output=Re(imagepad)[options$fft$cropx,options$fft$cropy]
  }
  if(sky!=0){
    output=output+sky
  }
  
  if(plot){
	  magimage(output, ...)
	}
	
  return=output
}