# based on kdipw2ds.py from 2014 May 1
from preamble import *
from scipy.special  import fresnel

rcParams['font.size'] = 11

#------------------------------------------------------------------------------

def gen_gaussian_lens(philens, xcen, ycen, wex, wey , dx, dy, xW, yW):
    """
    Generates phase screen for a Gaussian lens centered.

    philens         = maximum phase 
    xcen, ycen      = center position of lens 
    wex, wey        = 1/e widths of lens
    dx, dy          = x and y intervals
    xW, yW          = width of screen in x and y 
    """ 

    # x,y grid
    nx = int(xW/dx)
    ny = int(yW/dy)
    xvec = (arange(0.,nx)-nx/2+1)*dx
    yvec = (arange(0.,ny)-ny/2+1)*dy

    glens = zeros((xvec.size, yvec.size))
    
    for i, x in enumerate(xvec):
        for j, y in enumerate(yvec):
             arg = ((x-xcen)/wex)**2 +  ((y-ycen)/wey)**2
             glens[i,j] = exp(-arg)
    glens *= philens
    return glens

#------------------------------------------------------------------------------

def gen_pl_phase_screen_1d(si, phiF, rF, inner, outer, dx, xW, normfres=True, debug = False):
    """
    Generates a phase screen comprising a realization of power-law noise 
    with spectral index si,  inner and outer scales as specified.

    Length scales can be in any units so long as they are all the same units
    (rF, inner, outer, dx, xW)

    Input:
       si   	= spectral index (e.g. 8/3 for a 1d Kolmogorov spectrum)
       phiF 	= rms phase on the Fresnel scale
       rF   	= Fresnel scale  
       inner 	= inner scale
       outer 	= outer scale
       dx 	= spatial sample interval
       xW 	= width of screen
       normfres = True   to set the rms phase equal to phiF
                  False  for the screen to have unit variance 

    returns:
       xvec = x axis along screen (perpendicular to line of sight) 
       xseries = screen phase 
       xeries_norm = screen phase scaled to input rms phase on the Fresnel scale
       qvec = wavenumber axis
       shape = sqrt(shape of wavenumber spectrum)

    """
    # specified diffraction and refraction scales
    ld = rF / phiF 
    lr = rF * phiF 

    nx = int(xW/dx)
    if debug: print 'targeted number of x samples = ', nx
    xvec = (arange(0.,nx)-nx/2+1)*dx

    dq = 2.*pi /xW 
    #qmax = (2.*pi) / (2.*dx)
    #nq = 2*int(qmax/dq)

    #if debug:print 'targeted number of q samples = ', nq

    #if nq != nx: 
        #print "Forcing nq = nx = ", nx
        #nq = nx
    nq = nx
    qvec = (arange(0.,nq)-nq/2+1)*dq
    qvec = roll(qvec,nq/2+1)
    qin = 2.*pi / inner
    qout = 2.*pi / outer
    shape = (qout**2 + qvec**2)**(-si/4.) * exp(-abs(qvec)/(2.*qin)) 
    npoints = size(shape)

    print si, inner, outer, dx, npoints
    print dq, qin, qout

    xformr=randn(npoints)*shape
    xformi=randn(npoints)*shape
    xform = xformr + 1j*xformi
    spectrum=real(xform*conj(xform))
    xseries = real(ifft(xform))

    if normfres:
       frindx = int(rF/dx)
       var_fres_in = var(xseries[0:size(xseries)-frindx]-xseries[frindx:])
       xseries_norm = xseries * phiF/ sqrt(var_fres_in) 
       var_fres_out = var(xseries_norm[0:size(xseries_norm)-frindx]-xseries_norm[frindx:])
       print "index of fresnel scale = ", frindx
       print var_fres_in, var_fres_out

    return xvec, xseries, xseries_norm, qvec, shape

#------------------------------------------------------------------------------

def gen_pl_phase_screen_2d(si, phiF, rF, inner, outer, dx, dy, xwidth, ywidth, apply_inner=False,  apply_outer=False, normfres=True, debug=False):
    """
    Generates npoints of a realization of power-law noise with unit 
    variance with spectral index si and inner and outer scales as 
    specified for a sample interval dx. 


    input:
    si = spectral index of power-law wavenumber spectrum
    phiF = rms phase at Fresnel scale (rad)
        length scales: all dimensionless:
        rF      = Fresnel scale
        inner, outer    = inner and outer scales
        dx, dy      = sample intervals
        xwidth, ywidth  = screen extent
    logical:
            normfres    = True implies normalization to phiF

        Definition of Fresnel scale: r_F^2 = \lambda D / 2\pi

    returns:
       xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape
       (xvec, yvec) coordinates perpendicular to line of sight 
       xseries = screen phase 
       xseries_norm = screen phase scaled to input rms phase on Fresnel scale
       qxvec, qyvec  = wavenumber axes
       qshape = sqrt(shape of wavenumber spectrum)

    NB. previous name: def gen_powphase2d( ... )

    """

    global spectrum, qshape_rolloff

    nx = int(xwidth/dx)
    #ny = nx
    ny = int(ywidth/dy)
    print 'targeted number of x,y samples = ', nx,ny
    xvec = (arange(0.,nx)-nx/2+1)*dx
    yvec = (arange(0.,ny)-ny/2+1)*dy

    dqx = 2.*pi / xwidth
    dqy = 2.*pi / ywidth
    qmaxx = (2.*pi) / (2.*dx)
    qmaxy = (2.*pi) / (2.*dy)

    nqx = 2*int(qmaxx/dqx)
    nqy = 2*int(qmaxy/dqy)
    print 'targeted number of q samples = ', nqx, nqy 
    if nqx != nx: 
        print "Forcing nqx = nx = ", nx
        nqx = nx
    if nqy != ny: 
        print "Forcing nqy = ny = ", ny
        nqy = ny
    qxvec = (arange(0.,nqx)-nqx/2+1)*dqx
    qxvec = roll(qxvec,nqx/2+2)
    qyvec = (arange(0.,nqy)-nqy/2+1)*dqy
    qyvec = roll(qyvec,nqy/2+2)

    qin = 2.*pi / inner
    qout = 2.*pi / outer
    qshape = zeros((nqx, nqy))
    qshape_rolloff = zeros((nqx, nqy))      # for upper wavenumber rolloff
    if apply_outer: 
        print "Applying outer-scale rolloff"
    
    # 2016 Jan 1: put in upper wavenumber cutoff at array size
    # to avoid aliasing
    qmax = qxvec.max()/2.
    for i, qxi in enumerate(qxvec):
       for j, qyj in enumerate(qyvec):
          qsq = qxi**2 + qyj**2
          qshape[i,j] = (qout**2 + qsq)**(-si/4.) * exp(-qsq/(2.*qmax**2)) 
          qshape_rolloff[i,j] = exp(-qsq / (2.*qin**2))
          if apply_outer: 
             qshape[i,j] *= exp(-qout**2 / (2.*qsq))

    npoints = size(qshape)
    print si, inner, outer, dx, npoints
    print dqx, dqy, qin, qout

    ## new 2016 Jan 1:  create real white noise in x domain and FFT
    ## to get Hermitian noise
    ##xx = randn(nqx, nqy)      # real
    ##xq = fft(xx)          # Hermitian 
    xformr=randn(nqx, nqy)*qshape
    xformi=randn(nqx, nqy)*qshape
    xform = xformr + 1j*xformi
    ##xform = xq * qshape           # Hermitian
    spectrum = abs(xform)**2
    ##xseries = ifft2(xform).real
    xseries = real(ifft2(xform))

    # Normalization factor needs to be calculated on pure power-law spectrum
    # before any rolloff at the refraction scale
    if normfres:
       frindx = int(rF/dx)
       x1dcut = xseries[0,:]
       var_fres_in = var(x1dcut[0:size(x1dcut)-frindx]-x1dcut[frindx:])
       norm_factor = phiF / sqrt(var_fres_in)
       xseries_norm = xseries * norm_factor
       xn1dcut = xseries_norm[0,:]
       var_fres_out = var(xn1dcut[0:size(xn1dcut)-frindx]-xn1dcut[frindx:])
       print "index of fresnel scale = ", frindx
       print var_fres_in, var_fres_out

    # applying inner scale now an option with apply_inner = True (2015 Dec 28)
    # needs to be applied *after* normalization!
    # now need to recalculate the realization and apply norm_factor 
    if apply_inner: 
       print "Applying inner-scale rolloff"

       # recalculate 
       xform *= qshape_rolloff
       spectrum = abs(xform)**2
       xseries = real(ifft2(xform))
       xseries_norm = xseries * norm_factor

    return xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape

#------------------------------------------------------------------------------

def gen_powphase2d_orig(si, phiF, rF, inner, outer, dx, dy, xwidth, ywidth, apply_inner=False,  apply_outer=False, normfres=True, debug=False):
    """
    Generates npoints of a realization of power-law noise with unit 
    variance with spectral index si and inner and outer scales as 
    specified for a sample interval dx. 


    input:
    si = spectral index of power-law wavenumber spectrum
    phiF = rms phase at Fresnel scale (rad)
        length scales: all dimensionless:
        rF      = Fresnel scale
        inner, outer    = inner and outer scales
        dx, dy      = sample intervals
        xwidth, ywidth  = screen extent
    logical:
            normfres    = True implies normalization to phiF

        Definition of Fresnel scale: r_F^2 = \lambda D / 2\pi

    returns:
       xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape
       (xvec, yvec) coordinates perpendicular to line of sight 
       xseries = screen phase 
       xseries_norm = screen phase scaled to input rms phase on Fresnel scale
       qxvec, qyvec  = wavenumber axes
       qshape = sqrt(shape of wavenumber spectrum)

    """

    global spectrum, qshape_rolloff

    nx = int(xwidth/dx)
    #ny = nx
    ny = int(ywidth/dy)
    print 'targeted number of x,y samples = ', nx,ny
    xvec = (arange(0.,nx)-nx/2+1)*dx
    yvec = (arange(0.,ny)-ny/2+1)*dy

    dqx = 2.*pi / xwidth
    dqy = 2.*pi / ywidth
    qmaxx = (2.*pi) / (2.*dx)
    qmaxy = (2.*pi) / (2.*dy)

    nqx = 2*int(qmaxx/dqx)
    nqy = 2*int(qmaxy/dqy)
    print 'targeted number of q samples = ', nqx, nqy 
    if nqx != nx: 
        print "Forcing nqx = nx = ", nx
        nqx = nx
    if nqy != ny: 
        print "Forcing nqy = ny = ", ny
        nqy = ny
    qxvec = (arange(0.,nqx)-nqx/2+1)*dqx
    qxvec = roll(qxvec,nqx/2+2)
    qyvec = (arange(0.,nqy)-nqy/2+1)*dqy
    qyvec = roll(qyvec,nqy/2+2)

    qin = 2.*pi / inner
    qout = 2.*pi / outer
    qshape = zeros((nqx, nqy))
    qshape_rolloff = zeros((nqx, nqy))      # for upper wavenumber rolloff
    if apply_outer: 
        print "Applying outer-scale rolloff"
    
    for i, qxi in enumerate(qxvec):
       for j, qyj in enumerate(qyvec):
          qsq = qxi**2 + qyj**2
          qshape[i,j] = (qout**2 + qsq)**(-si/4.) 
          qshape_rolloff[i,j] = exp(-qsq / (2.*qin**2))
          if apply_outer: 
             qshape[i,j] *= exp(-qout**2 / (2.*qsq))

    npoints = size(qshape)
    print si, inner, outer, dx, npoints
    print dqx, dqy, qin, qout

    xformr=randn(nqx, nqy)*qshape
    xformi=randn(nqx, nqy)*qshape
    xform = xformr + 1j*xformi
    spectrum = abs(xform)**2
    xseries = real(ifft2(xform))

    # Normalization factor needs to be calculated on pure power-law spectrum
    # before any rolloff at the refraction scale
    if normfres:
       frindx = int(rF/dx)
       x1dcut = xseries[0,:]
       var_fres_in = var(x1dcut[0:size(x1dcut)-frindx]-x1dcut[frindx:])
       norm_factor = phiF / sqrt(var_fres_in)
       xseries_norm = xseries * norm_factor
       xn1dcut = xseries_norm[0,:]
       var_fres_out = var(xn1dcut[0:size(xn1dcut)-frindx]-xn1dcut[frindx:])
       print "index of fresnel scale = ", frindx
       print var_fres_in, var_fres_out

    # applying inner scale now an option with apply_inner = True (2015 Dec 28)
    # needs to be applied *after* normalization!
    # now need to recalculate the realization and apply norm_factor 
    if apply_inner: 
       print "Applying inner-scale rolloff"

       # recalculate 
       xform *= qshape_rolloff
       spectrum = abs(xform)**2
       xseries = real(ifft2(xform))
       xseries_norm = xseries * norm_factor

    return xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape

#------------------------------------------------------------------------------

def gen_powphase2d_old(si, phiF, rF, inner, outer, dx, dy, xW, yW,  normfres=True, debug=True):
    """
    Generates npoints of a realization of power-law noise with unit 
    variance with spectral index si and inner and outer scales as specified 
    for a sample interval dx. 

    normfres = True implies normalization so that the rms value on the
    Fresnel scale is set to rF

    returns:
       (xvec, yvec) = x,y axes in screen (perpendicular to line of sight) 
       xyseries = screen phase 
       xyeries_norm = screen phase scaled to input rms phase on Fresnel scale
       qxvec, qyvec  = wavenumber axes
       qshape = sqrt(shape of wavenumber spectrum)

    """
    # specified diffraction and refraction scales
    ld = rF / phiF  
    lr = rF * phiF 

    nx = int(xW/dx)
    ny = nx
    if debug: print 'targeted number of x,y samples = ', nx,ny
    xvec = (arange(0.,nx)-nx/2+1)*dx
    yvec = (arange(0.,ny)-ny/2+1)*dy

    dqx = 2.*pi / xW 
    dqy = 2.*pi / yW
    qmaxx = (2.*pi) / (2.*dx)
    qmaxy = (2.*pi) / (2.*dy)

    nqx = 2*int(qmaxx/dqx)
    nqy = 2*int(qmaxy/dqy)
    if debug: print 'targeted number of q samples = ', nqx, nqy 
    if nqx != nx: 
        print "Forcing nqx = nx = ", nx
        nqx = nx
    if nqy != ny: 
        print "Forcing nqy = ny = ", ny
        nqy = ny
    qxvec = (arange(0.,nqx)-nqx/2+1)*dqx
    qxvec = roll(qxvec,nqx/2+1)
    qyvec = (arange(0.,nqy)-nqy/2+1)*dqy
    qyvec = roll(qyvec,nqy/2+1)

    qin = 2.*pi / inner
    qout = 2.*pi / outer
    qshape = zeros((nqx, nqy))
    
    for i, qxi in enumerate(qxvec):
       for j, qyj in enumerate(qyvec):
          qsq = qxi**2 + qyj**2
          qshape[i,j] = (qout**2 + qsq)**(-si/4.) 
          #qshape[i,j] = (qout**2 + qsq)**(-si/4.) * exp(-(qsq/(2.*qin**2))) 
    npoints = size(qshape)

    if debug:
        print si, inner, outer, dx, npoints
        print dqx, dqy, qin, qout

    xformr=randn(nqx, nqy)*qshape
    xformi=randn(nqx, nqy)*qshape
    xform = xformr + 1j*xformi
    spectrum=real(xform*conj(xform))
    xseries = real(ifft2(xform))

    if normfres:
       frindx = int(rF/dx)
       x1dcut = xseries[0,:]
       var_fres_in = var(x1dcut[0:size(x1dcut)-frindx]-x1dcut[frindx:])
       xseries_norm = xseries * rF / sqrt(var_fres_in) 
       xn1dcut = xseries_norm[0,:]
       var_fres_out = var(xn1dcut[0:size(xn1dcut)-frindx]-xn1dcut[frindx:])
       #var_fres_out = var(xseries_norm[0:size(xseries_norm)-frindx]-xseries_norm[frindx:])
       print "index of fresnel scale = ", frindx
       print var_fres_in, var_fres_out

    return xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape

#------------------------------------------------------------------------------
