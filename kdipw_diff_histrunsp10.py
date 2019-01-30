# 2018 Jan 6:  select functions etc. from here for book code
# move phase generation functions to en_phase_screen.py
# begun 2015 Dec 17 based on kdipw2ds_large.py (c. Jan 2015) to study RISS

from preamble import *
from scipy.special  import fresnel

import scipy.ndimage as  ndimage
from scipy.stats import skew, kurtosis

# based on kdipw2ds.py from 2014 May 1

import datetime

import os.path

rcParams['font.size'] = 11

np.set_printoptions(threshold='nan')

#BELOW ARE THE THREE FUNCTIONS FROM JIM'S ORIGINAL CODE THAT RELATE TO THE HISTOGRAM

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
        si       = spectral index (e.g. 8/3 for a 1d Kolmogorov spectrum)
        phiF     = rms phase on the Fresnel scale
        rF       = Fresnel scale
        inner     = inner scale
        outer     = outer scale
        dx     = spatial sample interval
        xW     = width of screen
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
    #if debug: print 'targeted number of x samples = ', nx
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
    
    #print si, inner, outer, dx, npoints
    #print dq, qin, qout
    
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
        #print "index of fresnel scale = ", frindx
        #print var_fres_in, var_fres_out

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
    #print 'targeted number of x,y samples = ', nx,ny
    xvec = (arange(0.,nx)-nx/2+1)*dx
    yvec = (arange(0.,ny)-ny/2+1)*dy
    
    dqx = 2.*pi / xwidth
    dqy = 2.*pi / ywidth
    qmaxx = (2.*pi) / (2.*dx)
    qmaxy = (2.*pi) / (2.*dy)
    
    nqx = 2*int(qmaxx/dqx)
    nqy = 2*int(qmaxy/dqy)
    #print 'targeted number of q samples = ', nqx, nqy
    if nqx != nx:
        #print "Forcing nqx = nx = ", nx
        nqx = nx
    if nqy != ny:
        #print "Forcing nqy = ny = ", ny
        nqy = ny
    qxvec = (arange(0.,nqx)-nqx/2+1)*dqx
    qxvec = roll(qxvec,nqx/2+2)
    qyvec = (arange(0.,nqy)-nqy/2+1)*dqy
    qyvec = roll(qyvec,nqy/2+2)

    qin = 2.*pi / inner
    qout = 2.*pi / outer
    qshape = zeros((nqx, nqy))
    qshape_rolloff = zeros((nqx, nqy))      # for upper wavenumber rolloff
    #if apply_outer:
        #print "Applying outer-scale rolloff"

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
    #print si, inner, outer, dx, npoints
    #print dqx, dqy, qin, qout
    
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
        #print "index of fresnel scale = ", frindx
        #print var_fres_in, var_fres_out

# applying inner scale now an option with apply_inner = True (2015 Dec 28)
# needs to be applied *after* normalization!
# now need to recalculate the realization and apply norm_factor
    if apply_inner:
        #print "Applying inner-scale rolloff"
            
            # recalculate
        xform *= qshape_rolloff
        spectrum = abs(xform)**2
        xseries = real(ifft2(xform))
        xseries_norm = xseries * norm_factor

    return xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape


#------------------------------------------------------------------------------

#INSPIRED BY JIM'S CODE

# Main
#print "start: ", datetime.datetime.now()

do_gaussian = True
do_gaussian = False

do_diffraction = False
do_diffraction = True

# Gaussian parameters
phiblob = -1.        # < 0 ==> diverging
phiblob = -3.        # < 0 ==> diverging
phiblob =  0.3        # > 0 ==> converging
phiblob = -20.        # < 0 ==> diverging
phiblob = -30.        # < 0 ==> diverging
phiblob = -10.        # < 0 ==> diverging
phiblob = -0.3        # < 0 ==> diverging
phiblob = 0.

wex = 1.
wey = 2.

wex = 0.3
wey = 0.15

wex = 2
wey = 1

#DIFFRACTION PARAMETERS

si_kol_d = 11./3.
#phiF_d = 0.01
#phiF_d = 0.025
#phiF_d = 1.
#phiF_d = 2.
#phiF_d = 5.
phiF_d = 10.


rF = 1.             # Fresnel scale by definition
#rF = np.sqrt((37.5*1)/(2*np.pi))


# Note that array sizes scale as xwidth x ywidth
#xwidth = ywidth = 5.
xwidth = ywidth = 10.
#xwidth = ywidth = 20.
#xwidth = ywidth = 30.
#xwidth = ywidth = 200.
#xwidth = ywidth = 100.
#xwidth = ywidth = 75.
#xwidth = ywidth = 50.
#xwidth = ywidth = 64.

xmax = xwidth/2
ymax = ywidth/2

# get diffraction phase screen:
#define diffractive and refractive scales:
lscale_index = 2. / (si_kol_d-2.)
u = max(1., phiF_d**lscale_index)
ld = rF / u
lr = rF * u
inners = 2.*ld
inners = ld/2.
# outer scale = refraction scale for diffracting screen
# this isolates diffraction from refraction:
outers = lr

# Nyquist sample the Fresnel function and resolve the diffraction scale:
dx = min(pi * rF**2 / xmax, ld/2.)
dy = min(pi * rF**2 / ymax, ld/2.)
dx /= 2
dy /= 2

#print "diffractive phase: ", datetime.datetime.now()

# Diffractive phase:
xvec, yvec, xseries, phid, qxvec, qyvec, qshape = \
    gen_pl_phase_screen_2d(si_kol_d,phiF_d,rF,inners,outers,\
                                 dx,dy,xwidth,ywidth,normfres=True)

# subtract mean phase from screen (it doesn't really matter)
phid -= phid.mean()

nqx, nqy = shape(qshape)

gausstype = 'diverge'
if phiblob > 0.:  gausstype='converge'
plotfilepre = 'kdipw_riss_2016_2018' + \
    '%4.2f'%(si_kol_d) + \
    '_'+str(nqx)+'x'+str(nqy) + \
    '_rmsd_' + str(phiF_d*do_diffraction) + \
    '_g_' + gausstype + '_' + str(abs(phiblob)*do_gaussian)

# screen x,y array:
xy = meshgrid(xvec, yvec)
rsqvec = xy[0]**2 + xy[1]**2

# Gaussian blob at center:
gb = gen_gaussian_lens(phiblob, 0., 0., wex, wey, dx, dy, xwidth, ywidth)

xvecobs = copy(xvec)
xvecobshalf = xvec[where(abs(xvec) < xmax/2.)]

#print "propagation: ", datetime.datetime.now()

# propagator to observer:
# note definition of Fresnel scale implied is r_F^2 = \lambda D / 2\pi
#kernel0 = exp(1j*pi*rsqvec/rF**2)
kernel0 = exp(1j*rsqvec/(2.*rF**2))

# screen:
phi = zeros(shape(phid))
if do_gaussian:
    phi += gb

if do_diffraction:
    phi += phid

screen0 = exp(1j*phi)

kernel0fft = fft2(kernel0)
screen0fft= fft2(screen0)
field0fft = kernel0fft * screen0fft
field0 = ((dx*dy)/(2.*pi*rF**2)) * fftshift(ifft2(field0fft))
intensity0 = abs(field0)**2

#print "intensity0: ", intensity0

#intensity0 = abs(field0**2)
# nb: get same intensity0 using above two methods; kind of surprising
obsphase = arctan2(field0.imag, field0.real)
image = fftshift(abs(fft2(field0))) / shape(field0)[0] / pi  # norm'n = guess!
image /= image.max()
intensity0fft = fft2(intensity0)
acfI = fftshift(ifft2(abs(intensity0fft)**2)) / size(intensity0)
acfI = acfI.real
acfI -= 1.                # removes acf asymptotic value


#HISTOGRAM OF INTENSITIES

fig=figure()
xxx = copy(intensity0)
xxx.shape = (size(intensity0))
hist(xxx, bins=100, log=True, normed=True)
axis(ymin=1.e-10)
skewxxx = skew(xxx)

kurtxxx = kurtosis(xxx)
meanxxx = mean(xxx)
medianxxx = median(xxx)
stdxxx = std(xxx)
mIxxx = stdxxx / meanxxx
#print "Intensity stats:         mean        median        rms             skew           kurt             mI"
#print "                 ", meanxxx, medianxxx, stdxxx, skewxxx, kurtxxx, mIxxx
aa,bb,cc,dd = axis()
xpdf = arange(aa,bb, 0.01)

# Theoretical log-normal PDF
#lognormtheory_lab = r'$\rm Log-normal\ (theory)$'
#lnptheory = ((phiF_r, -phiF_r**2/2., 1.))
#lnpshape_theory = phiF_r
lnpshape_theory = mIxxx        # use measured mod index instead of theoretical
lnploc_theory = -lnpshape_theory**2/2.
lognormpdf_theory = (2.*pi*lnpshape_theory**2)**(-0.5) * (xpdf+1.e-10)**(-1.) * exp(-(log(xpdf+1.e-10)-lnploc_theory)**2/(2.*lnpshape_theory**2))
lntheory_line, = plot(xpdf, lognormpdf_theory, 'g-', lw=3)

#print "lognormpdf_theory: ", lognormpdf_theory

# Exponential PDF (Chi^2 with df = 2):
explab = r'$\rm Exponential \ (\chi_2^2)$'
expline, = plot(xpdf, exp(-xpdf), c = 'm', lw=2)

xlabel(r'$\rm Intensity$')
ylabel(r'$\rm PDF \ f_I(I)$')
#title(r'$\rm Histogram \ of \ Intensities \ from \  Screen \ \ \ DISS, Lens = (%d, %d)$'%(int(do_diffraction), int(do_gaussian)), fontsize=13)
title(r'$\rm Phi_d = 1.0 \ \ \ 1 \ Realization$')
annotate(r'$\rm {\bf Diffraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(phiF_d*float(do_diffraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)

#legend(plotlines, plotlabs, loc=4, prop={'size':12}, handlelength=3, labelspacing=0.01)
#axis(ymax=8)

#TURN ON OR OFF IF YOU WANT TO SEE THE GRAPHIC FOR A SINGLE REALIZATION

#show()

#SMOOTHING OF INTENSITY

lsmooth = lr / (2.*dx)        # smoothing length = one sigma
intensity0_smoothed = ndimage.gaussian_filter(intensity0, sigma=(lsmooth, lsmooth), order=0, mode='wrap')

#POTENTIAL PRINT STATEMENTS FOR A SINGLE REALIZATION OF INTENSITY STATISTICS

"""print " "
print "Raw intensity:"
print "      Max Intensity = ", intensity0.max()
print "     Mean intensity = ", intensity0.mean()
print "      RMS intensity = ", intensity0.std()
print "     Modulation index = ", intensity0.std() / intensity0.mean()

print " "
print "Smoothing scale = ", lsmooth, " samples"
print "Smoothed intensity"
print "      Max Intensity = ", intensity0_smoothed.max()
print "     Mean intensity = ", intensity0_smoothed.mean()
print "      RMS intensity = ", intensity0_smoothed.std()
print "     Modulation index = ", intensity0_smoothed.std() / intensity0_smoothed.mean()
print " "
"""
#print "Intensity stats:     max             max_smooth             mean               mean_smooth         rms             rms_smooth"
#print "                  ", intensity0.max(), intensity0_smoothed.max(), intensity0.mean(), intensity0_smoothed.mean(), intensity0.std(), intensity0_smoothed.std()

Intensity_stats = np.array([intensity0.max(), intensity0_smoothed.max(), intensity0.mean(), intensity0_smoothed.mean(), intensity0.std(), intensity0_smoothed.std()])

#SAVE DATA TO RESPECTIVE FILES

########################################################  INTENSITY STATS

if (os.path.isfile("Intensity_stats_25p10.txt")):
    with open("Intensity_stats_25p10.txt", "a") as stream:
        for i in range(1):
            out_string = "\n"
            out_string += str(Intensity_stats)
            stream.write(out_string)
else:
    with open("Intensity_stats_25p10.txt", "w") as stream:
        for i in range(1):
            out_string = ""
            out_string += str(Intensity_stats)
            stream.write(out_string)

###############################################################  XPDF

if (os.path.isfile("xpdf_25p10.txt")):
    with open("xpdf_25p10.txt", "a") as stream:
        for i in range(1):
            out_string = "\n"
            out_string += str(xpdf)
            stream.write(out_string)
else:
    with open("xpdf_25p10.txt", "w") as stream:
        for i in range(1):
            out_string = ""
            out_string += str(xpdf)
            stream.write(out_string)

################################################################  LOGNORMPDF_THEORY

if (os.path.isfile("lognormpdf_theory_25p10.txt")):
    with open("lognormpdf_theory_25p10.txt", "a") as stream:
        for i in range(1):
            out_string = "\n"
            out_string += str(lognormpdf_theory)
            stream.write(out_string)
else:
    with open("lognormpdf_theory_25p10.txt", "w") as stream:
        for i in range(1):
            out_string = ""
            out_string += str(lognormpdf_theory)
            stream.write(out_string)

#################################################################  XXX HISTOGRAM


if (os.path.isfile("xxx_25p10.txt")):
    with open("xxx_25p10.txt", "a") as stream:
        for i in range(1):
            out_string = "\n"
            out_string += str(xxx)
            stream.write(out_string)
else:
    with open("xxx_25p10.txt", "w") as stream:
        for i in range(1):
            out_string = "\n"
            out_string += str(xxx)
            stream.write(out_string)



#close('all')
