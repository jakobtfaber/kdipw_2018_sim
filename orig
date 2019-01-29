# 2018 Jan 6:  select functions etc. from here for book code
# move phase generation functions to en_phase_screen.py
# begun 2015 Dec 17 based on kdipw2ds_large.py (c. Jan 2015) to study RISS

from preamble import *
from scipy.special  import fresnel

import scipy.ndimage as  ndimage
from scipy.stats import skew, kurtosis

import gen_phase_screen as genps

#from make_2d import *

#from rebin import *



rcParams['font.size'] = 11
 
import datetime
"""
def make_2d(x,y):
"""
"""Change from 1-d indexing to 2-d indexing
(translated from IDL to Python).
Convert an N element X vector, and an M element Y vector, into
N x M arrays giving all possible combination of X and Y pairs.
Useful for obtaining the X and Y positions of each element of
a regular grid.
CALLING SEQUENCE:
xx,yy = make_2d.make_2d(x,y)
INPUTS:
x - N element vector of X positions
y - M element vector of Y positions
RETURNS:
xx - N x M element array giving the X position at each pixel
yy - N x M element array giving the Y position of each pixel
If only 2 parameters are supplied then X and Y will be
updated to contain the output arrays
EXAMPLE:
To obtain the X and Y position of each element of a 30 x 15 array
import make_2d
x = numpy.arange(30)  ;  y = numpy.arange(15)
xx,yy = make_2d.make_2d( x, y )
REVISION HISTORY:
Written                     Wayne Landsman,ST Systems Co.    May,            1988
Added /NOZERO keyword       W. Landsman                      March,          1991
Converted to IDL V5.0       W. Landsman                      September,      1997
Improved speed              P. Broos                         July,           2000
Converted to Python         D. Jones                         January,        2014
"""
"""


ny = len(y)
nx = len(x)
xx = x.reshape(1,nx)
yy = y.reshape(ny,1)

xx = rebin.rebin(xx, [ny, nx])
yy = rebin.rebin(yy, [ny, nx])

return(xx,yy)

def rebin(a, new_shape):
    
M, N = a.shape
    m, n = new_shape
        
        if m<M:
return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
else:
return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)
"""

def make_2d(x,y):
    xx = np.zeros((x.shape[0],x.shape[0])) #initializes xx to proper shape but full of 0s
    xx[:,:]=x #sets the columns of xx equal to x
    yy = np.zeros((y.shape[0],y.shape[0])) #initializes yy to proper shape but full of 0s
    y = y.tolist()
    y.reverse()
    y = np.asarray(y)
    y.shape = (y.shape[0],1) #turns y into a (2d) column
    yy[:,:]=y #sets the rows of yy equal to y
    return xx,yy

#------------------------------------------------------------------------------

def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return nr, radial_prof

#------------------------------------------------------------------------------

def plot_wavenumber_spectrum():

    """
    Plots the 2D wavenumber spectrum of the phase 
    """
    #global qxvec, qyvec, qshape

    # plot imposed power spectrum of phase (qshape**2)
    # roll spectral shape so origin is in middle:
    qqxvec = copy(qxvec)
    qqyvec = copy(qyvec)
    qqshape = copy(qshape)
    ncen = size(qqxvec)/2 -1
    qqxvec = roll(qqxvec, ncen)
    qqyvec = roll(qqyvec, ncen)
    qqshape = roll(roll(qqshape, ncen, axis=0), ncen, axis=1)

    fig=figure()
    imshow(2.*log10(qqshape), origin='lower', extent = (qqxvec.min(), qqxvec.max(), qqyvec.min(), qqyvec.max()), interpolation='nearest')
    xlabel(r'$ q_x$', fontsize=15)
    ylabel(r'$ q_y$', fontsize=15)
    colorbar()
    title(r'$\rm 2D \ Phase \ Power \ Spectrum \ (log10 \ color \ levels)\ \ \ \beta = %4.2f$'%(si_kol_d))

    #inset to show spectrum slice on log-log scale
    a = axes([0.22e0,0.650, 0.20, 0.225], frameon=True)
    indices = where(qqxvec > 0.)[0]
    plot(qqxvec[indices], qqshape[indices, ncen]**2, 'k.')
    xscale('log')
    yscale('log')
    xlabel(r'$q_x$', color='w', fontsize=14)
    ylabel(r'$P_{\delta n_e}$', color='w', fontsize=14)
    tick_params(axis='both', which='minor', bottom='off', left='off', top='off', right='off')
    tick_params(axis='both', which='major', labelcolor='white', labelsize=8)
    show()
    savefig('kdi_spectrum_'+str(si_kol_d)+'_'+str(nqx)+'x'+str(nqy)+'.pdf')
    return

#------------------------------------------------------------------------------

def plot_fourframeI():

    """
    Plotting: four frame I:   
        diffraction phase 
        refraction phase 
        real part of screen kernel
        intensity pattern 
    """

    # Plotting: four frame I
    fig=figure()
    subplots_adjust(hspace=0.2)

    subplot(221)
    imshow(phid, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    ylabel(r'$Y$')
    tick_params(axis='x', labelbottom='off')
    title(r'$\rm Diffraction \ Phase \ \  \phi_d$')
    tick_params(axis='both', labelsize='11')
    colorbar()

    subplot(222)
    imshow(phir, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    tick_params(axis='x', labelbottom='off')
    title(r'$\rm Refraction \  Phase \ \  \phi_r$')
    tick_params(axis='both', labelsize='11')
    colorbar()

    subplot(223)
    imshow(screen0.real,aspect='auto',origin='lower',extent=(-xmax,xmax,-ymax,ymax), interpolation='nearest')
    xlabel(r'$X$')
    ylabel(r'$Y$')

    phi_string = ''
    n_in_string = 0
    if do_diffraction: 
       n_in_string += 1
       phi_string += '\phi_d'
    if do_refraction: 
       n_in_string += 1
       if n_in_string > 1: phi_string += '+'
       phi_string += '\phi_r'
    if do_gaussian: 
       n_in_string += 1
       if n_in_string > 1: phi_string += '+'
       phi_string += '\phi_g'
    if n_in_string > 1: phi_string = '(' + phi_string + ')' 

    title(r'$\rm Screen \ \  {\cal R} \{ e^{i %s}\}$'%(phi_string))
    tick_params(axis='both', labelsize='11')
    colorbar()

    subplot(224)
    imshow(log10(intensity0), aspect='auto', origin='lower',extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    xlabel(r'$X$')
    title(r'$\rm log10 \ (Observed \ Intensity)$')
    tick_params(axis='both', labelsize='11')
    colorbar()

    annotate(r'$\rm {\bf Diffraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Refraction}: \ \phi_F = %5.2f \ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(phiF_d*float(do_diffraction), phiF_r*float(do_refraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)
    
    show()
    savefig(plotfilepre + '_fourframesI' + '.pdf')
    return

#------------------------------------------------------------------------------

def plot_fourframeII():

    """
    Plotting: four frame II:   
         screen phase 
         intensity pattern
         image 
         intensity ACF
    """

    fig=figure()
    subplots_adjust(hspace=0.35, wspace=0.28)

    subplot(221)
    imshow(phi, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    xlabel(r'$X$')
    ylabel(r'$Y$')
    #tick_params(axis='x', labelbottom='off')
    phi_string = ''
    n_in_string = 0
    if do_diffraction: 
       n_in_string += 1
       phi_string += '\phi_d'
    if do_refraction: 
       n_in_string += 1
       if n_in_string > 1: phi_string += '+'
       phi_string += '\phi_r'
    if do_gaussian: 
       n_in_string += 1
       if n_in_string > 1: phi_string += '+'
       phi_string += '\phi_g'
    if n_in_string > 1: phi_string = '(' + phi_string + ')' 

    #title(r'$\rm Screen \ \  {\cal R} \{ e^{i %s}\}$'%(phi_string))
    title(r'$\rm Screen\  Phase \ \ \phi$')
    tick_params(axis='both', labelsize='11')
    colorbar()

    subplot(222)
    imshow(log10(intensity0), aspect='auto', origin='lower',extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    xlabel(r'$X$')
    title(r'$\rm log10 \ (Observed \ Intensity)$')
    tick_params(axis='both', labelsize='11')
    colorbar()

    subplot(223)
    imshow(image, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    xlabel(r'$\theta_x$')
    ylabel(r'$\theta_y$')
    #tick_params(axis='x', labelbottom='off')
    title(r'$\rm Image \ \ I(\theta)$')
    tick_params(axis='both', labelsize='11')
    colorbar()

    subplot(224)
    # zoom in to 1/4 of frame
    zoom = 16 
    nxselect = where(abs(xvec)<=xmax/zoom)[0]	# where gives  tuple 
    nyselect = where(abs(yvec)<=ymax/zoom)[0] 
    nzoom = size(nxselect)
    acfIplot = zeros((nzoom, nzoom))
    """for nzx, nx in enumerate(nxselect):
       for nzy, ny in enumerate(nyselect):"""
    nx, ny = make_2d(nxselect, nyselect)
    #acfIplot[nzx,nzy] = acfI[nx,ny]
    acfIplot = acfI
    #acfIplot = acfI[where(abs(xvec)<=xmax/zoom), where(abs(xvec)<=xmax/zoom)]
    imshow(acfIplot,aspect='auto',origin='lower',extent=(-xmax/zoom,xmax/zoom,-ymax/zoom,ymax/zoom), interpolation='nearest')
    xlabel(r'$\rm Lag \ \ \delta X$')
    ylabel(r'$\rm Lag \ \ \delta Y$')

    title(r'$\rm \delta(Intensity) \ ACF $')
    tick_params(axis='both', labelsize='11')
    colorbar()

    annotate(r'$\rm {\bf Diffraction+Refraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(phiF_r*float(do_refraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)

    show()
    savefig(plotfilepre + '_fourframesII' + '.pdf')
    #raw_input('hit return')
    #close()
    return

#------------------------------------------------------------------------------

def plot_Ihistogram(dofits=False):

    # histogram of intensities
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
    print "Intensity stats:         mean        median        rms             skew           kurt             mI"
    print "                 ", meanxxx, medianxxx, stdxxx, skewxxx, kurtxxx, mIxxx
    aa,bb,cc,dd = axis()
    xpdf = arange(aa,bb, 0.01)

    # Theoretical log-normal PDF
    lognormtheory_lab = r'$\rm Log-normal\ (theory)$'
    lnptheory = ((phiF_r, -phiF_r**2/2., 1.))
    #lnpshape_theory = phiF_r
    lnpshape_theory = mIxxx		# use measured mod index instead of theoretical
    lnploc_theory = -lnpshape_theory**2/2.
    lognormpdf_theory = (2.*pi*lnpshape_theory**2)**(-0.5) * (xpdf+1.e-10)**(-1.) * exp(-(log(xpdf+1.e-10)-lnploc_theory)**2/(2.*lnpshape_theory**2))
    lntheory_line, = plot(xpdf, lognormpdf_theory, 'g-', lw=3)

    # Exponential PDF (Chi^2 with df = 2):
    explab = r'$\rm Exponential \ (\chi_2^2)$'
    expline, = plot(xpdf, exp(-xpdf), 'c-', lw=2)

    if dofits:
        # Best-fit Chi^2 PDF:
        chi2parms = chi2.fit(xxx)
        chi2xxx = chi2.pdf(xpdf, chi2parms[0], chi2parms[1], chi2parms[2])
        chi2lab = r'$\rm \chi^2 $'
        chi2line, = plot(xpdf, chi2xxx, 'r-', lw=2)
        print " "
        print "Chi2parms: ", chi2parms

        # Best-fit log-normal PDF:
        lognparms = lognorm.fit(xxx)
        lognormpdf = lognorm.pdf(xpdf, lognparms[0], lognparms[1], lognparms[2])
        lognormlab = r'$\rm Log-normal \ (fit)$'
        lognormline, = plot(xpdf, lognormpdf, 'g-',  lw=2)
        print " "
        print "Lognormparms: ", lognparms


        # Best-fit Normal PDF:
        normparms = norm.fit(xxx)
        normpdf = norm.pdf(xpdf, normparms[0], normparms[1])
        normlab = r'$\rm Normal $'
        normline, = plot(xpdf, normpdf, 'g--',  lw=2)
        print " "
        print "Normparms: ", normparms

    #
    # Best fit parameters for chosen pdf:
    #    norm.fit(xxx)
    #    lognorm.fit(xxx)
    #    chi2.fit(xxx)  returns shape, loc, scale
    #                   (4.4515167090247481, -0.0018131689491289039, 0.22495835196637748) 
    #    xvec = arange(0., 10., 0.01)
    #    chi2xxx = chi2.pdf(xvec, 4.4515167090247481, -0.0018131689491289039, 0.22495835196637748)
    #    plot(xvec, chi2xxx)
    #    seems to match up well with histogram
    # Note: plot(xvec, chi2.pdf(xvec, 2., 0., 0.5)) gives an exponential PDF
    #       not sure why the scale = last parameter = 0.5
    #       looks like the first parameter = degrees of freedom
    # Best fit log-normal:
    muhat = mean(log(xxx))
    varhat = mean((log(xxx)-muhat)**2)
    # or
    # see help(stats.lognorm):
    #     sigma = shape parameter
    #     exp(mu) = scale parameter
    # see help(lognorm.fit)
    # lognorm.fit(xxx) returns MLE for shape, location, and scale parameters
    #     shape = sigma
    #     location = mu
    #     scale parameter = ?
    # can calculate analytical PDF using:
    # logpdf=lognorm.logpdf(xvec, sqrt(varhat), muhat)   (not sure if sqrt is right)
    # or
    # pdf=lognorm.pdf(xvec, varhat, muhat)
    # to get moments of chi2 PDF:
    # chi2.stats(1,2,moments='mvsk')
    # returns (array(1.0), array(4.0), array(0.0), array(0.0))
    # see http://www.johndcook.com/blog/distributions_scipy/
    xlabel(r'$\rm Intensity$')
    ylabel(r'$\rm PDF \ f_I(I)$')
    title(r'$\rm Histogram \ of \ Intensities \ from \  Screen \ \ \ DISS, RISS, Lens = (%d, %d, %d)$'%(int(do_diffraction), int(do_refraction), int(do_gaussian)), fontsize=13)
    annotate(r'$\rm {\bf Diffraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Refraction}: \ \phi_F = %5.2f \ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(phiF_d*float(do_diffraction), phiF_r*float(do_refraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)

    if dofits:
       plotlines=((chi2line, expline, lognormline, lntheory_line, normline))
       plotlabs=((chi2lab, explab, lognormlab, lognormtheory_lab, normlab))
    else:
       plotlines=((expline, lntheory_line))
       plotlabs=(( explab, lognormtheory_lab))

    legend(plotlines, plotlabs, loc=4, prop={'size':12}, handlelength=3, labelspacing=0.01)
    show()
    savefig(plotfilepre + '_histogram' + '.pdf')
    return

#------------------------------------------------------------------------------

def plot_intensity_slice():

    # Intensity slice
    fig=figure()
    plot(xvec, intensity0[shape(intensity0)[0]/2])
    xlabel(r'$\rm X$')
    ylabel(r'$\rm Intensity$')
    title(r'$\rm Slice \ Through \ Spatial \ Intensity \ Pattern \  from \  Screen \ \ \ DISS, RISS, Lens = (%d, %d, %d)$'%(int(do_diffraction), int(do_refraction), int(do_gaussian)), fontsize=13)
    annotate(r'$\rm {\bf Diffraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Refraction}: \ \phi_F = %5.2f \ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(phiF_d*float(do_diffraction), phiF_r*float(do_refraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)
    show()
    savefig(plotfilepre + '_intensity_slice' + '.pdf')
    #raw_input('hit return')
    #close()
    return

#------------------------------------------------------------------------------

def plot_wavefield_phase():

    # Wavefield phase at screen
    fig=figure()
    plot(xvec, obsphase[shape(obsphase)[0]/2])
    xlabel(r'$\rm X$')
    ylabel(r'$\rm Wavefield \ Phase \ \ (rad)$')
    title(r'$\rm Slice \ Through \ Phase \ of\  Wavefield\ Pattern  \ \ \ DISS, RISS, Lens = (%d, %d, %d)$'%(int(do_diffraction), int(do_refraction), int(do_gaussian)), fontsize=13)
    annotate(r'$\rm {\bf Diffraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Refraction}: \ \phi_F = %5.2f \ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(phiF_d*float(do_diffraction), phiF_r*float(do_refraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)
    show()
    savefig(plotfilepre + '_wavefield_phase' + '.pdf')
    return

#------------------------------------------------------------------------------

def plot_1D_intensity_spectrum():

    fig=figure()             

    # this isn't quite right: summing over 1d spectra and multiplying
    # by 1d wavenumber; the low wavenumber portion of the spectrum is then
    # too steep (linear instead of q**(1/3.)) but I don't know why.
    # all told my worries about inconsistency with theoretical spectra
    # seem to be due to the azimuthal summing of the spectrum. 
    plot(qxvec[0:size(qxvec)/2], qxvec[0:size(qxvec)/2]*Ispectrum[0:size(qxvec)/2])
    #plot(qxvec[0:size(qxvec)/2], Ispectrum[0:size(qxvec)/2])
    xscale('log')
    yscale('log')
    axis(ymin=Ispectrum.max()/1.e10)
    xlabel(r'$\rm Wavenumber \ q \ \ (rad \ per \  Fresnel \ scale) $')
    ylabel(r'$\rm Spectrum $')
    show()
    savefig(plotfilepre + '_intensity_power_spectrum' + '.pdf')
    return

# -----------------------------------------------------------------------------

def plot_spectral_shapes():

    figure()
    plot(qxvec[0:size(qxvec)/2], qshape[0:size(qxvec)/2, 0]**2)
    plot(qxvec[0:size(qxvec)/2], genps.qshape_rolloff[0:size(qxvec)/2, 0]**2)
    plot(qxvec[0:size(qxvec)/2], (qshape[0:size(qxvec)/2, 0]*genps.qshape_rolloff[0:size(qxvec)/2, 0])**2, lw=2)
    plot(qxvec[0:size(qxvec)/2], genps.spectrum[0:size(qxvec)/2, 0])
    xscale('log')
    yscale('log')
    axis(ymin=qshape.max()/1.e10)
    show()
    return

# -----------------------------------------------------------------------------


# Main

if __name__ == '__main__':

    print "start: ", datetime.datetime.now()

    do_gaussian = True
    do_gaussian = False

    do_diffraction = False
    do_diffraction = True

    do_refraction = True
    do_refraction = False

    apply_inner_refraction = True
    apply_inner_refraction = False

    apply_outer_refraction = True
    apply_outer_refraction = False

    # Gaussian parameters
    phiblob = -1.		# < 0 ==> diverging
    phiblob = -3.		# < 0 ==> diverging
    phiblob =  0.3		# > 0 ==> converging
    phiblob = -20.		# < 0 ==> diverging
    phiblob = -30.		# < 0 ==> diverging
    phiblob = -10.		# < 0 ==> diverging
    phiblob = -0.3		# < 0 ==> diverging
    phiblob = 0.

    wex = 1.
    wey = 2. 

    wex = 0.3
    wey = 0.15

    wex = 2
    wey = 1

    # Diffraction parameters
    si_kol_d = 11./3.
    phiF_d = 0.01
    phiF_d = 0.025
    phiF_d = 2.
    phiF_d = 5.
    phiF_d = 10.

    # Refraction parameters
    si_kol_r = 8./3.
    si_kol_r = 5./3.
    si_kol_r = 11./3.
    phiF_r = 1. 
    phiF_r = 0.75 
    phiF_r = 0.4
    phiF_r = 0.05 
    phiF_r = 0.01 
    phiF_r = 0.2
    phiF_r = 0.1
    phiF_r = 10.
    phiF_r = 0.7
    phiF_r = 1.
    phiF_r = 0.5 
    phiF_r = 0.3
    phiF_r = 2. 
    phiF_r = 5.
    phiF_r = 10.
    phiF_r = 1.

    rF = 1. 			# Fresnel scale by definition

    # Note that array sizes scale as xwidth x ywidth
    xwidth = ywidth = 5.
    xwidth = ywidth = 10.
    xwidth = ywidth = 20.
    xwidth = ywidth = 30.
    xwidth = ywidth = 200.
    xwidth = ywidth = 100.
    xwidth = ywidth = 75.
    xwidth = ywidth = 50.

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

    print "diffractive phase: ", datetime.datetime.now()

    # Diffractive phase:
    xvec, yvec, xseries, phid, qxvec, qyvec, qshape = \
         genps.gen_pl_phase_screen_2d(si_kol_d,phiF_d,rF,inners,outers,\
         dx,dy,xwidth,ywidth,normfres=True)

    # subtract mean phase from screen (it doesn't really matter)
    phid -= phid.mean()

    nqx, nqy = shape(qshape)

    gausstype = 'diverge'
    if phiblob > 0.:  gausstype='converge'
    plotfilepre = 'kdipw_riss_2016_2018' + \
       '%4.2f'%(si_kol_d) + '_' + '%4.2f'%(si_kol_r) + \
       '_'+str(nqx)+'x'+str(nqy) + \
       '_rmsd_' + str(phiF_d*do_diffraction) + \
       '_rmsr_' + str(phiF_r*do_refraction) + \
       '_g_' + gausstype + '_' + str(abs(phiblob)*do_gaussian) 

    # get refraction phase screen:
    # define refraction scale using Fresnel phase for refracting screen
    # since we want to control refraction separately from diffraction

    #lr = rF * max(1., phiF_r)
    lscale_index = 2. / (si_kol_r-2.)
    u = max(1., phiF_r**lscale_index)
    ld = rF / u
    lr = rF * u

    #*** TEMPORARY ***
    inners = 5.*lr
    inners = 1.*lr
    inners = 3.*lr
    inners = 2.*lr
    #*** END TEMPORARY ***
    inners = lr
    outers = xmax/2.

    # phase on the lr scale:
    phi_lr = phiF_r * u**(5./6.)		# from phase SF 
    # gain on lr scale:
    g_lr = phi_lr / lr**2			# rough estimate of second derivative
    G_lr = sqrt(2.*g_lr**2 + g_lr**4)	# guess total approximate gain rms
    print "Phase on refraction scale = ", phi_lr
    print "1D gain on refraction scale = ", g_lr
    print "Total gain on refraction scale = ", G_lr

    print "refractive phase: ", datetime.datetime.now()

    xvec, yvec, xseries, phir, qxvec, qyvec, qshape = \
       genps.gen_pl_phase_screen_2d(si_kol_r,phiF_r,rF,inners,outers, \
       dx,dy,xwidth,ywidth, normfres=True, \
       apply_outer=apply_outer_refraction, \
       apply_inner=apply_inner_refraction)

    # subtract mean refraction phase
    phir -= phir.mean()

    # screen x,y array:
    xy = meshgrid(xvec, yvec)
    rsqvec = xy[0]**2 + xy[1]**2

    # Gaussian blob at center:
    gb = genps.gen_gaussian_lens(phiblob, 0., 0., wex, wey, dx, dy, xwidth, ywidth)

    xvecobs = copy(xvec)
    xvecobshalf = xvec[where(abs(xvec) < xmax/2.)]

    print "propagation: ", datetime.datetime.now()

    # propagator to observer: 
    # note definition of Fresnel scale implied is r_F^2 = \lambda D / 2\pi
    #kernel0 = exp(1j*pi*rsqvec/rF**2)
    kernel0 = exp(1j*rsqvec/(2.*rF**2))
    #kernel0 = exp(rsqvec/(2.*rF**2))

    # screen:
    phi = zeros(shape(phid))
    if do_gaussian:
       phi += gb
    if do_refraction:
       phi += phir
    if do_diffraction:
       phi += phid

    screen0 = exp(1j*phi)
    #screen0 = exp(phi)

    kernel0fft = fft2(kernel0)
    screen0fft= fft2(screen0)
    field0fft = kernel0fft * screen0fft
    field0 = ((dx*dy)/(2.*pi*rF**2)) * fftshift(ifft2(field0fft))
    intensity0 = abs(field0)**2 
    #intensity0 = abs(field0**2) 
    # nb: get same intensity0 using above two methods; kind of surprising
    obsphase = arctan2(field0.imag, field0.real)
    image = fftshift(abs(fft2(field0))) / shape(field0)[0] / pi  # norm'n = guess! 
    image /= image.max()
    intensity0fft = fft2(intensity0)
    acfI = fftshift(ifft2(abs(intensity0fft)**2)) / size(intensity0)
    acfI = acfI.real
    acfI -= 1.				# removes acf asymptotic value

    """
    ######################
    # also try with zero padding in both dimensions (==> 2x2 increase in size)
    # doesn't seem to get rid of some of the oscillations; are they artifacts
    # or real (looked at case with gaussian lens only with phiblob=-30
    # and xwidth = ywidth = 50; lens has 1/e widths (1,2)
    # in fact intensity2 looks wrong (too much ringing, even if 
    # a null screen is used.   Also, the unpadded result intensity0
    # doesn't give unit intensity when a null screen is used. 
    # it gives 0.97

    kernel2 = zeros((2*xvec.size, 2*yvec.size), dtype='complex')
    screen2 = zeros((2*xvec.size, 2*yvec.size), dtype='complex')

    kernel2[:kernel0.shape[0], :kernel0.shape[1]] = kernel0
    screen2[:screen0.shape[0], :screen0.shape[1]] = screen0

    kernel2fft = fft2(kernel2)
    screen2fft= fft2(screen2)
    field2fft = kernel2fft * screen2fft
    field2 = ((dx*dy)/(2.*pi*rF**2)) * ifft2(field2fft)
    #field2 = ((dx*dy)/(2.*pi*rF**2)) * fftshift(ifft2(field2fft))
    intensity2 = abs(field2)**2 
    obsphase2 = arctan2(field2.imag, field2.real)
    image2 = fftshift(abs(fft2(field2))) / shape(field2)[0] / pi  # norm'n = guess! 
    image2 /= image2.max()
    intensity2fft = fft2(intensity2)
    acfI2 = fftshift(ifft2(abs(intensity2fft)**2)) / size(intensity2)
    acfI2 = acfI.real
    acfI2 -= 1.				# removes acf asymptotic value
    ######################
    """

    print " "
    print "RMS Fresnel phase: = ", phiF_r

    lsmooth = lr / (2.*dx)		# smoothing length = one sigma
    intensity0_smoothed = ndimage.gaussian_filter(intensity0, sigma=(lsmooth, lsmooth), order=0, mode='wrap') 

    print " "
    print "Raw intensity:"
    print "     Mean intensity = ", intensity0.mean()
    print "      RMS intensity = ", intensity0.std()
    print "     Modulation index = ", intensity0.std() / intensity0.mean()

    print " "
    print "Smoothing scale = ", lsmooth, " samples"
    print "Smoothed intensity"
    print "     Mean intensity = ", intensity0_smoothed.mean()
    print "      RMS intensity = ", intensity0_smoothed.std()
    print "     Modulation index = ", intensity0_smoothed.std() / intensity0_smoothed.mean()
    print " "

    # calculate spectrum along one axis
    Ispatial = intensity0 - intensity0.mean()
    Ispectrum_slices = abs(fft(Ispatial))**2 / size(Ispatial)**2
    Ispectrum = average(Ispectrum_slices, axis=0)


    # This approach is to bin the 2D spectrum into scalar wavenumber bins
    # Doesn't work as well as the above approach that averages
    # 1D spectra and applies a q^2 weight when plotting. 

    Ispectrum2d = abs(fft2(Ispatial))**2 / size(Ispatial)**2

    # calculate binned spectrum 
    nqxhalf = shape(Ispatial)[0]/2
    nqyhalf = shape(Ispatial)[1]/2

    # long array:
    Ispectrum_scalar = zeros(nqxhalf*nqyhalf)
    qscalarvec = zeros(nqxhalf*nqyhalf)

    for i, qx in enumerate(qxvec[0:nqxhalf]):
        for j, qy in enumerate(qyvec[0:nqyhalf]):
            indxy = j*nqyhalf + i
            qscalar = sqrt(qx**2 + qy**2)
            qscalarvec[indxy] = qscalar
            Ispectrum_scalar[indxy] = Ispectrum2d[i,j]
                    
    qbins = copy(qxvec[0:nqxhalf])
    inds = digitize(qscalarvec, qbins)
    Ispectrum_binned = zeros(inds.max())

    for i in range(inds.max()):
       Ispectrum_binned[i]  = sum(Ispectrum_scalar[where(inds==i)])
    #"""
    nr, radial_prof = azimuthalAverage(fftshift(Ispectrum2d))

    # -----------------------------------------------------------------------------

    print "plotting: ", datetime.datetime.now()

    plot_wavenumber_spectrum()

    plot_fourframeI()

    plot_fourframeII()	    # screen phase, intensity, image, intensity ACF

    plot_Ihistogram()       # histogram of spatial intensity 

    plot_intensity_slice()

    plot_wavefield_phase()

    plot_1D_intensity_spectrum()

    plot_spectral_shapes()


    raw_input('hit return')
    close('all')
