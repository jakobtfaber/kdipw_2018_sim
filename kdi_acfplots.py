import scipy.ndimage as  ndimage
from scipy.stats import skew, kurtosis

#rcParams['font.size'] = 11

import datetime

import numpy as np

from preamble import *
from scipy.special  import fresnel

import matplotlib.mlab as mlab

import seaborn as sns

#------------------------------------------------------------------------------

def gen_pl_phase_screen_2d(si, phiF, rF, inner, outer, dx, dy, xwidth, ywidth, normfres=True):
    
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

    # 2016 Jan 1: put in upper wavenumber cutoff at array size
    # to avoid aliasing
    qmax = qxvec.max()/2.
    for i, qxi in enumerate(qxvec):
        for j, qyj in enumerate(qyvec):
            qsq = qxi**2 + qyj**2
            qshape[i,j] = (qout**2 + qsq)**(-si/4.) * exp(-qsq/(2.*qmax**2))
            qshape_rolloff[i,j] = exp(-qsq / (2.*qin**2))

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

    return xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape


def plot_fourframeII():

    """
    Plotting: four frame II:
    screen phase
    intensity pattern
    image
    intensity ACF
    """

    fig = plt.figure()
    ax1 = fig.add_subplot(221)

    """subplot(121)
    imshow(phi, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    xlabel(r'$X$')
    ylabel(r'$Y$')
    #tick_params(axis='x', labelbottom='off')

    #title(r'$\rm Screen \ \  {\cal R} \{ e^{i %s}\}$'%(phi_string))
    title(r'$\rm Screen\  Phase \ \ \phi$')
    tick_params(axis='both', labelsize='11')
    colorbar()"""
    
    #imshow(image, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    imshow(image, aspect='auto', origin='lower', interpolation='nearest')
    ax1.set_xlabel(r'$\theta_x$')
    ax1.set_ylabel(r'$\theta_y$')
    plt.title(r'$\rm Image \ \ I(\theta)$')
    #tick_params(axis='both', labelsize='11')
    #ax1.colorbar()
    
    """imshow(log10(intensity0), aspect='auto', origin='lower',extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    xlabel(r'$X$')
    title(r'$\rm log10 \ (Observed \ Intensity)$')
    tick_params(axis='both', labelsize='11')
    colorbar()"""
    
    ax2 = fig.add_subplot(222)
    
    #plot(imx, imy)
    #imshow(image, aspect='auto', origin='lower', extent=(-xmax, xmax, -ymax, ymax), interpolation='nearest')
    #sns.distplot(imx1, kde_kws={"color": "r", "lw": 1.5, "label": "Slice 1"}, hist=False)
    #sns.distplot(imx2, kde_kws={"color": "m", "lw": 1.5, "label": "Slice 2"}, hist=False)
    #sns.distplot(imx3, kde_kws={"color": "b", "lw": 1.5, "label": "Slice 3"}, hist=False)
    plt.plot(imx1, lw = 1.5, label = "Slice 1")
    #sns.distplot(sig, kde_kws={"color": "k", "lw": 1.5, "label": "SD"}, hist=False)
    #plt.scatter(sig,imx)
    #sns.distplot(negsig, kde_kws={"color": "b", "lw": 1.5, "label": "KDE 1"}, hist=False)
    #hist(imx, bins = 5000)
    #plt.plot(sigma, -(sigma))
    #mean = np.mean(imxx)
    #plot(imx, mlab.normpdf(imx, mean, sigma))
    ax2.set_ylabel(r'$\rm Image \ \ I(\theta)$')
    ax2.set_xlabel(r'$\theta_x$')
    #ax1.set_xlim([0,1e-8])
    #ax1.set_ylim([10e-10,10e8])
    #ax1.set_yscale("log", nonposy='clip')
    #tick_params(axis='x', labelbottom='off')
    plt.title(r'$\rm Image \ \ I(\theta) \ \ Slice$')
    #tick_params(axis='both', labelsize='11')
    #colorbar()

    ax3 = fig.add_subplot(223)

    sns.distplot(acf1, kde_kws={"color": "r", "lw": 1.5, "label": "Slice 1"}, hist=False)
    sns.distplot(acf2, kde_kws={"color": "m", "lw": 1.5, "label": "Slice 2"}, hist=False)
    sns.distplot(acf3, kde_kws={"color": "b", "lw": 1.5, "label": "Slice 3"}, hist=False)
    ax3.set_xlabel(r'$\rm \delta\theta_x $')
    ax3.set_ylabel(r'$\rm R_I (\delta\theta) $')
    #ax3.set_ylim([0,25])
    #ax3.set_yscale("log", nonposy='clip')
    #ax3.set_ylim([0,10**8])
    plt.title(r'$\rm R_I (\delta\theta) \ ACF \ \ Slice$')

    ax4 = fig.add_subplot(224)
    
    sns.distplot(intensity01, kde_kws={"color": "r", "lw": 1.5, "label": "Slice 1 (Intensity0)"}, hist=False)
    sns.distplot(intensity02, kde_kws={"color": "m", "lw": 1.5, "label": "Slice 2 (Intensity0)"}, hist=False)
    sns.distplot(intensity03, kde_kws={"color": "b", "lw": 1.5, "label": "Slice 3 (Intensity0)"}, hist=False)
    #sns.distplot(field001, kde_kws={"color": "m", "lw": 1.5, "label": "Slice 2 (Field00)"}, hist=False)
    #ax4.set_xlim([-0.5,10])
    ax4.set_ylabel(r'$\rm \Gamma(\delta X)$')
    ax4.set_xlabel(r'$\rm \delta X $')
    plt.title(r'$\rm \Gamma(\delta X) \ ACF \ \ Slice$')

    """subplot(122)
    ylim(0, 0.02)
    #sns.distplot(imxx) #, rug=True, hist=False)
    sns.distplot(imxx, kde_kws={"color": "k", "lw": 1, "label": "KDE"}) #, hist_kws={"histtype": "step", "linewidth": 2,"alpha": 1, "color": "g"})
    #ylim(0,10e7)
    #xlim(0,5e-8)"""
    
    """
    # zoom in to 1/4 of frame
    zoom = 16
    nxselect = where(abs(xvec)<=xmax/zoom)[0]    # where gives  tuple
    nyselect = where(abs(yvec)<=ymax/zoom)[0]
    nzoom = size(nxselect)
    acfIplot = zeros((nzoom, nzoom))
    for nzx, nx in enumerate(nxselect):
        for nzy, ny in enumerate(nyselect):
            acfIplot[nzx,nzy] = acfI[nx,ny]
    #acfIplot = acfI[where(abs(xvec)<=xmax/zoom), where(abs(xvec)<=xmax/zoom)]
    imshow(acfIplot,aspect='auto',origin='lower',extent=(-xmax/zoom,xmax/zoom,-ymax/zoom,ymax/zoom), interpolation='nearest')
    xlabel(r'$\rm Lag \ \ \delta X$')
    ylabel(r'$\rm Lag \ \ \delta Y$')

    title(r'$\rm \delta(Intensity) \ ACF $')
    tick_params(axis='both', labelsize='11')
    colorbar()

    #annotate(r'$\rm {\bf Diffraction+Refraction}: \ \phi_F = %5.2f\ rad \ \ \ {\bf Gaussian \ Lens}: \ \phi_G = %5.2f \ rad $'%(float(do_refraction), phiblob*float(do_gaussian)) , xy=(0.5, 0.95), xycoords='figure fraction', ha='center', fontsize=12)
    """
    show()
    #savefig(plotfilepre + '_fourframesII' + '.pdf')
    #raw_input('hit return')
    #close()
    return

if __name__ == '__main__':

    print "start: ", datetime.datetime.now()
    
    # Diffraction parameters
    si_kol_d = 11./3.
    phiF_d = 0.01
    phiF_d = 0.025
    phiF_d = 2.
    #phiF_d = 5.
    #phiF_d = 10.
    
    
    rF = 1.             # Fresnel scale by definition
    
    # Note that array sizes scale as xwidth x ywidth
    #xwidth = ywidth = 5.
    #xwidth = ywidth = 10.
    xwidth = ywidth = 20.
    #xwidth = ywidth = 30.
    #xwidth = ywidth = 200.
    #xwidth = ywidth = 100.
    #xwidth = ywidth = 75.
    #xwidth = ywidth = 50.
    
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
        gen_pl_phase_screen_2d(si_kol_d,phiF_d,rF,inners,outers,\
                                     dx,dy,xwidth,ywidth,normfres=True)
    
    xy = meshgrid(xvec, yvec)
    rsqvec = xy[0]**2 + xy[1]**2

    
    xvecobs = copy(xvec)
    xvecobshalf = xvec[where(abs(xvec) < xmax/2.)]

    kernel0 = exp(1j*rsqvec/(2.*rF**2))
    
    # screen:
    phi = zeros(shape(phid))
    phi += phid
    
    screen0 = exp(1j*phi)
    
    #kernel0fft = fft2(kernel0)
    kernel0fft = fft2(kernel0)
    #screen0fft= fft2(screen0)
    screen0fft= fft2(screen0)
    field0fft = kernel0fft * screen0fft
    #field0 = field0fft
    #field0 = ((dx*dy)/(2.*pi*rF**2)) * fftshift(ifft2(field0fft))
    field0 = ((dx*dy)/(2.*pi*rF**2)) * fftshift(ifft2(field0fft))
    field00 = abs(field0fft)**2
    field00 = ifft2(field00)
    field00_flat = field00.flatten()
    field001 = field00[1,:]
    #field002 = field00[2,:]
    #intensity0 = abs(field0)**2
    intensity0 = abs(field0)**2
    intensity01 = intensity0[1,:]
    intensity02 = intensity0[2,:]
    intensity03 = intensity0[3,:]
    intensity0_flat = intensity0.flatten()
    
    #print "intensity: ", intensity0
    
    #intensity0 = abs(field0**2)
    # nb: get same intensity0 using above two methods; kind of surprising
    obsphase = arctan2(field0.imag, field0.real)
    #image = fftshift(abs(fft2(field0))) / shape(field0)[0] / pi  # norm'n = guess!
    image = fftshift(abs(fft2(field0))) / shape(field0)[0] / pi
    image /= image.max()
    #image = image[0,:]
    imx = image.flatten()
    imx1 = image[1,:]
    imx2 = image[2,:]
    imx3 = image[3,:]
    #imxx = imxx.flatten()
    #imx = image[0,:]
    imy = image[:,0]
    sigma = np.std(imx)
    sig = [x + sigma for x in imx]
    sig = np.asarray(sig)

    
    #intensity0fft = fft2(intensity0)
    intensity0fft = fft2(intensity0)
    #acfI = fftshift(ifft2(abs(intensity0fft)**2)) / size(intensity0)
    acfI = fftshift(ifft2(abs(intensity0fft)**2)) / size(intensity0)
    acf_flat = acfI.flatten()
    acf1 = acfI[1,:]
    acf2 = acfI[2,:]
    acf3 = acfI[50,:]
    acfI = acfI.real
    acfI -= 1.                # removes acf asymptotic value

    print "plotting: ", datetime.datetime.now()
    
    
    #print(imy)
    #print(sig)
    

    plot_fourframeII()
