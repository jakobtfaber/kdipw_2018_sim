import numpy as np

#
# wavefront definitions
#

def wavefront_initialize(pixelsize_h=1e-6,pixelsize_v=1e-6,npixels_h=1024,npixels_v=1024,amplitude_value=0.0):
    #
    #create array at object (aperture) plane
    #
    amplitude = np.zeros((npixels_h,npixels_v))  # amplitude map
    amplitude += amplitude_value
    p_i_h = np.arange(npixels_h) * pixelsize_h
    p_x = (p_i_h - 0.5 * (p_i_h[-1] - p_i_h[0]) )
    p_i_v = np.arange(npixels_v) * pixelsize_v
    p_y = (p_i_v - 0.5 * (p_i_v[-1] - p_i_v[0]) )
    return p_x,p_y,amplitude

def wavefront_aperture(p_x,p_y,amplitude,diameter=40e-6,type=0):
    # aperture_type: 0=circular, 1=Square, 2=Gaussian
    p_xx = p_x[:, np.newaxis]
    p_yy = p_y[np.newaxis, :]

    filter = np.zeros_like(amplitude)
    if type == 0:  # Circular aperture
        radius = (diameter/2)
        print("radius=%f um"%(1e6*radius))
        filter_illuminated_indices = np.where(p_xx**2 + p_yy**2 < radius**2)
        if filter_illuminated_indices[0].size ==0:
            print("Warning: wavefront_aperture(): Nothing goes trough the aperture")
        else:
            filter[filter_illuminated_indices] = 1.0
    elif type == 1:  # square
        radius = (diameter/2)
        print("radius=%f um"%(1e6*radius))
        filter_illuminated_indices = np.where( (np.abs(p_xx) < radius) & (np.abs(p_yy) < radius))
        if filter_illuminated_indices[0].size ==0:
            print("Warning: wavefront_aperture(): Nothing goes trough the aperture")
        else:
            filter[filter_illuminated_indices] = 1.0
    elif type == 2:  # Gaussian
        sigma = diameter/2.35
        print("source sigma=%f um"%(1e6*sigma))
        rho2 = p_xx**2 + p_yy**2
        #TODO: add Gaussian amplitude
        filter = np.sqrt(np.exp(-rho2/2/sigma**2)) # Gaussian in intensity, so srrt for amplitude
        filter = np.exp(-rho2/2/sigma**2) # Gaussian amplitude
    else:
        raise ValueError("Aperture type (shape) not valid")

    return p_x,p_y,amplitude*filter

#
# tools
#
def propagator2d(x,y,z,method="fraunhofer",wavelength=1e-10,propagation_distance=1.0,return_angles=0):
    #
    # interface to different propagators
    #
    from timeit import default_timer as timer

    t_start = timer()
    if method == "fraunhofer":
        x1,y1,z1 = propagator2d_fraunhoffer(x,y,z,wavelength=wavelength)
        if return_angles:
            pass
        else:
            x1 *= propagation_distance
            y1 *= propagation_distance
    elif method == "fourier_convolution":
        x1,y1,z1 = propagator2d_fourier_convolution(x,y,z,propagation_distance=propagation_distance,wavelength=wavelength)
        if return_angles:
            x1 /= propagation_distance
            y1 /= propagation_distance
    elif method == "integral":
        x1,y1,z1 = propagator2d_integral(x,y,z,propagation_distance=propagation_distance,wavelength=wavelength)
        if return_angles:
            x1 /= propagation_distance
            y1 /= propagation_distance
    elif method == "srw":
        x1,y1,z1 = propagator2d_srw(x,y,z,propagation_distance=propagation_distance,wavelength=wavelength)
        if return_angles:
            x1 /= propagation_distance
            y1 /= propagation_distance
    else:
        raise Exception("method %s not implemented"%method)
    t_end = timer()
    print("Elapsed time in propagation calculations: %5.3f ms"%((t_end-t_start)*1e3))
    print("Shapes in propagation calculations: before: ",z.shape," after: ",z1.shape)
    print("Limits in propagation calculations H: before: ",x[0],x[-1]," after: ",x1[0],x1[-1]," points: ",x.shape)
    print("Limits in propagation calculations V: before: ",y[0],y[-1]," after: ",y1[0],y1[-1]," points: ",y.shape)
    return x1,y1,z1

def propagator2d_srw(p_x,p_y,amplitude,propagation_distance=1.0,wavelength=1e-10):
    #
    # convolving with the Fresnel kernel via SRW package
    #
    import srwlib
    from NumpyToSRW import numpyArrayToSRWArray, SRWWavefrontFromElectricField, SRWEFieldAsNumpy

    # srw_amplituder = numpyArrayToSRWArray(amplitude)
    # print(type(srw_amplituder))



    srw_wfr = SRWWavefrontFromElectricField(p_x[0], p_x[-1], amplitude,
                                  p_y[0], p_y[-1], np.zeros_like(amplitude),
                                  12396.0/(wavelength*1e10), 1.0, 1.0, 1e-3, 1.0, 1e-3)

    print(type(srw_wfr))


    #
    # propagation
    #
    optDrift = srwlib.SRWLOptD(propagation_distance) #Drift space

    #                 0  1  2   3  4  5   6   7   8   9 10 11
    # propagParDrift = [1, 1, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]
    propagParDrift = [0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]
    #Wavefront Propagation Parameters:
    #[0]: Auto-Resize (1) or not (0) Before propagation
    #[1]: Auto-Resize (1) or not (0) After propagation
    #[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
    #[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[5]: Horizontal Range modification factor at Resizing (1. means no modification)
    #[6]: Horizontal Resolution modification factor at Resizing
    #[7]: Vertical Range modification factor at Resizing
    #[8]: Vertical Resolution modification factor at Resizing
    #[9]: Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
    optBL = srwlib.SRWLOptC([optDrift], [propagParDrift]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)

    print('   Simulating Electric Field Wavefront Propagation bu SRW ... ', end='\n')
    srwlib.srwl.PropagElecField(srw_wfr, optBL)


    amplitude2 = SRWEFieldAsNumpy(srw_wfr)
    amplitude2 = amplitude2[0,:,:,0]
    print("Amplitude shape before:",amplitude.shape,"; after: ",amplitude2.shape)



    # fft = np.fft.fft2(amplitude)
    #
    # # frequency for axis 1
    # pixelsize = p_x[1] - p_x[0]
    # npixels = p_x.size
    # freq_nyquist = 0.5/pixelsize
    # freq_n = np.linspace(-1.0,1.0,npixels)
    # freq_x = freq_n * freq_nyquist
    # # freq = freq * wavelength
    #
    # # frequency for axis 2
    # pixelsize = p_y[1] - p_y[0]
    # npixels = p_y.size
    # freq_nyquist = 0.5/pixelsize
    # freq_n = np.linspace(-1.0,1.0,npixels)
    # freq_y = freq_n * freq_nyquist
    # # freq_y = freq_y * wavelength
    #
    # freq_xy = np.array(np.meshgrid(freq_y,freq_x))
    #
    # fft *= np.exp((-1.0j) * np.pi * wavelength * propagation_distance *
    #               np.fft.fftshift(freq_xy[0]*freq_xy[0] + freq_xy[1]*freq_xy[1]) )
    # ifft = np.fft.ifft2(fft)

    p_x2 = np.linspace(srw_wfr.mesh.xStart, srw_wfr.mesh.xFin, srw_wfr.mesh.nx)
    p_y2 = np.linspace(srw_wfr.mesh.yStart, srw_wfr.mesh.yFin, srw_wfr.mesh.ny)

    return p_x2,p_y2,amplitude2


def propagator2d_integral(p_x,p_y,amplitude,propagation_distance=1.0,wavelength=1e-10,shuffle_interval=1e-5):
    #
    # Fresnel-Kirchhoff integral (neglecting inclination factor)
    #
    det_x = p_x.copy()
    det_y = p_y.copy()

    #
    # manual
    #
    # p_xy = np.zeros((2,p_x.size,p_y.size))
    # det_xy = np.zeros((2,det_x.size,det_y.size))
    # for i in range(p_x.size):
    #     for j in range(p_y.size):
    #         p_xy[0,i,j] = p_x[i]
    #         p_xy[1,i,j] = p_y[j]
    #         det_xy[0,i,j] = det_x[i]
    #         det_xy[1,i,j] = det_y[j]

    #
    # np
    #
    p_xy = np.array(np.meshgrid(p_y,p_x))
    det_xy = np.array(np.meshgrid(det_y,det_x))


    amplitude_propagated = np.zeros_like(amplitude,dtype='complex')

    wavenumber = 2 * np.pi / wavelength

    for i in range(det_x.size):
        for j in range(det_y.size):
            if shuffle_interval == 0:
                rd_x = 0.0
                rd_y = 0.0
            else:
                rd_x = (np.random.rand(p_x.size,p_y.size)-0.5)*shuffle_interval
                rd_y = (np.random.rand(p_x.size,p_y.size)-0.5)*shuffle_interval

            r = np.sqrt(    np.power(p_xy[0,:,:] + rd_x - det_xy[0,i,j],2) +
                            np.power(p_xy[1,:,:] + rd_y - det_xy[1,i,j],2) +
                            np.power(propagation_distance,2) )
            amplitude_propagated[i,j] = (amplitude / r * np.exp(1.j * wavenumber *  r)).sum()

    return det_x,det_y,amplitude_propagated # .reshape((det_x.size,det_y.size))

def propagator2d_fourier_convolution(p_x,p_y,image,propagation_distance=1.0,wavelength=1e-10):
    #
    # convolving with the Fresnel kernel via FFT multiplication
    #
    fft = np.fft.fft2(image)

    # frequency for axis 1
    pixelsize = p_x[1] - p_x[0]
    npixels = p_x.size
    freq_nyquist = 0.5/pixelsize
    freq_n = np.linspace(-1.0,1.0,npixels)
    freq_x = freq_n * freq_nyquist
    # freq = freq * wavelength

    # frequency for axis 2
    pixelsize = p_y[1] - p_y[0]
    npixels = p_y.size
    freq_nyquist = 0.5/pixelsize
    freq_n = np.linspace(-1.0,1.0,npixels)
    freq_y = freq_n * freq_nyquist
    # freq_y = freq_y * wavelength

    freq_xy = np.array(np.meshgrid(freq_y,freq_x))

    fft *= np.exp((-1.0j) * np.pi * wavelength * propagation_distance *
                  np.fft.fftshift(freq_xy[0]*freq_xy[0] + freq_xy[1]*freq_xy[1]) )

    # fft = np.fft.fftshift(fft)
    # fft *= np.exp((-1.0j) * np.pi * wavelength * propagation_distance *
    #               (freq_xy[0]*freq_xy[0] + freq_xy[1]*freq_xy[1]) )
    # fft = np.fft.ifftshift(fft)

    ifft = np.fft.ifft2(fft)

    return p_x.copy(),p_y.copy(),ifft


def propagator2d_fraunhoffer(p_x,p_y,image,wavelength=1e-10):
    """
    Fraunhoffer propagator
    :param x: x array of spatial coordinates in meters
    :param y: y array of spatial coordinates in meters
    :param complax_amplitude array: shape: [n_points_x,n_points_y]
    :param wavelength: photon wavelength in meters
    :return: three arrays with the propagated pattern : angle_x [rad], angle_y, complex_amplitude.
    """

    #
    #compute Fourier transform
    #
    F1 = np.fft.fft2(image)  # Take the fourier transform of the image.
    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2 = np.fft.fftshift( F1 )

    # frequency for axis 1
    pixelsize = p_x[1] - p_x[0]
    npixels = p_x.size
    freq_nyquist = 0.5/pixelsize
    freq_n = np.linspace(-1.0,1.0,npixels)
    freq_x = freq_n * freq_nyquist
    freq_x *= wavelength

    # frequency for axis 2
    pixelsize = p_y[1] - p_y[0]
    npixels = p_y.size
    freq_nyquist = 0.5/pixelsize
    freq_n = np.linspace(-1.0,1.0,npixels)
    freq_y = freq_n * freq_nyquist
    freq_y *= wavelength

    return freq_x,freq_y,F2

def line_image(image,horizontal_or_vertical='H':
    if horizontal_or_vertical == "H":
        npixels = image.shape[0]
        tmp = image[:,image.shape[(1)]/(2)]
    else:
        npixels = image.shape[1]
        tmp = image[image.shape[0]/2,:]
    return tmp

def line_fwhm(line):
    #
    #CALCULATE fwhm in number of abscissas bins (supposed on a regular grid)
    #
    tt = np.where(line>=max(line)*0.5)
    if line[tt].size > 1:
        # binSize = x[1]-x[0]
        FWHM = (tt[0][-1]-tt[0][0])
        return FWHM
    else:
        return -1

#
# plotting tools
#
def plot_show():

    import matplotlib.pylab as plt

    plt.show()

def plot_image(mymode,theta,psi,title="TITLE",xtitle=r"X [$\mu m$]",ytitle=r"Y [$\mu m$]",cmap=None,show=1):

    import matplotlib.pylab as plt

    fig = plt.figure()

    # cmap = plt.cm.Greys
    plt.imshow(mymode.T,origin='lower',extent=[theta[0],theta[-1],psi[0],psi[-1]],cmap=cmap)
    plt.colorbar()
    ax = fig.gca()
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    plt.title(title)

    if show: plt.show()

def plot(*positional_parameters,title="",xtitle="",ytitle="",show=1,legend=None,color=None):

    import matplotlib.pylab as plt

    n_arguments = len(positional_parameters)
    if n_arguments == 0:
        return

    fig = plt.figure()

    if n_arguments == 1:
        y = positional_parameters[0]
        x = np.arange(y.size)
        plt.plot(x,y,label=legend)
    elif n_arguments == 2:
        x = positional_parameters[0]
        y = positional_parameters[1]
        plt.plot(x,y,label=legend,color=color)
    elif n_arguments == 4:
        x1 = positional_parameters[0]
        y1 = positional_parameters[1]
        x2 = positional_parameters[2]
        y2 = positional_parameters[3]
        if legend != None:
            legend1 = legend[0]
            legend2 = legend[1]
        else:
            legend1 = None
            legend2 = None
        if color != None:
            color1 = color[0]
            color2 = color[1]
        else:
            color1 = None
            color2 = None
        plt.plot(x1,y1,label=legend1,color=color1)
        plt.plot(x2,y2,label=legend2,color=color2)
    else:
        "Incorrect number of arguments, plotting only two first arguments"
        x = positional_parameters[0]
        y = positional_parameters[1]
        plt.plot(x,y,label=legend)

    if legend != None:
        ax = plt.subplot(111)
        ax.legend(bbox_to_anchor=(1.1, 1.05))

    plt.title(title)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)


    if show:
        plt.show()

def main():
    #
    # inputs (in SI)
    #

    wavelength        = 1.24e-10

    aperture_diameter = 40e-6 # if Gaussian, aperture_diameter = 2.35*sigma
    aperture_type     = 2     # 0=circular, 1=Square, 2=Gaussian (sigma = diameter/2.35)

    pixelsize_x = 1e-6
    pixelsize_y = pixelsize_x
    npixels_x =  1024    # 200 #
    npixels_y =  npixels_x # 50  #

    propagation_distance = 30.0 

    # method = "fourier_convolution"
    method = "fraunhofer"
    # method = "integral"
    # method = "srw"

    #
    # calculations
    #

    # get a wavefront
    p_x,p_y,amplitude = wavefront_initialize(pixelsize_x,pixelsize_y,npixels_x,npixels_y,amplitude_value=1.0)
    # set aperture
    p_x,p_y,amplitude = wavefront_aperture(p_x,p_y,amplitude,diameter=aperture_diameter,type=aperture_type)
    #plot aperture
    plot_image(np.abs(amplitude)**2,p_x*1e6,p_y*1e6, show=0,
               title="aperture intensity, Diameter=%5.1f um"%(1e6*aperture_diameter),xtitle="X [um]",ytitle="Y [um]")

    #
    # propagation
    #
    angle_x, angle_y, amplitude_propagated = propagator2d(p_x,p_y,amplitude,method=method,wavelength=wavelength,
                                        propagation_distance=propagation_distance,return_angles=1)

    # angle_x, angle_y, amplitude_propagated = propagator2d_fraunhoffer(p_x,p_y,amplitude,wavelength=wavelength)
    if method == "fraunhofer":
        print("Fraunhoffer diffraction valid for distances > > a^2/lambda = %f m"%((aperture_diameter/2)**2/wavelength))

    plot_image(np.abs(amplitude_propagated)**2,angle_x*1e6,angle_y*1e6, show=0,
               title="Diffracted intensity (%s)"%method,xtitle="X [urad]",ytitle="Y [urad]")

    #
    # extract profiles and calculate theoretical ones
    #

    # retrieve H and V profiles
    horizontal_intensity_profile = line_image(np.abs(amplitude_propagated)**2,horizontal_or_vertical='H')
    horizontal_intensity_profile /= horizontal_intensity_profile.max()

    vertical_intensity_profile = line_image(np.abs(amplitude_propagated)**2,horizontal_or_vertical='V')
    vertical_intensity_profile /= vertical_intensity_profile.max()

    # theoretical profile
    if aperture_type == 0: #circular, also display analytical values
        from scipy.special import jv
        x = (2*np.pi/wavelength) * (aperture_diameter/2) * angle_x
        y = (2*np.pi/wavelength) * (aperture_diameter/2) * angle_y
        U_vs_theta_x = 2*jv(1,x)/x
        U_vs_theta_y = 2*jv(1,y)/y
        I_vs_theta_x = U_vs_theta_x**2
        I_vs_theta_y = U_vs_theta_y**2
    elif aperture_type == 1: # square
        x = (2*np.pi/wavelength) * (aperture_diameter/2) * angle_x
        y = (2*np.pi/wavelength) * (aperture_diameter/2) * angle_y
        U_vs_theta_x = 2*np.sin(x)/x
        U_vs_theta_y = 2*np.sin(y)/y
        I_vs_theta_x = U_vs_theta_x**2
        I_vs_theta_y = U_vs_theta_y**2
        I_vs_theta_x /= I_vs_theta_x.max()
        I_vs_theta_y /= I_vs_theta_y.max()
    elif aperture_type == 2: #Gaussian
        sigma = aperture_diameter/2.35
        sigma_ft = 1.0 / sigma * wavelength / (2.0 * np.pi)
        # Factor 2.0 is because we wwant intensity (amplitude**2)
        I_vs_theta_x = np.exp( -2.0*(angle_x**2/sigma_ft**2/2) )
        I_vs_theta_y = np.exp( -2.0*(angle_y**2/sigma_ft**2/2) )

    fwhm_intensity_profile_horizontal = line_fwhm(horizontal_intensity_profile) * (angle_x[1]-angle_x[0])
    fwhm_intensity_profile_vertical = line_fwhm(vertical_intensity_profile) * (angle_y[1]-angle_y[0])
    fwhm_theoretical_profile_horizontal = line_fwhm(I_vs_theta_x) * (angle_x[1]-angle_x[0])
    fwhm_theoretical_profile_vertical = line_fwhm(I_vs_theta_y) * (angle_y[1]-angle_y[0])

    #
    # calculate widths
    #
    print("HORIZONTAL FWHM (%s) : %f urad, FWHM theoretical: %f urad, 1.22*wavelength/Diameter: %f urad"%(
        method,
        1e6*fwhm_intensity_profile_horizontal,1e6*fwhm_theoretical_profile_horizontal,1e6*1.22*wavelength/aperture_diameter))
    print("VERTICAL FWHM (%s) : %f urad, FWHM theoretical: %f urad, 1.22*wavelength/Diameter: %f urad"%(
        method,
        1e6*fwhm_intensity_profile_vertical,1e6*fwhm_theoretical_profile_vertical,1e6*1.22*wavelength/aperture_diameter))
    print("HORIZONTAL (4pi/lambda) sigma sigma' : (%s): %f, theoretical: %f "%(method,
        4*np.pi / wavelength * fwhm_intensity_profile_horizontal/2.35 * aperture_diameter/2.35,
        4*np.pi / wavelength * fwhm_theoretical_profile_horizontal/2.35 * aperture_diameter/2.35  ))



    # plot profiles
    plot( angle_x*1e6, horizontal_intensity_profile, angle_x*1e6, I_vs_theta_x, show=0,
          legend=["profile","theory"],color=["red","black"],
          title="Horizontal profile of diffracted intensity (%s)"%method,xtitle='theta [urad]',ytitle='Diffracted intensity [a.u.]')
    plot( angle_y*1e6, vertical_intensity_profile, angle_y*1e6, I_vs_theta_y, show=1,
          legend=["profile","theory"],color=["red","black"],
          title="Vertical profile of diffracted intensity (%s)"%method,xtitle='theta [urad]',ytitle='Diffracted intensity [a.u.]')

if __name__ == "__main__":
    main()
