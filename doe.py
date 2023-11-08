#!/usr/bin/python3
#
# Detection Of Extrema (doe)
# 
# Code designed to detect the position of one or multiple blended peaks in the input data
# either in emission i.e. cross-correlation function mode (-t max, default)
# either in absorption i.e. normalized spectrum mode (-t min)
#
# The standard output lies on one line with:
# - name of the input file
# - number of peaks detected
# - the positions of the peaks and the associated errors (in the same unit as the input data)
#
# Input:   *.dat     (ASCII file with 2 columns [abscissa ordinates] or FITS file (with CRVAL1 and CDELT1 keywords))
# Outputs: *_sd.dat  (selected successives derivatives [abscissa ordinates 1st 2nd 3rd])
#          *_xp.dat  (peaks properties [x x_err f width xwmin xwmax sigma it]
#                     if in addition, a fit of the components is performed, then:
#                        with -G option: [x x_err f width xwmin xwmax sigma it, Gaussian fit (mu, mu_err, sigma, n) + offset] 
#                        with -L option: [x x_err f width xwmin xwmax sigma it, Lorentzian fit (mu, mu_err, gamma, n) + offset]
#                        with -V option: [x x_err f width xwmin xwmax sigma it, Voigtian fit (mu, mu_err, alpha, gamma, n) + offset]
#                        with -R option: [x x_err f width xwmin xwmax sigma it, rotational fit (mu, mu_err, epsilon, vsini, n) + offset]
#          *.png     (if -p or -pp is given, plot of the CCF and its successive derivatives)
#
# History:
# 20231106: Set the main part of the code in a function run_doe to be used in import
# 20221106: Add the values of velocities, errors and widths on the plot
#         : Add the '-m' option to force the number of components to fit
#         : Add the internal LEASTSQ parameter to compare Gaussian fitting procedure of scipy.optimize
#         : Refine the plot
#         : in gaussian_filter1d was change from mode='reflect' to mode='nearest' => better behaviour for large sigma values  
#         : break back compatibility for fitting parmaeters: -g => -G, -l => -L, -vf => -V, -r => -R, --VELTOL => COMPTOL 
# 20200609: Missing 0.5 in the exponential of the Gaussian mixture model fit function gmm() (Thanks Kate!)
# 20200713: Add the bounds for the fit and the --velrange input parameter
# 20170608: original version

import sys
import numpy as np
import argparse as ap
import os.path as op
from scipy import ndimage  # for calculating the derivatives
from scipy import special  # for fitting a Voigt profile
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit, leastsq

def find_widths(x, xp, ggf):
   ''' Find width for a given peak using the second derivative
   '''
   z_sign = np.sign(ggf).astype(int)
   idx = np.where(np.diff(z_sign))[0]
   
   xw0 = x[idx]

   xwi = np.zeros(len(xw0))

   #Interpolation of the abscissae
   for ind, i in enumerate(idx):
      # interp need increasing abscissae
      if np.all(np.diff((ggf[i], ggf[i+1])) > 0):
         xwi[ind] = np.interp(0.0, (ggf[i], ggf[i+1]), (x[i], x[i+1]))
      else:
         xwi[ind] = np.interp(0.0, (ggf[i+1], ggf[i]), (x[i+1], x[i]))

   xw = []

   #Select only xwi around selected peaks
   for j in range(len(xp)):
      #for i in range(0, len(xw0)-1,2):
      for i in range(0, len(xw0)-1):
         if xp[j] > xwi[i] and xp[j] < xwi[i+1]:
            xw.append(xwi[i])
            xw.append(xwi[i+1])

   #xw = np.array(sorted(set(xw)))
   xw = np.array(xw)
   #print('xw=', xw)

   #Compute the width for each peak
   wp = np.zeros(len(xp))
   wp = xw[1::2] - xw[0::2]
   #print('wp=', wp)

   #print(len(xw0), len(xw), xw)

   #for ind, i in enumerate(xp):
   #   #dw = np.array(sorted(i-xw))
   #   dw = i - xw
   #   dw_plus = dw[ dw > 0. ]
   #   dw_minus = dw[ dw < 0. ]
   #   print(ind, dw_plus, dw_minus)
   #   if xw.size:
   #      wp[ind] = min(dw_plus) - max(dw_minus)

   return xw, wp

def find_peaks(x, f, ggf, gggf):
   ''' Detect positive sign change on the third derivative z.
       Return positions of the sign change on x, f and ggf 
   '''
   z_sign = np.sign(gggf).astype(int)

   # Example:
   #   -1 -1 -1 +1 +1 +1 -1 -1  np.roll(z_sign, 1)
   # -(-1 -1 +1 +1 +1 -1 -1 -1) z_sign
   #    0  0 -2  0  0 +2  0  0
   # -2 detects the positive sign change while 
   # +2 detects the negative sign change

   zeros = ((np.roll(z_sign, 1) - z_sign) == -2)

   zeros[0] = 0
   msk = zeros == 1
   xp = x[msk]
   fp = f[msk]
   ggfp = ggf[msk]

   idx = np.where(msk)[0]

   # Interpolation of the values 
   xpi = np.zeros(len(idx))
   fpi = np.zeros(len(idx))
   ggfpi = np.zeros(len(idx))

   for ind, i in enumerate(idx):
      xpi[ind] = np.interp(0., (gggf[i-1], gggf[i]), (x[i-1], x[i]))
      fpi[ind] = np.interp(xpi[ind], (x[i-1], x[i]), (f[i-1], f[i]))
      ggfpi[ind] = np.interp(xpi[ind], (x[i-1], x[i]), (ggf[i-1], ggf[i]))

   if False:
      print("fp =", fp)
      print("fpi =", fpi)
      print("ggfp =", ggfp)
      print("ggfpi =", ggfpi)
      print("xp = ", xp)
      print("xpi = ", xpi)

   return xpi, fpi, ggfpi 

def gmm(x, *params):
    '''Gaussian mixture model 
       First element of params is the offset
       Block of three elements gives mu (mean), sig (sigma), n (normalisation factor)
    '''
    y = np.zeros_like(x)
    for i in range(1, len(params), 3):
        mu = params[i]
        sig = params[i+1]
        n = params[i+2]
        y += n * np.exp( -0.5*((x - mu)/sig)**2)
    y += params[0]
    return y

def lmm(x, *params):
    '''Lorentzian mixture model 
       First element of params is the offset
       Block of three elements gives mu (mean), gam (gamma), n (normalisation factor)
    '''
    y = np.zeros_like(x)
    for i in range(1, len(params), 3):
        mu = params[i]
        gam = params[i+1]
        n = params[i+2]
        y += n * gam / np.pi / ((x - mu)**2 + gam**2)
    y += params[0]
    return y

def vmm(x, *params):
    '''Voigtian mixture model
       First element of params is the offset
       Block of four elements gives mu (mean), alp (alpha), gam (gamma), n (normalisation factor)
    '''    
    y = np.zeros_like(x)
    for i in range(1, len(params), 4):
        mu = params[i]
        alp = params[i+1]
        gam = params[i+2]
        n = params[i+3]
        sig = alp / np.sqrt(2 * np.log(2))
        y += n * np.real(special.wofz(((x - mu) + 1j * gam) / sig / np.sqrt(2)))
    y += params[0]
    return y

def rmm(x, *params):
    ''' Rotational mixture model
        First element of params is the offset (currently not used in the fitting process)
        Block of four elements gives mu (mean), eps (epsilon), vr (rotational velocity), n (normalisation factor)
    '''
    y = np.zeros_like(x)
    for i in range(1, len(params), 4):
        mu = params[i]
        eps = params[i+1]
        vr = params[i+2]
        n = params[i+3]
        xx = ((x-mu)/vr)**2
        #xx = xx.compress(xx < 1.)
        r0 = (2+(np.pi/2-2)*eps) / (np.pi*vr*(1-eps/3))
        xx[abs(xx)>1] = 1
        y += n/r0 * (2*(1-eps)*(1-xx)**0.5 + 0.5*np.pi*eps*(1-xx)) / (np.pi*vr*(1-eps/3.))
    y += params[0]
    return y

def load_fits(finp, extn, verb):
    ''' Load and extract data from a FITS file
        finp: input filename
        extn: extension number
        verb: verbosity
    '''
    from astropy.io import fits

    try:
        hdulist = fits.open(finp)
    except:
        print("Input filename does not exist or is not in FITS format.")
        quit(1)

    head = hdulist[extn].header
    y = np.array(hdulist[extn].data).flatten()
    hdulist.close()

    try: 
        crpix = head["CRPIX1"]
    except KeyError: 
        crpix = 1.
    
    try: 
        crval = head["CRVAL1"]
    except KeyError: 
        print("Error: CRVAL1 is undefined in "+finp)
        quit(1)
   
    try: 
        cdelt = head["CDELT1"]
    except KeyError: 
        if verb:
            print("Warning: CDELT1 is undefined in "+finp)
        try: 
            cdelt = head["CD1_1"] # ESO keyword used if CDELT1 is undefined
        except KeyError:
            print("Error: CDELT1 or CD1_1 is undefined in "+finp)
            quit(1)

    try: 
        naxis = head["NAXIS1"]
    except KeyError: 
        naxis = numpy.size(y)

    x_start = crval - (crpix - 1)*cdelt
    x_final = x_start + (naxis-1)*cdelt
    x = np.linspace(x_start, x_final, naxis)
    x_delta = x[1]-x[0]
    if abs(x_delta - cdelt) > 1e-8:
        print("Error: interval problem.")
        quit(1)

    try:
        if head['CTYPE1']=='log(wavelength)':
            x = np.exp(x)
    except KeyError:
        pass

    if verb:
        print("Number of points = ", naxis)
        print("Initial x        = ", x[0])
        print("Final x          = ", x[naxis-1])
        print("Step             = ", cdelt, x_delta)

    return x, y

#def method_err(SIGMA, dx, popt):
#  return abs(popt[0] + popt[1]*SIGMA + popt[2]*dx + popt[3]*SIGMA*dx)

def run_doe(IFN, BLA, PLT, TYP, ONE_PASS, THRES0, THRES2, SIGMA0, N, N_OFFSET, XMIN, XMAX, YMIN, YMAX, NOT, LAB, PURE_DER, EXT_NUMBER, Gaussian_FIT, Lorentzian_FIT, Voigtian_FIT, ROTATIONAL_FIT, COMPTOL, NO_OUTPUTS, NO_OUTPUT2, NC):
  if Gaussian_FIT:
     MODEL_FIT = "Gaussian"
  elif Lorentzian_FIT:
     MODEL_FIT = "Lorentzian"
  elif Voigtian_FIT: 
      MODEL_FIT = "Voigtian"
  elif ROTATIONAL_FIT:
      MODEL_FIT = "rotational"
  else:
      MODEL_FIT = None
  
  if TYP not in ('min', 'max'):
     print('The type option -t accepts only "min" or "max" values.')
     quit(1)
  
  #Load input data
  if IFN.endswith('.fits') or IFN.endswith('.fit'):
    x_inp, f_inp = load_fits(IFN, EXT_NUMBER, BLA)
  elif IFN.endswith('.npy'):
    idata = np.load(IFN)
    x_inp, f_inp = idata[0], idata[1]
  else:
    x_inp, f_inp = np.loadtxt(IFN, unpack=True, usecols=(0, 1))
  
  #Reverse the input data to detect only maxima
  if TYP == 'min':
     f_inp = 1 - f_inp
     unit = 'Å'
  else:
     unit = 'km/s'
  
  #Remove absolute very high value probably coming from badly formatted input data
  idx, = np.where(abs(f_inp) < 1.e30)
  if list(idx) != []:
     x0 = x_inp[idx]
     f0 = f_inp[idx]
  else:
     x0 = x_inp
     f0 = f_inp
  
  #Threshold on the intensity
  thres0 = THRES0 * np.ptp(f0) + min(f0) 

  #imin and imax are the index of the min and max of the range of selected intensity
  NR = len(x0)
  idx_f0, = np.where(f0[:] > thres0) 
  
  if min(idx_f0):
      i_plus = 1
  else:
      i_plus = 0
  
  imin = min(idx_f0) - i_plus
  
  if max(idx_f0) < NR-1:
      i_plus = 2
  else:
      i_plus = 0   
  
  imax = max(idx_f0) + i_plus
  
  #Selected data
  if imin-N_OFFSET > 0 and imax+N_OFFSET < max(idx):
     x = x0[imin-N_OFFSET:imax+N_OFFSET]
     f = f0[imin-N_OFFSET:imax+N_OFFSET]
  else:
     x = x0[imin:imax]
     f = f0[imin:imax]
  
  #Interpolation of selected data if required
  if N:
     if N < f.size:
        print(f'The number of interpolated points {N} should be larger than the selection {x.size} (change -n option).')
        quit(1)
     xi = np.linspace(min(x), max(x), int(N/x.size)*len(x))
     fi = np.interp(xi, x, f) 
     x, f = xi, fi
  
  #Normalization:
  dx = x[1] - x[0] # use np.diff(x) if x is not uniform
  # case where the abscissa are in decreasing order
  if dx < 0.:
     x = x[::-1]
     f = f[::-1]
     dx = x[1] - x[0]
  dxdx = dx**2
  dxdxdx = dx**3
  
  #Normal derivation without smoothing
  if PURE_DER:
     f0d = np.diff(f0,1)
     f0dd = np.diff(f0d,1)
     f0ddd = np.diff(f0dd,1)
     x0d = x0[:x0.size-1]
     x0dd = x0d[:x0d.size-1]
     x0ddd = x0dd[:x0dd.size-1] 
  
  #Range of selected data
  Dx = max(x)-min(x)
  
  if not SIGMA0:
     if TYP == 'max':
        #SIGMA0 = 2*len(x)/Dx                  # Fully empirical
        #SIGMA0 = np.ptp(x)/4
        SIGMA0 = 3*dx 
     elif TYP == 'min':
        SIGMA0 = 20*len(x)/(Dx*3.e5/x.mean()) # Fully empirical
     #ONE_PASS = True
  
  #Output file names
  parstr = 'doe_'+format(THRES0, '3.2f')+'_'+format(THRES2, '3.2f')+'_'+format(SIGMA0, '3.2f')+'_'+format(N_OFFSET, '0>3.0f')
  OFN1 = op.basename(IFN[:-4])+parstr+'_xp.dat' # Store the positions of the detected peaks
  OFN2 = op.basename(IFN[:-4])+parstr+'_sd.dat' # Store the selected input data and their 3 successives derivatives
  if PLT:
     OFN3 = op.basename(IFN[:-4])+parstr+'.png' 
  
  #Initialization of the counter and the SIGMA value.
  #The loop is performed on SIGMA adding DSIGMA at each iteration.
  count = 0
  SIGMA = SIGMA0
  
  while True:
     # We use here the convolution with a Gaussian kernel (ndimage.Gaussian_filter1d) (f*g)
     # which is equivalent to a smoothed version of df = np.diff(f) / dx
     # Morevover,(f*g)' = f'*g = f*g'
     # 
     #First derivatives:
     gf = ndimage.gaussian_filter1d(f, sigma=SIGMA, order=1, mode='nearest') / dx
     #Second derivatives:
     ggf = ndimage.gaussian_filter1d(f, sigma=SIGMA, order=2, mode='nearest') / dxdx
     #Third derivatives:
     gggf = ndimage.gaussian_filter1d(f, sigma=SIGMA, order=3, mode='nearest') / dxdxdx
  
     
     xp_tot, fp_tot, ggfp_tot = find_peaks(x, f, ggf, gggf)
     
     #Threshold on the second derivative
     thres2 = THRES2 * (max(ggf)-min(ggf))
  
     if len(xp_tot) < 1:
        warn_text = "No peaks detected"
        xp_sel = np.zeros(0)
        xw = np.zeros(0)
        NO_OUTPUTS = True
        break
  
     #Selection of the peak according THRES0 and THRES2
     #crit = np.logical_and(ggfp_tot < -thres2, fp_tot > thres0)
     #Selection of the peak according THRES2 alone
     crit = ggfp_tot < -thres2
     idx = np.where(crit)[0]
     
     xp_sel = xp_tot[idx]
     fp_sel = fp_tot[idx]
     ggfp_sel = ggfp_tot[idx]
  
     if len(xp_sel) < 1:
        warn_text = "No peaks selected"
        xw = np.zeros(0)
        NO_OUTPUTS = True      
        break
  
     #Estimation of the width of each selected peak
     xw, wp = find_widths(x, xp_sel, ggf)
  
     #print(xp_sel, wp, fp_sel)
     #print(xwmin, xwmax)

     xwmin = xw[::2]
     xwmax = xw[1::2]

     if xw.size == 0 or wp.size == 0:
        warn_text = "Warning: DOE cannot measure the width of the second derivative! Maybe the SIGMA parameter is too high? Or use '-n' option to interpolate on more absicssa."
        print(warn_text)
        xw = np.nan*np.ones(xp_sel.size)
        wp = np.nan*np.ones(xp_sel.size)
        xwmin = np.nan*np.ones(xp_sel.size)
        xwmax = np.nan*np.ones(xp_sel.size)        
        #quit(1)

     #print(f'xwmin: {xwmin}')
     #print(f'xwmax: {xwmax}')
  
     # Iteration if the number of selected peak is not equal to the number of width in the 2nd derivative
     if wp.size == 0 or len(set(wp)) == len(xp_sel) or ONE_PASS:
        break
     else:
        count += 1
        SIGMA += DSIGMA
        if False:
           print("count = ", count, ", SIGMA = ", SIGMA)


  # Error of the method (valid for one component only)
  # Lower value for more than one component
  xp_sel_err1, xp_sel_err2 = [], []
  for i in range(xp_sel.size):
    #xp_sel_err.append(method_err(SIGMA, dx, POPT_ERR))
    xp_sel_err1.append((x0[1]-x0[0])/2)
    #Using Zucker 2003, MNRAS, 342, 1291 Eq. page 1293 
    #print(f.size, fp_sel[i], ggfp_sel[i])
    #print(-( f.size * ggfp_sel[i]/fp_sel[i] * 1/(1-1/fp_sel[i]**2) )**(-1))
    xp_sel_err2.append(np.sqrt( abs( -(f.size * ggfp_sel[i]/fp_sel[i] * 1/(1-1/fp_sel[i]**2) )**(-1)) ))

  #Select the type of error
  xp_sel_err = xp_sel_err1  

  if BLA:
     print(f"\nDetection Of Extrema - Version {version} - Thibaut Merle - Verbose mode activated \n")
     if TYP == 'max':
         print(f"Cross-correlation function-like input (TYP = {TYP})")
     else:
        print(f"Spectrum-like input (TYP = {TYP})")
     print(f"Number of the abscissa points: NR = {NR}")
     print(f"Step: dx = {dx:6.4f} {unit}")     
     print("Relative threshold on the CCF/spectrum:    THRES0   = ", THRES0*100, "% => Absolute threshold:", format(thres0, '10.3e'))
     print(f" Number of abscissa points selected:                           {imax-imin:>4d}")
     print(f" Number of additional abscissa points on each side: N_OFFSET = {N_OFFSET:>4d}")
     print(f" Total number of abscissa selected for the derivatives:        {imax-imin+2*N_OFFSET:>4d}")
     if N:
        print(f" Interpolation activated, number of abscissa:              N = {N:>4d}")
     print(f" Interval of abscissa:    [{min(x)}, {max(x)}] {unit}, range of {np.ptp(x):8.3f} {unit}")
     print("Initial derivative smoothing parameter:    SIGMA0   = ", format(SIGMA0, '4.1f'),' ', unit)
     if count:
        print(f"Automatic iteration performed to smooth the derivatives:")
        print(f" Number of iterations:                 IT       = {count}")
        print(f" Step in SIGMA:                        DSIGMA   = {DSIGMA}")
        print(f' Final derivative smoothing parameter: SIGMA    = {SIGMA:5.1f} {unit}')
     print(f"\nNumber of ascending zeros on the third derivative:     {xp_tot.size}")   
     print("\nRelative threshold on the 2nd derivative:  THRES2   = ", THRES2*100, "% => Absolute threshold:", format(-thres2, '10.3e'))
     print(f"\nNumber of RV components kept on the second derivative: {xp_sel.size}")
     for i, j in zip(xp_sel, xp_sel_err1, xp_sel_err2): 
        print('{:8.3f} \u00B1 {:5.3f} {:5.3f}'.format(i, j), end=' ')
     print()


  # fit with model functions of the detected peaks
  if MODEL_FIT:
    if BLA:
        print(f'Fit the CCF with a {MODEL_FIT} mixture model.')
        if NC:
            print(f'/!\ User requires to fit {NC} components.')

  if MODEL_FIT and NC == 0:     
     if MODEL_FIT in ["Gaussian", "Lorentzian"]: 
         npar = 3 # number of parameter for one profile
         guess = zip(xp_sel, wp/2, fp_sel)
         lbounds = zip(xp_sel-COMPTOL, 0*wp, -np.inf*fp_sel)
         ubounds = zip(xp_sel+COMPTOL, +np.inf*wp, +np.inf*fp_sel)
         # [(x0, y0, z0), (x1, y1, z1), ...]
     elif MODEL_FIT == "Voigtian":
         npar = 4
         guess = zip(xp_sel, wp, wp, fp_sel)
         lbounds = zip(xp_sel-COMPTOL, -np.inf*wp, -np.inf*wp, -np.inf*fp_sel)
         ubounds = zip(xp_sel+COMPTOL, +np.inf*wp, +np.inf*wp, +np.inf*fp_sel)
     elif MODEL_FIT == "rotational":
         npar = 4
         guess = zip(xp_sel, 0.6*np.ones(wp.size), wp, fp_sel)
         lbounds = zip(xp_sel-COMPTOL, -np.inf*wp, -np.inf*wp, -np.inf*fp_sel)
         ubounds = zip(xp_sel+COMPTOL, +np.inf*wp, +np.inf*wp, +np.inf*fp_sel)

  if MODEL_FIT and NC > 0:
    vguess = [x[np.argmax(f)] + i*15 for i in range(NC)] #Each component are initially separated by 5 or 15 km/s
    #fguess = [max(f)/(2*i+1) for i in range(NC)] #Each component has an intensity divided by a factor of 3,5,7, etc.   
    fguess = [max(f)/(i+1) for i in range(NC)] #Each component has an intensity divided by a factor of 2,3,4, etc.   
    if MODEL_FIT in ["Gaussian", "Lorentzian"]:
        npar = 3 # number of parameter for one profile
        guess = zip(vguess, NC*[5.], fguess)
        lbounds = zip(vguess-COMPTOL*np.ones(NC), 0*np.ones(NC), 0*np.ones(NC))
        ubounds = zip(vguess+COMPTOL*np.ones(NC), +np.inf*np.ones(NC), +np.inf*np.ones(NC))
    elif MODEL_FIT in ["Voigtian", "rotational"]:
        npar = 4
        vguess = [x[np.argmax(f)] + i*3 for i in range(NC)]
        guess = zip(vguess, NC*[0.6], NC*[5.], fguess)
        lbounds = zip(vguess-COMPTOL*np.ones(NC), -np.inf*np.ones(NC), 0*np.ones(NC), -np.inf*np.ones(NC))
        ubounds = zip(vguess+COMPTOL*np.ones(NC), +np.inf*np.ones(NC), +np.inf*np.ones(NC), +np.inf*np.ones(NC))                    
  
  if not MODEL_FIT and NC > 0:
    print(f"'-m' option runs with the fitting options [-G, -L, -V, -R]")
    quit(1)

  #Finalization of the guess and boundaries
  if MODEL_FIT:
     # [c, x0, y0, z0, x1, y1, z1, ...]
     guess = [min(f0)] + [i for sublist in guess for i in sublist]
     lbounds = [min(f0) - 0.1*np.ptp(f0)] + [i for sublist in lbounds for i in sublist]
     ubounds = [min(f0) + 0.3*np.ptp(f0)] + [i for sublist in ubounds for i in sublist]
     
     if MODEL_FIT == "Gaussian":
        mmf = gmm   # mmf = mixture model function
     elif MODEL_FIT == "Lorentzian":
        mmf = lmm
     elif MODEL_FIT == "Voigtian":
        mmf = vmm
     elif MODEL_FIT == "rotational":
        mmf = rmm
    
     try:

        popt, pcov = curve_fit(mmf, x, f, p0=guess, bounds=(lbounds, ubounds))
        perr = np.sqrt(np.diag(pcov))
     except:
        print("MODEL_FIT FAILED")
        MODEL_FIT = None
  
     #Another way to fit data
     if LEASTSQ:    
        errfunc = lambda guess, x, f: gmm(x, *guess) - f
        try:
            popt2, pcov2, infodict, mesg, success = leastsq(errfunc, guess.copy(), args=(x, f), full_output=1, ftol = 1.e-3, xtol=1.e-3)
            perr2 = np.sqrt(np.diag(pcov2))
        except ValueError:
            popt2 = np.nan*np.array(guess)
            perr2 = np.nan*np.array(guess)
        #print(p1, cov_x, infodict, mesg, success)


  if MODEL_FIT:
     fit = mmf(x, *popt)
     c_mmf = popt[0] # offset of the Gaussian/Lorentzian/Voigtian mixture
     mu_mmf = popt[1::npar]   # radial velocities
     mu_mmf_err = perr[1::npar] # error on radial velocities
     n_mmf = popt[npar::npar]   # normalisation factor
     idx = mu_mmf.argsort() # Index of the sorted increasing radial velocities     
     mu_mmf  = mu_mmf[idx]
     n_mmf = n_mmf[idx] 

     #For the Gaussian fit with LEASTSQ
     if LEASTSQ:
        fit2 = gmm(x, *popt2) 
        mu_gmm2 = popt2[1::npar]   # radial velocities     
        mu_gmm_err2 = perr2[1::npar] # error on radial velocities
  
     if MODEL_FIT in ["Gaussian", "Lorentzian"]: 
        sig_mmf = popt[2::npar] # standard deviation of the Gaussian profile / HWHM of the Lorentzian profile 
        sig_mmf = sig_mmf[idx]
     elif MODEL_FIT == "Voigtian":
        alp_mmf = popt[2::4] # HWHM of the Gaussian profile
        gam_mmf = popt[3::4] # HWHM of the Lorentzian profile
        alp_mmf = alp_mmf[idx]
        gam_mmf = gam_mmf[idx]
     elif MODEL_FIT == 'rotational':
        eps_mmf = popt[2::4] # Limb-darkening coefficient (set to 0.6 for solar-like stars)
        vr_mmf  = popt[3::4] # Rotational velocity
        eps_mmf = eps_mmf[idx]
        vr_mmf  = vr_mmf[idx]

  
  #Write output files
  
  if xp_sel.size and wp.size and xwmin.size and xwmax.size:
     # Only if 2/3 peaks has a common width 
     if len(xw)/2 != len(xp_sel):
        wp = np.array(len(xp_sel)*list(wp))
        xwmin =  np.array(len(xp_sel)*list(xwmin))
        xwmax =  np.array(len(xp_sel)*list(xwmax))

  #Case where we force the number of components to fit and that number is larger to what DOE get from the derivatives
  if (xp_sel.size == 1) & (NC > 1):
    xp_sel = np.array(NC*list(xp_sel))
    xp_sel_err = np.array(NC*list(xp_sel_err))
    fp_sel = np.array(NC*list(fp_sel))
    wp = np.array(NC*list(wp))
    xwmin =  np.array(NC*list(xwmin))
    xwmax =  np.array(NC*list(xwmax))
  elif (xp_sel.size > 1) & (NC > xp_sel.size >1):
    print('Warning: output files will not reflect the number of components reflected for the fit.')

  if TYP == 'min':
    fp_sel = 1 - fp_sel

  if NO_OUTPUTS:
      pass
  else:
     if MODEL_FIT: # Add the optimal parameters determined by the Gaussian fit process at the end of OFILE1
        if MODEL_FIT == "Gaussian":
           np.savetxt(OFN1, list(zip(xp_sel, xp_sel_err, fp_sel, wp, xwmin, xwmax, SIGMA*np.ones(wp.size), count*np.ones(wp.size), mu_mmf, mu_mmf_err, sig_mmf, n_mmf, c_mmf*np.ones(wp.size))), fmt='%8.3f %8.3f %10.5f %7.3f %8.3f %8.3f %8.3f %3.0f %8.3f %5.3f %7.3f %6.3f %6.4f', header='x x_err f width xwmin xwmax sigma it, Gaussian fit (mu, mu_err, sigma, n) + offset')
        elif MODEL_FIT == "Lorentzian":
           np.savetxt(OFN1, list(zip(xp_sel, xp_sel_err, fp_sel, wp, xwmin, xwmax, SIGMA*np.ones(wp.size), count*np.ones(wp.size), mu_mmf, mu_mmf_err, sig_mmf, n_mmf, c_mmf*np.ones(wp.size))), fmt='%8.3f %8.3f %10.5f %7.3f %8.3f %8.3f %8.3f %3.0f %8.3f %5.3f %7.3f %6.3f %6.4f', header='x x_err f width xwmin xwmax sigma it, Lorentzian fit (mu, mu_err, gamma, n) + offset')
        elif MODEL_FIT == "Voigtian":
           np.savetxt(OFN1, list(zip(xp_sel, xp_sel_err, fp_sel, wp, xwmin, xwmax, SIGMA*np.ones(wp.size), count*np.ones(wp.size), mu_mmf, mu_mmf_err, alp_mmf, gam_mmf, n_mmf, c_mmf*np.ones(wp.size))), fmt='%8.3f %10.5f %7.3f %7.3f %8.3f %8.3f %8.3f %3.0f %8.3f %5.3f %7.3f %7.3f %6.3f %6.4f', header='x x_err f width xwmin xwmax sigma it, Voigtian fit (mu, mu_err, alpha, gamma, n) + offset')
        elif MODEL_FIT == 'rotational':
           np.savetxt(OFN1, list(zip(xp_sel, xp_sel_err, fp_sel, wp, xwmin, xwmax, SIGMA*np.ones(wp.size), count*np.ones(wp.size), mu_mmf, mu_mmf_err, eps_mmf, vr_mmf, n_mmf, c_mmf*np.ones(wp.size))), fmt='%8.3f %8.3f %10.5f %7.3f %8.3f %8.3f %8.3f %3.0f %8.3f %5.3f %7.3f %7.3f %6.3f %6.4f', header='x x_err f width xwmin xwmax sigma it, rotational fit (mu, mu_err, epsilon, vsini, n) + offset')
     else:
        #print(xp_sel, xp_sel_err, fp_sel, wp, xwmin, xwmax, SIGMA*np.ones(wp.size), count*np.ones(wp.size))
        np.savetxt(OFN1, list(zip(xp_sel, xp_sel_err, fp_sel, wp, xwmin, xwmax, SIGMA*np.ones(wp.size), count*np.ones(wp.size))), fmt='%8.3f %8.3f %9.3f %7.3f %8.3f %8.3f %8.3f %3.0f', header='x x_err f width xwmin xwmax sigma it')
  
     if NO_OUTPUT2:
        pass
     else:
        np.savetxt(OFN2, list(zip(x, f, gf, ggf, gggf)), fmt='%10.3f %10.3e %10.3e %10.3e %10.3e', header='x f gf ggf gggf')
  
  if BLA:
    if MODEL_FIT:
       print(f'Tolerance for each component: COMPTOL = {COMPTOL} {unit}')
       #print('guess: ', guess)
       #print('popt: ', popt)
       #print('pcov:', pcov)
       for i, j in zip(mu_mmf, mu_mmf_err):
           print('{:8.3f} \u00B1 {:5.3f}'.format(i, j), end=' ')    
       print()
    if NO_OUTPUTS:
       print('No output required.')
    elif NO_OUTPUT2:
       print("Outputs:")
       print(' Successive derivatives output:', OFN2) 
       print(' Components output:            ', OFN1)
    if PLT:
       print(' Control plot:                 ', OFN3, end='\n')
  
  if NC:
    npeaks = NC  
  else:
    npeaks = xp_sel.size

  #Standard output
  if not BLA:
      #if MODEL_FIT:
      #    if NC:
      #      npeaks = NC 
      #    else:
      #      npeaks = xp_sel.size 
      #else:
      #    npeaks = xp_sel.size 
      
      print(IFN, npeaks, end=' ')
      
      if xp_sel.size == 0:
         print(warn_text+' (IT = '+str(count)+')')
      if not MODEL_FIT:
         for i, j in zip(xp_sel, xp_sel_err): 
            print('{:8.3f} \u00B1 {:5.3f}'.format(i, j), end=' ')
      else:
         for i, j in zip(mu_mmf, mu_mmf_err):
            print('{:8.3f} \u00B1 {:5.3f}'.format(i, j), end=' ')
      print()
  
  #Plotting:
  
  if PLT:
  
     title = op.basename(IFN) + ' (' + str(len(x_inp)) + ') ' + '\n'
     title += '[THRES0 = '+str(int(THRES0*100))+'%, THRES2 = '+str(int(THRES2*100))+'%, SIGMA = '+format(SIGMA, '4.1f')+' '+unit

     if N:
        title += ', N = '+str(N)
     if not ONE_PASS:
        title += ', IT = '+str(count)

     if NC:
        title += ', NC = '+str(NC)
  
     title = title + ']'
  
     plt.figure(figsize=(8, 7), dpi=100, facecolor='w', edgecolor='k')
  
     plt.subplots_adjust(hspace=0.05)
     plt.rc('font', size=13)
     
     plt1 = plt.subplot(411)
     plt2 = plt.subplot(412, sharex=plt1)
     plt3 = plt.subplot(413, sharex=plt1)
     plt4 = plt.subplot(414, sharex=plt1)
  
     for i in (plt1, plt2, plt3, plt4):
        for j in xp_tot:
           i.axvline(x=j, ls='solid', color='0.95', lw=0.2)
  
     if TYP == 'max':
        plt1.plot(x0, f0, '-', lw=1, color='0.50', label='original')
        plt1.plot(x, f, 'ko', ms=1.5, label='selection')
        if MODEL_FIT:
           #plt1.plot(x, fit, 'b-', label=MODEL_FIT + ' fit')
           plt1.plot(x0, mmf(x0, *popt), 'b-', label=MODEL_FIT + ' fit')
           if npeaks > 1:
              #Plot individual components
              if MODEL_FIT == 'Gausian':
                 for mu, sigma, n, offset in zip(mu_mmf, sig_mmf, n_mmf, c_mmf*np.ones(wp.size)):
                    ytemp = n*np.exp(-(x0-mu)**2/(2*sigma**2)) + offset
                    plt1.plot(x0, ytemp, '-b', lw=1, alpha=0.3)
              if MODEL_FIT == 'Lorentzian':
                 for mu, gam, n, offset  in zip(mu_mmf, sig_mmf, n_mmf, c_mmf*np.ones(wp.size)):             
                    ytemp = n * gam / np.pi / ((x0 - mu)**2 + gam**2) + offset
                    plt1.plot(x0, ytemp, '-b', lw=1, alpha=0.3)
              if MODEL_FIT == 'Voigtian':
                 for mu, sig, gam, n, offset in zip(mu_mmf, alp_mmf, gam_mmf, n_mmf, c_mmf*np.ones(wp.size)): 
                    ytemp = n * np.real(special.wofz(((x0 - mu) + 1j * gam) / sig / np.sqrt(2))) + offset       
                    plt1.plot(x0, ytemp, '-b', lw=1, alpha=0.3)
              if MODEL_FIT == 'rotational':
                 for mu, eps, vr, n, offset in zip(mu_mmf, eps_mmf, vr_mmf, n_mmf, c_mmf*np.ones(wp.size)):
                    xx = ((x0-mu)/vr)**2
                    r0 = (2+(np.pi/2-2)*eps) / (np.pi*vr*(1-eps/3))
                    xx[abs(xx)>1] = 1
                    ytemp = n/r0 * (2*(1-eps)*(1-xx)**0.5 + 0.5*np.pi*eps*(1-xx)) / (np.pi*vr*(1-eps/3.))
                    plt1.plot(x0, ytemp, '-b', lw=1, alpha=0.3)
           if LEASTSQ:    
             plt1.plot(x, fit2, 'g-', label='Gaussian fit (leastsq)')          
        #plt1.axhline(y=min(f0), ls='dotted', color='0.00')
        plt1.axhline(y=thres0, ls='solid', color='r', lw=0.5)
     elif TYP == 'min':
        plt1.plot(x0, 1-f0, '-', lw=1, color='0.50', label='original')
        plt1.plot(x, 1-f, 'ko', ms=1.5, label='selection')
        if MODEL_FIT:
           plt1.plot(x, 1-fit, 'b-', label=MODEL_FIT + ' fit')      
        plt1.axhline(y=1, ls='dotted', color='0.00', alpha=0.4, lw=0.5)
        plt1.axhline(y=1-thres0,ls='solid', color='r', lw=0.5)
  
    
     if LAB:
        plt1.set_title(LAB,fontsize=15)
     elif not NOT:
        plt1.set_title(title, fontsize=8) 
  
  
     plt2.plot(x, gf, 'k-o', ms=1.5, label='1st smoothed')
     
     plt3.plot(x, ggf, 'k-o', ms=1.5, label='2nd smoothed')
     plt3.axhline(y=-thres2, ls='solid', color='r', lw=0.5)
     for i in range(0, len(xw)-1, 2):
        plt3.plot(xw[i:i+2], np.array([0,0]), 'k-', lw=2)
  
     plt4.plot(x, gggf, 'k-o', ms=1.5, label='3rd smoothed')
  
     if TYP == 'max':
        plt4.set_xlabel('$v$ [km/s]', size=15)
     elif TYP == 'min':
        plt4.set_xlabel('$\\lambda$ [Å]', size=15)

     if PURE_DER:
        plt2.plot(x0d, f0d, '-', color='0.8', lw=0.5, label='1st derivative', zorder=0)
        plt3.plot(x0dd, f0dd, '-', color='0.8', lw=0.5, label='2nd derivative', zorder=0)
        plt4.plot(x0ddd, f0ddd, '-', color='0.8', lw=0.5, label='3rd derivative', zorder=0)
  
     for i in (plt1, plt2, plt3, plt4):
        i.legend(frameon=True, loc=1, fontsize=10)
        i.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1,2,5,9,10], prune='lower'))
        i.ticklabel_format(axis='y', style='plain', scilimits=(-2,2), useOffset=False)
        i.axvline(x=0, ls='dotted', color='0.00', alpha=0.4, lw=0.5)
        i.axhline(y=0, ls='dotted', color='0.00', alpha=0.4, lw=0.5)

        for j in xp_sel:
           i.axvline(x=j, ls='solid', color='0.00', lw=0.7)
  
        if MODEL_FIT:
          for j in mu_mmf:
             plt1.axvline(x=j, ls='solid', color='b', lw=0.7)
             if LEASTSQ:
                plt1.axvline(x=j, ls='solid', color='g', lw=0.7)
        
     for i in (plt1, plt2, plt3):
        i.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, labelbottom=False)
     plt4.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, labelbottom=True)

  
     if XMIN and XMAX:
        plt1.set_xlim([XMIN, XMAX])
  
     if XMIN and not XMAX:
        plt1.set_xlim([XMIN, -XMIN])
  
     if XMAX and not XMIN:
        plt1.set_xlim([-XMAX, XMAX])
  
     if not XMIN and not XMAX:
        if TYP == 'max':
           plt1.set_xlim([min(x)-1.2*Dx, max(x)+1.2*Dx])  
           #plt1.set_xlim([min(x), max(x)])
        else:
           plt1.set_xlim([min(x), max(x)])
  
     if YMIN and not YMAX:
        plt1.set_ylim([YMIN, max(f)])
     
     if YMAX and not YMIN:
        plt1.set_ylim([min(f), YMAX])
  
     if YMIN and YMAX:
        plt1.set_ylim([YMIN, YMAX])
  
     #if PURE_DER:
     #   plt2.set_ylim([-5e-3,5e-3])
     #   plt3.set_ylim([-3e-4,3e-4])
     #   plt4.set_ylim([-5e-5,5e-5])
  
  
     if xp_sel.size > NP_MAX:
        plt4.text(0.05, 0.1, str(xp_sel.size)+' components', transform=plt1.transAxes) #, fontsize=11)
     else:
  
        if TYP == 'max':
           base = 'v'
        else:
           base = 'l'      
  
        #Print velocities from third derivatives and widths
        for ind, (itm1, itm2) in enumerate(zip(xp_sel[:NP_MAX], xp_sel_err[:NP_MAX])):
           textv = '$v_{}$ = {:7.3f} \u00B1 {:.3f} km/s'.format(ind+1, itm1, itm2)
           try:
              textw = '$w_{}$ = {:7.3f} km/s'.format(ind+1, wp[ind]) 
           except IndexError:
              pass           
           ycord = 0.8-ind/6.
           plt4.text(0.02, ycord, textv, transform=plt4.transAxes, fontsize=9)
           try:
              plt3.text(0.02, ycord, textw, transform=plt3.transAxes, fontsize=9)
           except NameError:
              pass

        #Print velocities from fit
        for ind in range(max(NC, xp_sel.size)):

           if MODEL_FIT:
              if MODEL_FIT in ['Gaussian', 'Lorentzian']:
                width = 2*sig_mmf[ind]
              elif MODEL_FIT == 'Voigtian':
                width = 2*alp_mmf[ind] / np.sqrt(2 * np.log(2))
              elif MODEL_FIT == 'rotational':
                width = vr_mmf[ind] 
              try:
                 textv_gf = '$v_{}$={:7.3f}\u00B1{:.3f} ($w_{}$={:.3f}) km/s'.format(ind+1, mu_mmf[ind], mu_mmf_err[ind], ind+1, width)
              except IndexError:
                 pass
              if LEASTSQ:
                try:
                    textv_gf2 = '$v_{}$={:7.3f}\u00B1{:.3f} km/s'.format(ind+1, mu_gmm2[ind], mu_gmm_err2[ind])
                except IndexError:
                    pass                

              ycord = 0.8-ind/6.
              plt1.text(0.02, ycord, textv_gf, transform=plt1.transAxes, color='blue', fontsize=9)
              if LEASTSQ:
                ycord2 = 0.4-ind/6. 
                plt1.text(0.02, ycord2, textv_gf2, transform=plt1.transAxes, color='green', fontsize=9)

  
     #plt3.set_xticklabels([])
     #plt4.set_xticklabels(['', -75, -50, -25, 0, 25, 50, 75])

     plt.savefig(OFN3, dpi=150, orientation='landscape', bbox_inches='tight')    
  
     if PLT == 2:
        plt.show()
     else:
        plt.draw()


#######################################################################################


if __name__ == '__main__':

  #Default parameters
  version = 2.0 # DOE version
  TYP = 'max'   # Kind of detection (min for absorption spectrum, max for cross-correlation function)
  THRES0 = 0.3  # Default threshold on the CCF (in % of the full amplitude)
  THRES2 = 0.1  # Default threshold on the 2nd derivative (in % of the full amplitude)
  PLT = False
  MODEL_FIT = None # Name of the mixture model function use to fit the peaks
  NP_MAX = 5 # Maximum number of peak positions written on the plot
  # 
  DSIGMA = 2    # [km/s or Å] Step on the SIGMA parameter use for iteration  (when ONE_PASS = False)
  N_OFFSET = 10 # Number of points selected on the left and right of the selected CCF and derivatives in addition to imin and imax
  #POPT_ERR = [0.0000, -0.0000, -0.0000,  0.0000] # Internal + CCF sampling error estimation using script_fit_bias_precision.py
  
  LEASTSQ = False #Internal option to compare with a Gaussian fit from the scipy.optimize.leastsq

  np.set_printoptions(linewidth=sys.maxsize)

  #Command line arguments
  parser = ap.ArgumentParser(description=f'Detection Of Extrema (DOE). Version {version}. Use successive derivatives to find blended components.', epilog='2022-11-06 Thibault Merle')
  #Positional
  parser.add_argument('ifn', help='Input filename [x, f] in unicode format or FITS (with CRVAL1 and CDELT1 defined).')
  #Optional -- general behaviour
  group1 = parser.add_argument_group('General behaviour')
  group1.add_argument('-v', '--verbose', action='store_true', default=False, help='verbosity')
  group1.add_argument('-t', '--type', default=TYP, help='Extrema to find: min for absorption spectrum, max for cross-correlation function (default = '+TYP+').')
  group1.add_argument('-p', '--plot', action='count', default=None, help='"-p" write a plot, "-pp" and display it.')
  group1.add_argument('-e', '--extension_number', type=int, default=0, help='Input FITS extension number (default: 0)')
  group1.add_argument('-i', '--ignore_outputs', action='store_true', default=False, help='Do not write the outputs *_xp.dat and *_sd.dat') 
  group1.add_argument('-i2', '--ignore_output2', action='store_true', default=False, help='Do not write the output *_sd.dat') 
  #Optional -- Code parameters
  group2 = parser.add_argument_group('Code parameters')
  group2.add_argument('-c', '--thres0', type=float, default=THRES0, help='threshold on the intensity in percent of the full scale (default = '+str(THRES0)+').')
  group2.add_argument('-d', '--thres2', type=float, default=THRES2, help='threshold on the 2nd derivative of the intensity in percent of the full scale (default = '+str(THRES2)+').')
  group2.add_argument('-s', '--sigma', type=float, default=None, help='Standard deviation for the Gaussian kernel in km/s or Å used to smooth the derivatives.')
  group2.add_argument('-n', type=int, default=None, help='Number of interpolated points required for the input data (should be larger than the number of selected x points.)')
  group2.add_argument('-no', '--n_offset', type=int, default=N_OFFSET, help='Number of points below/above the CCF threshold taken on the left and right of the selected CCF')
  group2.add_argument('-1', '--one_pass', action='store_true', default=False, help='One pass without automatic adjustment of SIGMA (no iteration)')
  #Optional -- Fitting parameters
  group4 = parser.add_argument_group('Fitting parameters')
  group4.add_argument('--comptol', type=float, default=20, help='Fit tolerance on x position for each x component in km/s or Å (default: 100).')
  group4.add_argument('-G', '--Gaussian_fit', action='store_true', default=False, help='Fit the detected peaks with Gaussian function and give the velocities fitted on the line command.')
  group4.add_argument('-L', '--Lorentzian_fit', action='store_true', default=False, help='Fit the detected peaks with Lorentzian function and give the velocities fitted on the line command. /!\ only tested with TYP = "max"')
  group4.add_argument('-V', '--Voigtian_fit', action='store_true', default=False, help='Fit the detected peaks with Voigt function and give the velocities fitted on the line command. /!\ only tested with TYP = "max"')
  group4.add_argument('-R', '--rotational_fit', action='store_true', default=False, help='Fit the detected peaks with rotational function and give the velocities fitted on the line command. /!\ only tested with TYP = "max"')
  group4.add_argument('-m', type=int, default=0, help='Number of components to fit (it forces DOE to fit this number of RV components).')
  #Optional -- Graphical parameters
  group3 = parser.add_argument_group('Graphical parameters')
  group3.add_argument('-xmin', '--xmin', type=float, default=None, help='x minimum value for the plot')
  group3.add_argument('-xmax', '--xmax', type=float, default=None, help='x maximum value for the plot')
  group3.add_argument('--ymin', type=float, default=None, help='y minimum value for the CCF/flux of the plot')
  group3.add_argument('--ymax', type=float, default=None, help='y maximum value for the CCF/flux of the plot')
  group3.add_argument('-not', '--no_title', action='store_true', default=False, help='Do not display the title over the plot.')
  group3.add_argument('-lab', '--label', type=str, default=None, help='Display a label in place of title.')
  group3.add_argument('--pure_derivative', action='store_true', default=False, help='Show the pure derivative without smoothing.')
  
  args = parser.parse_args()
  
  #Read input line command arguments 
  IFN = args.ifn
  BLA = args.verbose
  PLT = args.plot
  TYP = args.type
  ONE_PASS = args.one_pass # If True no iteration on SIGMA allowed
  THRES0 = args.thres0
  THRES2 = args.thres2
  SIGMA0 = args.sigma
  N = args.n
  N_OFFSET = args.n_offset
  XMIN = args.xmin
  XMAX = args.xmax
  YMIN = args.ymin
  YMAX = args.ymax
  NOT = args.no_title
  LAB = args.label
  PURE_DER = args.pure_derivative
  EXT_NUMBER = args.extension_number
  Gaussian_FIT = args.Gaussian_fit
  Lorentzian_FIT = args.Lorentzian_fit
  Voigtian_FIT = args.Voigtian_fit
  ROTATIONAL_FIT = args.rotational_fit
  COMPTOL = args.comptol
  NO_OUTPUTS = args.ignore_outputs
  NO_OUTPUT2 = args.ignore_output2
  NC = args.m
  

  run_doe(IFN, BLA, PLT, TYP, ONE_PASS, THRES0, THRES2, SIGMA0, N, N_OFFSET, XMIN, XMAX, YMIN, YMAX, NOT, LAB, PURE_DER, EXT_NUMBER,\
        Gaussian_FIT, Lorentzian_FIT, Voigtian_FIT, ROTATIONAL_FIT, COMPTOL, NO_OUTPUTS, NO_OUTPUT2, NC)

