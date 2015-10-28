from __future__ import division, print_function
import numpy as np
import os
import scipy.interpolate as spi
import scipy.signal as sps
import sys
import matplotlib.mlab as mplm
from t_getconsts import t_getconsts
from t_vuf import t_vuf
import matplotlib as mpl
import matplotlib.dates as dates

np.set_printoptions(precision=8, suppress=True)


def fourpad(conin):
    conin = np.array(conin, dtype='|S4')
    for i, con in enumerate(conin):
        conin[i] = con.ljust(4)

    return conin


def constituents(minres, constit, shallow, infname, infref, centraltime):
    """[name,freq,kmpr]=constituents(minres,infname) loads tidal constituent
     table (containing 146 constituents), then picks out only the '
     resolvable' frequencies (i.e. those that are MINRES apart), base on
     the comparisons in the third column of constituents.dat. Only
     frequencies in the 'standard' set of 69 frequencies are actually used.
     Also return the indices of constituents to be inferred.
     If we have the mat-file, read it in, otherwise create it and read
     it in!
     R Pawlowicz 9/1/01
     Version 1.0

        19/1/02 - typo fixed (thanks to  Zhigang Xu)
     Compute frequencies from astronomical considerations.
    """
    if minres > 1 / (18.6*365.25*24):
        # Choose only resolveable pairs for short
        const, sat, cshallow = t_getconsts(centraltime)
        # Time series
        ju = np.flatnonzero(const['df'] >= minres)
    else:
        # Choose them all if > 18.6 years.
        const, sat, cshallow = t_get18consts(centraltime)

        ju = np.array([range(2,
                             (max(const['freq'].shape)+1))]).reshape(1, -1).T
        # Skip Z0
        for ff in range(1, 3):
            # loop twice to make sure of neightbouring pairs
            jck = np.flatnonzero(diff(const['freq'][ju]) < minres)
            if (max(jck.shape) > 0):
                jrm = jck
                jrm = jrm + (abs(const['doodsonamp'][ju[(jck+1-1)]]) <
                             abs(const['doodsonamp'][ju[(jck-1)]]))
# disp('  Warning! Following constituent pairs violate Rayleigh criterion')
#               for ick in range(1, (max(jck.shape) +1)):
#                    disp('     ' + const.name(ju[(jck[(ick -1)] -1)], :)
# + ' vs ' + const.name(ju[(jck[(ick -1)] + 1 -1)], :) + ' - not using ' +
# const.name(ju[(jrm[(ick -1)] -1)], :))
                ju[(jrm-1)] = np.array([])

    if constit.size != 0:
        # Selected if constituents are specified in input.
        ju = np.array([], dtype=int)
        for k in range(0, (constit.shape[0])):
            j1 = np.where(const['name'] == constit[k])[0]
            if (j1.size == 0):
                print("Can't recognize name " +
                      constit[k].decode(sys.stdout.encoding) +
                      ' for forced search')
            else:
                ju = np.append(ju, j1)

        # sort in ascending order of frequency.
        II = np.argsort(const['freq'][ju])
        ju = ju[II]

    # cout
    # disp(['   number of standard constituents used: ',int2str(length(ju))])
    if shallow.size != 0:
        # Add explictly selected shallow water constituents.
        for k in range(1, (shallow.shape[0]+1)):
            j1 = strmatch(shallow[(k-1), :], const['name'])
            if (0 in j1.shape):
                disp("Can't recognize name " + shallow[(k-1), :] +
                     ' for forced search')
            else:
                if np.isnan(const['ishallow'][j1]):
                    disp(shallow[(k-1), :] +
                         ' Not a shallow-water constituent')
                disp('   Forced fit to ' + shallow[(k-1), :])
                ju = np.array([ju, j1]).reshape(1, -1)
    nameu = const['name'][ju]
    fu = const['freq'][ju]

    # Check if neighboring chosen constituents violate Rayleigh criteria.
    jck = np.flatnonzero(np.diff(fu) < minres)
    # cout
    # if (length(jck)>0)
    # disp(['  Warning! Following constituent
    # pairs violate Rayleigh criterion']);
    # for ick=1:length(jck);
    # disp(['     ',nameu(jck(ick),:),'  ',nameu(jck(ick)+1,:)]);
    # end;
    # end
    # For inference, add in list of components to be inferred.
    fi = np.array([])
    namei = np.array([])
    jinf = np.array([])
    jref = np.array([])
    if infname.size != 0:
        fi = np.zeros(shape=(infname.shape[0], 1), dtype='float64')
        namei = np.zeros(shape=(infname.shape[0], 4), dtype='float64')
        jinf = np.zeros(shape=(infname.shape[0], 1), dtype='float64') + NaN
        jref = np.zeros(shape=(infname.shape[0], 1), dtype='float64') + NaN
        for k in range(1, (infname.shape[0]+1)):
            j1 = strmatch(infname[(k-1), :], const.name)
            if (0 in j1.shape):
                disp("Can't recognize name" +
                     infname[(k-1), :] + ' for inference')
            else:
                jinf[(k-1)] = j1
                fi[(k-1)] = const['freq'][j1]
                namei[(k-1), :] = const['name'][j1, :]
                j1 = strmatch(infref[(k-1), :], nameu)
                if (0 in j1.shape):
                    disp("Can't recognize name " + infref[(k-1), :] +
                         ' for as a reference for inference')
                else:
                    jref[(k-1)] = j1
                    print('   Inference of ' + namei[(k-1), :] +
                          ' using ' + nameu[(j1-1), :] + '\\n')
        jinf[(isnan(jref)-1)] = NaN
    return nameu, fu, ju, namei, fi, jinf, jref


def fixgaps(x):
    """FIXGAPS: Linearly interpolates gaps in a time series
     YOUT=FIXGAPS(YIN) linearly interpolates over NaN in the input time
     series (may be complex), but ignores trailing and leading NaNs.
     R. Pawlowicz 11/6/99
     Version 1.0
    """
    y = x
    bd = np.isnan(x)
    gd = np.flatnonzero(~bd)
    idx = np.array([range(1, ((np.min(gd)-1)+1)),
                    range((np.max(gd) + 1), (len(gd)))]).reshape(1, -1)-1
    if idx.size != 0:
        bd[idx] = 0
        y[(bd-1)] = interp1(gd, x[(gd-1)], np.flatnonzero(bd))

    return y


def cluster(ain, clusang):
    """CLUSTER: Clusters angles in rows around the angles in the first
     column. CLUSANG is the allowable ambiguity (usually 360 degrees but
     sometimes 180).
    """

    makearray = (ain - np.repeat(ain[:, 0],
                 ain.shape[1]).reshape(ain.shape[0], ain.shape[1]))
    ii = makearray > clusang / 2
    ain[(ii)] = ain[(ii)] - clusang
    ii = (ain - np.repeat(ain[:, 0],
          ain.shape[1]).reshape(ain.shape[0], ain.shape[1])) < - clusang / 2
    ain[(ii)] = ain[(ii)] + clusang
    return ain


def noise_realizations(xres, fu, dt, nreal, errcalc):
    """NOISE_REALIZATIONS: Generates matrices of noise (with correct
     cross-correlation structure) for bootstrap analysis.

     R. Pawlowicz 11/10/00
     Version 1.0
    """
    if errcalc == 'cboot':
        fband, Pxrave, Pxiave, Pxcave = residual_spectrum(xres, fu, dt)
        aaa, bbb = Pxcave.shape
        Pxcave = np.zeros((aaa, bbb), dtype='float64')
        # For comparison with other technique!
        # print('**** Assuming no covariance between u and v errors!******\n');
    else:
        if errcalc == 'wboot':
            fband = np.array([0, 0.5]).reshape(1, -1)
            nx = max(xres.shape)
            A = np.cov(np.real(xres), np.imag(xres)) / nx
            Pxrave = np.array([A[0, 0]])
            Pxiave = np.array([A[1, 1]])
            Pxcave = np.array([A[0, 1]])
        else:
            sys.exit("Unrecognized type of bootstap analysis specified: '" +
                     errcalc + "'")

    nfband = fband.shape[0]
    Mat = np.zeros(shape=(4, 4, nfband), dtype='float64')
    for k in range(1, (nfband+1)):
        # The B matrix represents the covariance matrix for the vector
        # [Re{ap} Im{ap} Re{am} Im{am}]' where Re{} and Im{} are real and
        # imaginary parts, and ap/m represent the complex constituent
        # amplitudes for positive and negative frequencies when the input
        # is bivariate white noise. For a flat residual spectrum this works
        # fine.
        # This is adapted here for "locally white" conditions, but I'm still
        # not sure how to handle a complex sxy, so this is set to zero
        # right now.
        p = (Pxrave[(k-1)] + Pxiave[(k-1)]) / 2
        d = (Pxrave[(k-1)] - Pxiave[(k-1)]) / 2
        sxy = Pxcave[(k-1)]
        B = np.hstack([p, 0, d, sxy,
                       0, p, sxy, -d,
                       d, sxy, p, 0,
                       sxy, -d, 0, p]).reshape(4, 4)
        # Compute the transformation matrix that takes uncorrelated white
        # noise and makes noise with the same statistical structure as the
        # Fourier transformed noise.

        D, V = np.linalg.eigh(B)

        # next five lines are horrible coding/math someone should check it over
        # swap so the vectors match matlab, should check if this always holds
        V[[1, 0], [1, 3]] = V[[0, 1], [3, 1]]
        V[[3, 2], [1, 3]] = V[[2, 3], [3, 1]]
        # total cludge to deal with bad zeroing in eigh
        D[((D < 0) & (D > -0.00000000001))] = 0

        Mat[:, :, (k-1)] = np.dot(V, np.diag(np.sqrt(D)))
        # print Mat
    # Generate realizations for the different analyzed constituents.
    N = np.zeros(shape=(4, nreal), dtype='float64')
    NM = np.zeros(shape=(max(fu.shape), nreal), dtype='complex128')
    NP = np.zeros(shape=(max(fu.shape), nreal), dtype='complex128')

    for k in range(0, fu.shape[0]):
        l = np.squeeze(np.flatnonzero(np.all([fu[k] > fband[:, 0] , fu[k] < fband[:, 1]],axis=0)))
        N = np.hstack([np.zeros(4,).reshape(-1,1), np.dot(np.squeeze(Mat[:, :, l]), np.random.randn(4, nreal-1))])
        NP[(k), :] = (N[0, :]+1j*N[1, :])
        NM[(k), :] = (N[2, :]+1j*N[3, :])

    return NP, NM

def residual_spectrum(xres, fu, dt):
    """RESIDUAL_SPECTRUM: Computes statistics from an input spectrum over
     a number of bands, returning the band limits and the estimates for
     power spectra for real and imaginary parts and the cross-spectrum.

     Mean values of the noise spectrum are computed for the following
     8 frequency bands defined by their center frequency and band width:
     M0 +.1 cpd; M1 +-.2 cpd; M2 +-.2 cpd; M3 +-.2 cpd; M4 +-.2 cpd;
     M5 +-.2 cpd; M6 +-.21 cpd; M7 (.26-.29 cpd); and M8 (.30-.50 cpd).
     S. Lentz  10/28/99
     R. Pawlowicz 11/1/00
     Version 1.0
     Define frequency bands for spectral averaging.
    """
    fband = np.array([0.0001, 0.00417, 0.03192, 0.04859, 0.07218, 0.08884, 0.11243, 0.1291, 0.15269, 0.16936, 0.19295, 0.20961, 0.2332, 0.251, 0.26, 0.29, 0.3, 0.5]).reshape(9,2)

    # If we have a sampling interval> 1 hour, we might have to get
    # rid of some bins.
    #fband(fband(:,1)>1/(2*dt),:)=[];
    nfband = fband.shape[0]
    nx = max(xres.shape)

    # Spectral estimate (takes real time series only).
    # Matlab has changed their spectral estimator functions
    # To match the old code, I have to divide by 2*dt. This is because
    #
    #  PSD*dt  is two-sided spectrum in units of power per hertz.
    #
    #  PWELCH is the one-sided spectrum in power per hertz
    #
    #  So PWELCH/2 = PSD*dt
    #[Pxr,fx]=psd(real(xres),nx,1/dt);
    # Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
    #[Pxi,fx]=psd(imag(xres),nx,1/dt);
    # Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
    #[Pxc,fx]=csd(real(xres),imag(xres),nx,1/dt);
    # Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
    fx, Pxr = sps.welch(np.real(xres), window=np.hanning(nx), noverlap=np.ceil(nx / 2),nfft=nx,fs=1/dt,nperseg=nx) # nargout=2
    # Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
    Pxr = Pxr / 2 / dt
    fx, Pxi = sps.welch(np.imag(xres), window=np.hanning(nx), noverlap=np.ceil(nx / 2),nfft=nx,fs=1/dt,nperseg=nx) # nargout=2
    # Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
    Pxi = Pxi / 2 / dt
    Pxc,fx = mplm.csd(np.real(xres), np.imag(xres),nx,1 / dt) # nargout=2
    # Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
    #matlab cpsd returns only reals when given a real xres have to test for complex and maybe change to ifstatement
    Pxc=np.real(Pxc)
    Pxc = Pxc / 2 / dt
    df = fx[2] - fx[1]


    Pxr[np.around(fu / df).astype(int) ] = np.nan
    # Sets Px=NaN in bins close to analyzed frequencies
    Pxi[np.around(fu / df).astype(int)] = np.nan
    # (to prevent leakage problems?).
    Pxc[np.around(fu / df).astype(int)] = np.nan
    Pxrave = np.zeros(shape=(nfband, 1), dtype='float64')
    Pxiave = np.zeros(shape=(nfband, 1), dtype='float64')
    Pxcave = np.zeros(shape=(nfband, 1), dtype='float64')
    # Loop downwards in frequency through bands (cures short time series
    # problem with no data in lowest band).
    #
    # Divide by nx to get power per frequency bin, and multiply by 2
    # to account for positive and negative frequencies.
    #

    for k in range(nfband-1, -1, - 1):
        jband = np.flatnonzero(np.all(np.vstack([fx >= fband[(k), 0],fx <= fband[(k), 1] , np.isfinite(Pxr)]).T,axis=1))
        if any(jband):
            Pxrave[k] = np.dot(np.mean(Pxr[(jband)]), 2) / nx
            Pxiave[k] = np.dot(np.mean(Pxi[(jband)]), 2) / nx
            Pxcave[k] = np.dot(np.mean(Pxc[(jband)]), 2) / nx
        else:
            if k < nfband:
                Pxrave[k] = Pxrave[(k + 1)]
                # Low frequency bin might not have any points...
                Pxiave[k] = Pxiave[(k + 1)]
                Pxcave[k] = Pxcave[(k + 1)]

    return fband, Pxrave, Pxiave, Pxcave


def noise_stats(xres, fu, dt):
    """NOISE_STATS: Computes statistics of residual energy for all
     constituents (ignoring any cross-correlations between real and
     imaginary parts).
     S. Lentz  10/28/99
     R. Pawlowicz 11/1/00
     Version 1.0
    """
    fband, Pxrave, Pxiave, Pxcave = residual_spectrum(xres, fu, dt) # nargout=4
    nfband = fband.shape[0]
    mu = max(fu.shape)
    # Get the statistics for each component.
    ercx = np.zeros(shape=(mu, 1), dtype='float64')
    eicx = np.zeros(shape=(mu, 1), dtype='float64')
    for k1 in range(0, nfband):
        k = np.flatnonzero(np.all([fu >= fband[k1, 0] , fu <= fband[k1,1]],axis=0))
        ercx[(k)] = np.sqrt(Pxrave[(k1)])
        eicx[(k)] = np.sqrt(Pxiave[(k1)])
    return ercx, eicx

def errell(cxi, sxi, ercx, ersx, ercy, ersy):
    """[emaj,emin,einc,epha]=errell(cx,sx,cy,sy,ercx,ersx,ercy,ersy) computes
     the uncertainities in the ellipse parameters based on the
     uncertainities in the least square fit cos,sin coefficients.

      INPUT:  cx,sx=cos,sin coefficients for x
              cy,sy=cos,sin coefficients for y
              ercx,ersx=errors in x cos,sin coefficients
              ercy,ersy=errors in y cos,sin coefficients

      OUTPUT: emaj=major axis error
              emin=minor axis error
              einc=inclination error (deg)
              epha=pha error (deg)
     based on linear error propagation, with errors in the coefficients
     cx,sx,cy,sy uncorrelated.
     B. Beardsley  1/15/99; 1/20/99
     Version 1.0
    """
    r2d = 180.0 / np.pi
    cx = np.real(cxi[:])
    sx = np.real(sxi[:])
    cy = np.imag(cxi[:])
    sy = np.imag(sxi[:])
    ercx = ercx[:]
    ersx = ersx[:]
    ercy = ercy[:]
    ersy = ersy[:]
    rp = 0.5 * np.sqrt((cx + sy) ** 2 + (cy - sx) ** 2)
    rm = 0.5 * np.sqrt((cx - sy) ** 2 + (cy + sx) ** 2)
    ercx2 = ercx ** 2
    ersx2 = ersx ** 2
    ercy2 = ercy ** 2
    ersy2 = ersy ** 2
    # major axis error
    ex = (cx + sy) / rp
    fx = (cx - sy) / rm
    gx = (sx - cy) / rp
    hx = (sx + cy) / rm
    dcx2 = (0.25 * (ex + fx)) ** 2
    dsx2 = (0.25 * (gx + hx)) ** 2
    dcy2 = (0.25 * (hx - gx)) ** 2
    dsy2 = (0.25 * (ex - fx)) ** 2
    emaj = np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 + dcy2 * ercy2 + dsy2 * ersy2)
    # minor axis error
    dcx2 = (0.25 * (ex - fx)) ** 2
    dsx2 = (0.25 * (gx - hx)) ** 2
    dcy2 = (0.25 * (hx + gx)) ** 2
    dsy2 = (0.25 * (ex + fx)) ** 2
    emin = np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 + dcy2 * ercy2 + dsy2 * ersy2)
    # inclination error
    rn = np.dot(2.0, (cx * cy + sx * sy))
    rd = cx ** 2 + sx ** 2 - (cy ** 2 + sy ** 2)
    den = rn ** 2 + rd ** 2
    dcx2 = ((rd * cy - rn * cx) / den) ** 2
    dsx2 = ((rd * sy - rn * sx) / den) ** 2
    dcy2 = ((rd * cx + rn * cy) / den) ** 2
    dsy2 = ((rd * sx + rn * sy) / den) ** 2
    einc = r2d * np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 + dcy2 * ercy2 + dsy2 * ersy2)
    # phase error
    rn = np.dot(2.0, (cx * sx + cy * sy))
    rd = cx ** 2 - sx ** 2 + cy ** 2 - sy ** 2
    den = rn ** 2 + rd ** 2
    dcx2 = ((rd * sx - rn * cx) / den) ** 2
    dsx2 = ((rd * cx + rn * sx) / den) ** 2
    dcy2 = ((rd * sy - rn * cy) / den) ** 2
    dsy2 = ((rd * cy + rn * sy) / den) ** 2
    epha = r2d * np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 + dcy2 * ercy2 + dsy2 * ersy2)
    return emaj, emin, einc, epha


def style_check(style):

    if (style=='pandas'):
        try:
            import pandas as pd
        except ImportError:
            print('Cannot import pandas using classic output.')
            style='classic'


    return style


def classic_style(out):
        print('-----------------------------------')

        print('nobs = %d \nngood = %d \nrecord length (days) = %.2f' % (out['nobs'], out['ngood'], np.dot(max(out['xin'].shape), out['dt']) / 24))
        if ('stime' in out):
            print('start time: %s' % dates.num2date(out['stime']).strftime('%Y-%m-%d %H:%M:%S'))
        print('rayleigh criterion = %.1f\n' % out['ray'])
        print('%s' % out['nodcor'])
        print('x0= %.3g, x trend= %.3g' % ( np.real(out['z0']), np.real(out['dz0'])))
        print('var(x)= ' , out['varx'] , '   var(xp)= ' , out['varxp'] , '   var(xres)= ' , out['varxr'] , '')
        print()
        print('percent var predicted/var original= %.1f ' % (np.dot(100, out['varxp']) / out['varx']))
        print()
        if (not 'complex' in out['xin'].dtype.name):
            print('        tidal amplitude and phase with 95 % CI estimates')
            print(' tide      freq        amp      amp_err    pha      pha_err    snr')
            for k,fuk in enumerate(out['fu']):
                outstr=(out['nameu'][k].decode(sys.stdout.encoding), fuk, out['tidecon'][k,0], out['tidecon'][k,1], out['tidecon'][k,2], out['tidecon'][k,3], out['snr'][k])
                if out['snr'][k] > out['synth']:
                    print('* %s  %9.7f  %9.4f  %8.3f  %8.2f  %8.2f  %8.2g' % outstr)
                else:
                    print('  %s  %9.7f  %9.4f  %8.3f  %8.2f  %8.2f  %8.2g' % outstr)
        else:
            print('y0= %.3f, x trend= %.3f' % (np.imag(out['z0']), np.imag(out['dz0'])))
            print('var(y)= %f    var(yp)= %f  var(yres)= %f ' % (out['vary'],out['varyp'],out['varyr']))
            print()
            print('percent var predicted/var original=  %.1f  ' %( np.dot(100, out['varyp']) / out['vary']))
            print('total var= %f   pred var=  %f ' % (out['varx'] + out['vary'],out['varxp'] + out['varyp']))
            print('percent total var predicted/var original=  %.1f  ' % ( np.dot(100, (out['varxp'] + out['varyp'])) / (out['varx'] + out['vary'])))
            print()
            print('ellipse parameters with 95 % CI estimates')
            print(' tide     freq        major      emaj      minor      emin     inc      einc      pha       epha      snr')
            for k,fuk in enumerate(out['fu']):
                outstr=(out['nameu'][k].decode(sys.stdout.encoding), fuk, out['tidecon'][k,0], out['tidecon'][k,1], out['tidecon'][k,2], out['tidecon'][k,3], out['tidecon'][k,4], out['tidecon'][k,5], out['tidecon'][k,6], out['tidecon'][k,7], out['snr'][k])
                if out['snr'][k] > out['synth']:
                    print('* %s  %9.7f  %9.4f  %8.3f %9.4f  %8.3f %8.2f  %8.2f  %8.2f  %8.2f %8.2g' % outstr)
                else:
                    print('  %s  %9.7f  %9.4f  %8.3f %9.4f  %8.3f %8.2f  %8.2f  %8.2f  %8.2f %8.2g' % outstr)


def pandas_style(out):
    print
    if (not 'complex' in out['xin'].dtype.name):
        print('======================================================================')
    else:
        print('====================================================================================================================')

    print('Number of observations = %d \nNumber of observations used = %d \nRecord length (days) = %.2f' % (out['nobs'], out['ngood'], np.dot(max(out['xin'].shape), out['dt']) / 24))
    if ('stime' in out):
        print('Start time: %s\n' % dates.num2date(out['stime']).strftime('%Y-%m-%d %H:%M:%S'))
    print('%s\n' % out['nodcor'])
    print('x0= %.3g, xtrend= %.3g' % ( np.real(out['z0']), np.real(out['dz0'])))
    print('var(data)= %.2f   var(prediction)= %.2f   var(residual)= %.2f'% (out['varx'],out['varxp'],out['varxr']))
    print('var(prediction)/var(data) (%%) = %.1f\n' %( 100*(out['varxp']/out['varx'])))

    import pandas as pd
    if (out['xin'].dtype!=complex):
        df=pd.DataFrame(np.vstack([out['fu'],out['tidecon'].T,out['snr']]).T,np.array( [name.decode(sys.stdout.encoding) for name in out['nameu']]),['Freq','Amp','Amp Err','Phase','Phase Err','SNR'])
        df.index.name='Tide'
        #df=df.sort('Amp')
        fmt={'Freq':'{:,.6f}'.format,'Amp':'{:,.2f}'.format,'Amp Err':'{:,.2f}'.format,'Phase':'{:,.1f}'.format,'Phase Err':'{:,.1f}'.format,'SNR':'{:,.1g}'.format}
        print('            Tidal amplitude and phase with 95 % CI estimates')
        print(df.to_string(col_space=10,formatters=fmt))

    else:
        print('y0= %.3g, ytrend= %.3g' % ( np.imag(out['z0']), np.imag(out['dz0'])))
        print('var(data)= %.2f   var(prediction)= %.2f   var(residual)= %.2f'% (out['vary'],out['varyp'],out['varyr']))
        print('var(prediction)/var(data) (%%) = %.1f\n' %( 100*(out['varyp']/out['vary'])))

        df=pd.DataFrame(np.vstack([out['fu'],out['tidecon'].T,out['snr']]).T,np.array( [name.decode(sys.stdout.encoding) for name in out['nameu']]),['Freq','Major','Major Err','Minor','Minor Err','Inc','Inc Err','Phase','Phase Err','SNR'])
        df.index.name='Tide'
        #df=df.sort('Major')
        fmt={'Freq':'{:,.6f}'.format,'Major':'{:,.2f}'.format,'Major Err':'{:,.2f}'.format,'Minor':'{:,.2f}'.format,'Minor Err':'{:,.2f}'.format,'Inc':'{:,.1f}'.format,'Inc Err':'{:,.1f}'.format,'Phase':'{:,.1f}'.format,'Phase Err':'{:,.1f}'.format,'SNR':'{:,.1g}'.format}
        print('                                   Ellipse parameters with 95 % CI estimates')

        print(df.to_string(col_space=10,formatters=fmt))

    print
    if (not 'complex' in out['xin'].dtype.name):
        print('======================================================================')
    else:
        print('====================================================================================================================')
















