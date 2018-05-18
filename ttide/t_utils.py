from __future__ import division, print_function
import numpy as np
import scipy.signal as sps
import sys
from .t_getconsts import t_getconsts
from . import time

try:
    import pandas as pd
except ImportError:
    pd = None


np.set_printoptions(precision=8, suppress=True)

try:
    enc = sys.stdout.encoding
except AttributeError:
    enc = 'UTF-8'

if enc is None:
    # This is for piping/writing stdout to file, e.g. using '>'
    enc = 'UTF-8'


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
    if minres > 1 / (18.6 * 365.25 * 24):
        # Choose only resolveable pairs for short
        const, sat, cshallow = t_getconsts(centraltime)
        # Time series
        ju = np.flatnonzero(const['df'] >= minres)
    else:
        # Choose them all if > 18.6 years.
        const, sat, cshallow = t_get18consts(centraltime)

        ju = np.array([range(2,
                             (max(const['freq'].shape) + 1))]).reshape(1, -1).T
        # Skip Z0
        for ff in range(1, 3):
            # loop twice to make sure of neightbouring pairs
            jck = np.flatnonzero(np.diff(const['freq'][ju]) < minres)
            if (max(jck.shape) > 0):
                jrm = jck
                jrm = jrm + (abs(const['doodsonamp'][ju[(jck + 1 - 1)]]) <
                             abs(const['doodsonamp'][ju[(jck - 1)]]))
# disp('  Warning! Following constituent pairs violate Rayleigh criterion')
#               for ick in range(1, (max(jck.shape) +1)):
#                    disp('     ' + const.name(ju[(jck[(ick -1)] -1)], :)
# + ' vs ' + const.name(ju[(jck[(ick -1)] + 1 -1)], :) + ' - not using ' +
# const.name(ju[(jrm[(ick -1)] -1)], :))
                ju[(jrm - 1)] = np.array([])

    if constit.size:
        # Selected if constituents are specified in input.
        ju = np.array([], dtype=int)
        for k in range(0, (constit.shape[0])):
            j1 = np.where(const['name'] == constit[k])[0]
            if (j1.size == 0):
                print("Can't recognize name " +
                      constit[k].decode(enc) +
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
        for k in range(0, (shallow.shape[0])):
            j1 = np.where(const['name'] == shallow[k])[0]
            if (j1.size == 0):
                print("Can't recognize name " +
                      shallow[k].decode(enc) +
                      ' for forced search')
            else:
                ju = np.append(ju, j1)

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
        jinf = np.zeros(shape=(infname.shape[0], 1), dtype='float64') + np.NaN
        jref = np.zeros(shape=(infname.shape[0], 1), dtype='float64') + np.NaN
        for k in range(1, (infname.shape[0] + 1)):
            j1 = strmatch(infname[(k - 1), :], const.name)
            if (0 in j1.shape):
                print("Can't recognize name" +
                      infname[(k - 1), :] + ' for inference')
            else:
                jinf[(k - 1)] = j1
                fi[(k - 1)] = const['freq'][j1]
                namei[(k - 1), :] = const['name'][j1, :]
                j1 = strmatch(infref[(k - 1), :], nameu)
                if (0 in j1.shape):
                    print("Can't recognize name " + infref[(k - 1), :] +
                          ' for as a reference for inference')
                else:
                    jref[(k - 1)] = j1
                    print('   Inference of ' + namei[(k - 1), :] +
                          ' using ' + nameu[(j1 - 1), :] + '\\n')
        jinf[(isnan(jref) - 1)] = np.NaN
    return nameu, fu, ju, namei, fi, jinf, jref


def fixgaps(x):
    """FIXGAPS: Linearly interpolates gaps in a time series
     YOUT=FIXGAPS(YIN) linearly interpolates over NaN in the input time
     series (may be complex), but ignores trailing and leading NaNs.
     R. Pawlowicz 11/6/99
     Version 1.0
    """
    

    #find nans
    bd = np.isnan(x)

    #early exit if there are no nans  
    if not bd.any():
        return x
    
    #find nonnans index numbers
    gd = np.flatnonzero(~bd)

    #ignore leading and trailing nans
    bd[:gd.min()]=False
    bd[(gd.max()+1):]=False
    
    #interpolate nans
    x[bd] = np.interp(np.flatnonzero(bd),gd,x[gd])

    return x

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
    for k in range(1, (nfband + 1)):
        # The B matrix represents the covariance matrix for the vector
        # [Re{ap} Im{ap} Re{am} Im{am}]' where Re{} and Im{} are real and
        # imaginary parts, and ap/m represent the complex constituent
        # amplitudes for positive and negative frequencies when the input
        # is bivariate white noise. For a flat residual spectrum this works
        # fine.
        # This is adapted here for "locally white" conditions, but I'm still
        # not sure how to handle a complex sxy, so this is set to zero
        # right now.
        p = (Pxrave[(k - 1)] + Pxiave[(k - 1)]) / 2
        d = (Pxrave[(k - 1)] - Pxiave[(k - 1)]) / 2
        sxy = Pxcave[(k - 1)]
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

        Mat[:, :, (k - 1)] = np.dot(V, np.diag(np.sqrt(D)))
        # print Mat
    # Generate realizations for the different analyzed constituents.
    N = np.zeros(shape=(4, nreal), dtype='float64')
    NM = np.zeros(shape=(max(fu.shape), nreal), dtype='complex128')
    NP = np.zeros(shape=(max(fu.shape), nreal), dtype='complex128')

    for k in range(0, fu.shape[0]):
        l = np.squeeze(np.flatnonzero(np.all([fu[k] > fband[:, 0],
                                              fu[k] < fband[:, 1]], axis=0)))
        N = np.hstack([np.zeros(4,).reshape(-1, 1),
                       np.dot(np.squeeze(Mat[:, :, l]),
                              np.random.randn(4, nreal - 1))])
        NP[(k), :] = (N[0, :] + 1j * N[1, :])
        NM[(k), :] = (N[2, :] + 1j * N[3, :])

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
    fband = np.array([[0.0001, 0.00417],
                      [0.03192, 0.04859],
                      [0.07218, 0.08884],
                      [0.11243, 0.1291],
                      [0.15269, 0.16936],
                      [0.19295, 0.20961],
                      [0.2332, 0.251],
                      [0.26, 0.29],
                      [0.3, 0.5]])

    # If we have a sampling interval> 1 hour, we might have to get
    # rid of some bins.
    # fband(fband(:,1)>1/(2*dt),:)=[];
    nfband = fband.shape[0]
    nx = max(xres.shape)

    # Spectral estimate (takes real time series only).
    fx, Pxr = sps.welch(np.real(xres), window=np.hanning(nx),
                        noverlap=np.ceil(nx / 2), nfft=nx, fs=1 / dt, nperseg=nx)
    Pxr = Pxr / 2 / dt
    fx, Pxi = sps.welch(np.imag(xres), window=np.hanning(nx),
                        noverlap=np.ceil(nx / 2), nfft=nx, fs=1 / dt, nperseg=nx)
    Pxi = Pxi / 2 / dt
    #Pxc, fx = mplm.csd(np.real(xres), np.imag(xres), nx, 1 / dt)
    fx, Pxc = sps.csd(np.real(xres), np.imag(xres), fs=1 / dt, nperseg=nx, nfft=nx, )

    # matlab cpsd returns only reals when given a real xres have to
    # test for complex and maybe change to ifstatement
    Pxc = np.real(Pxc)
    Pxc = Pxc / 2 / dt
    df = fx[2] - fx[1]

    # Sets Px=NaN in bins close to analyzed frequencies
    # to prevent leakage problems?).
    Pxr[np.around(fu / df).astype(int)] = np.nan
    Pxi[np.around(fu / df).astype(int)] = np.nan
    Pxc[np.around(fu / df).astype(int)] = np.nan

    Pxrave = np.zeros(shape=(nfband, 1), dtype='float64')
    Pxiave = np.zeros(shape=(nfband, 1), dtype='float64')
    Pxcave = np.zeros(shape=(nfband, 1), dtype='float64')

    # Loop downwards in frequency through bands (cures short time series
    # problem with no data in lowest band).
    # Divide by nx to get power per frequency bin, and multiply by 2
    # to account for positive and negative frequencies.
    for k in range(nfband - 1, -1, - 1):
        jband = np.flatnonzero(np.all(np.vstack([fx >= fband[(k), 0],
                                                 fx <= fband[(k), 1],
                                                 np.isfinite(Pxr)]).T, axis=1))
        if any(jband):
            Pxrave[k] = 2 * np.mean(Pxr[(jband)]) / nx
            Pxiave[k] = 2 * np.mean(Pxi[(jband)]) / nx
            Pxcave[k] = 2 * np.mean(Pxc[(jband)]) / nx
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
    fband, Pxrave, Pxiave, Pxcave = residual_spectrum(xres, fu, dt)
    nfband = fband.shape[0]
    mu = max(fu.shape)
    # Get the statistics for each component.
    ercx = np.zeros(shape=(mu, 1), dtype='float64')
    eicx = np.zeros(shape=(mu, 1), dtype='float64')
    for k1 in range(0, nfband):
        k = np.flatnonzero(np.all([fu >= fband[k1, 0],
                                   fu <= fband[k1, 1]], axis=0))
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
    emaj = np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 +
                   dcy2 * ercy2 + dsy2 * ersy2)
    # minor axis error
    dcx2 = (0.25 * (ex - fx)) ** 2
    dsx2 = (0.25 * (gx - hx)) ** 2
    dcy2 = (0.25 * (hx + gx)) ** 2
    dsy2 = (0.25 * (ex + fx)) ** 2
    emin = np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 +
                   dcy2 * ercy2 + dsy2 * ersy2)
    # inclination error
    rn = 2.0 * (cx * cy + sx * sy)
    rd = cx ** 2 + sx ** 2 - (cy ** 2 + sy ** 2)
    den = rn ** 2 + rd ** 2
    dcx2 = ((rd * cy - rn * cx) / den) ** 2
    dsx2 = ((rd * sy - rn * sx) / den) ** 2
    dcy2 = ((rd * cx + rn * cy) / den) ** 2
    dsy2 = ((rd * sx + rn * sy) / den) ** 2
    einc = r2d * np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 +
                         dcy2 * ercy2 + dsy2 * ersy2)
    # phase error
    rn = 2.0 * (cx * sx + cy * sy)
    rd = cx ** 2 - sx ** 2 + cy ** 2 - sy ** 2
    den = rn ** 2 + rd ** 2
    dcx2 = ((rd * sx - rn * cx) / den) ** 2
    dsx2 = ((rd * cx + rn * sx) / den) ** 2
    dcy2 = ((rd * sy - rn * cy) / den) ** 2
    dsy2 = ((rd * cy + rn * sy) / den) ** 2
    epha = r2d * np.sqrt(dcx2 * ercx2 + dsx2 * ersx2 +
                         dcy2 * ercy2 + dsy2 * ersy2)

    return emaj, emin, einc, epha


def variance_str(out):
    '''Takes the out dictionary and prints the variance text'''

    x = np.var(out['xingd'].real, ddof=1)
    xp = np.var(out['xoutgd'].real, ddof=1)
    xr = np.var(out['xresgd'].real, ddof=1)
    z0r = out['z0'].real
    dz0r = out['dz0'].real

    outstr = 'x0= {:.3g}  xtrend= {:.3g}\n'.format(z0r, dz0r)
    outstr += ('var(data)= {:.2f}' + ' ' * 4 +
               'var(prediction)= {:.2f}' + ' ' * 4 +
               'var(residual)= {:.2f}\n').format(x, xp, xr)
    outstr += 'var(prediction)/var(data) (%%) = %.1f\n\n' % (100 * xp / x)

    if np.iscomplexobj(out['xin']):
        y = np.var(out['xingd'].imag, ddof=1)
        yp = np.var(out['xoutgd'].imag, ddof=1)
        yr = np.var(out['xresgd'].imag, ddof=1)
        z0r = out['z0'].imag
        dz0r = out['dz0'].imag

        outstr += 'y0= {:.3g}  ytrend= {:.3g}\n'.format(z0r, dz0r)
        outstr += ('var(data)= {:.2f}' + ' ' * 4 +
                   'var(prediction)= {:.2f}' + ' ' * 4 +
                   'var(residual)= {:.2f}\n').format(y, yp, yr)
        outstr += 'var(prediction)/var(data) (%) = {:.1f}\n\n'.format(
            100 * yp / y)

        outstr += 'total_var= {:f} pred_var=  {:f}\n'.format(
            (x + y), (xp + yp))
        outstr += 'total_var/pred_var (%) =  {:.1f}  \n'.format(
            100 * (xp + yp) / (x + y))
    return outstr


def classic_style(out):
    outstr = '-' * 35 + '\n'
    outstr += ('nobs = %d \nngood = %d \nrecord length (days) = %.2f' %
               (out['nobs'], out['ngood'],
                np.dot(max(out['xin'].shape), out['dt']) / 24)) + '\n'
    if ('stime' in out):
        outstr += 'start time: {}\n'.format(
            time.num2date(out['stime']).strftime('%Y-%m-%d %H:%M:%S'))
    outstr += ('rayleigh criterion = %.1f\n\n' % out['ray'])
    outstr += ('%s\n' % out['nodcor'])
    outstr += variance_str(out)

    if out['isComplex']:
        outstr += (' ' * 32 + 'ellipse parameters with 95 % CI estimates\n')
        outstr += (' tide     freq        major      emaj' +
                   '      minor      emin     inc      einc' +
                   '      pha       epha      snr\n')
        fmt = ('{star} {name}  {fuk:9.7f}  '
               '{c[0]:9.4f}  {c[1]:8.3f} {c[2]:9.4f}  {c[3]:8.3f} '
               '{c[4]:8.2f}  {c[5]:8.2f}  {c[6]:8.2f}  {c[7]:8.2f} '
               '{snr:8.2g}\n')
        for k, fuk in enumerate(out['fu']):
            if out['snr'][k] > out['synth']:
                star = '*'
            else:
                star = ' '
            outstr += fmt.format(star=star,
                                 name=out['nameu'][k].decode(enc),
                                 fuk=fuk,
                                 c=out['tidecon'][k],
                                 snr=out['snr'][k])
    else:
        outstr += '        tidal amplitude and phase with 95 % CI estimates\n'
        outstr += (' tide      freq        amp      amp_err' +
                   '   pha      pha_err    snr\n')
        fmt = ('{star} {name}  {fuk:9.7f}  '
               '{c[0]:9.4f}  {c[1]:8.3f}  {c[2]:8.2f}  {c[3]:8.2f}  '
               '{snr:8.2g}\n')
        for k, fuk in enumerate(out['fu']):
            if out['snr'][k] > out['synth']:
                star = '*'
            else:
                star = ' '
            outstr += fmt.format(star=star,
                                 name=out['nameu'][k].decode(enc),
                                 fuk=fuk,
                                 c=out['tidecon'][k],
                                 snr=out['snr'][k])
    return outstr


def pandas_style(out, dfTF=False):
    if pd is None:
        # Unable to import pandas.
        print("pandas is not available, falling back to out_style='classic'.")
        return classic_style(out)
    spacer = 70
    if out['isComplex']:
        spacer = 116
    outstr = '=' * spacer + '\n'

    outstr += ('Number of observations = %d' % out['nobs']) + '\n'
    outstr += ('Number of observations used = %d' % out['ngood']) + '\n'
    outstr += ('Record length (days) = %.2f' % (out['nobs'] * out['dt'] / 24.0)) + '\n'

    if ('stime' in out):
        outstr += ('Start time: %s\n' %
                   time.num2date(out['stime']).strftime('%Y-%m-%d %H:%M:%S')) + '\n'
    outstr += ('%s\n' % out['nodcor']) + '\n'
    outstr += variance_str(out)

    names = np.array([name.decode(enc) for name in out['nameu']])
    fmt = {'Freq': '{:,.6f}'.format, 'Major': '{:,.2f}'.format,
           'Major Err': '{:,.2f}'.format, 'Minor': '{:,.2f}'.format,
           'Minor Err': '{:,.2f}'.format, 'Inc': '{:,.1f}'.format,
           'Inc Err': '{:,.1f}'.format, 'Phase': '{:,.1f}'.format,
           'Phase Err': '{:,.1f}'.format, 'SNR': '{:,.1g}'.format,
           'Amp': '{:,.2f}'.format, 'Amp Err': '{:,.2f}'.format}

    if out['isComplex']:
        colnames = ['Freq', 'Major', 'Major Err', 'Minor', 'Minor Err',
                    'Inc', 'Inc Err', 'Phase', 'Phase Err', 'SNR']
        outstr += (' ' * 35 + 'Ellipse parameters with 95 % CI estimates') + '\n'
    else:
        colnames = ['Freq', 'Amp', 'Amp Err', 'Phase', 'Phase Err', 'SNR']
        outstr += (' ' * 12 + 'Tidal amplitude and phase with 95 % CI estimates') + '\n'

    dfdata = np.vstack([out['fu'], out['tidecon'].T, out['snr']]).T
    df = pd.DataFrame(dfdata, names, colnames)
    df.index.name = 'Tide'
    outstr += (df.to_string(col_space=10, formatters=fmt)) + '\n'

    outstr += ('=' * spacer) + '\n'
    if dfTF:
        return outstr, df
    else:
        return outstr
