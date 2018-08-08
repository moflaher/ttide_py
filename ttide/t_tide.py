from __future__ import division, print_function
import numpy as np
import scipy.interpolate as spi
from .t_vuf import t_vuf
from . import t_utils as tu
from .base import TTideCon, t_predic
import datetime
from . import time as tm


pi = np.pi


np.set_printoptions(precision=8, suppress=True)


def t_tide(xin, dt=1, stime=None, lat=None,
           out_style='classic',
           outfile=None,
           corr_fs=[0, 1e6], corr_fac=[1, 1],
           secular='mean',
           infiname=[], infirefname=[],
           ray=1,
           shallownames=[], constitnames=[],
           errcalc='cboot', synth=2,
           lsq='best'):
    """T_TIDE Harmonic analysis of a time series.

    Parameters
    ----------
    xin : array_like
       can be real (e.g. for elevations), or complex (U + 1j * V)
       for eastward velocity U and northward velocity V.

    dt : float
       Sampling interval in hours, default = 1.

    stime : float (mpl_datenum) or python datetime
       The start time of the series, in matplotlib_datenum format (default empty).

    lat : float
       decimal degrees (+north) (default: none).

    out_style : {None, 'classic', 'pandas'}
       where to send printed output
         None      - no printed output
         'classic' - to screen in classic-mode (default)
         'pandas'  - to screen in pandas-like mode.

    outfile : str or None
       The filename to write to (default: None, do not write to a
       file). This writes in format `out_style` ('classic' if
       `out_style` is None).

    corr_fs : array_like
        frequencies of the pre-filter transfer function (see note on
        pre-filtering)
    corr_fac : array_like (possible complex)
        correction factor magnitudes (see note on pre-filtering)

    secular : {'mean', 'linear'}
      Adjustment for long-term behavior ("secular" behavior).
           'mean'    - assume constant offset (default).
           'linear'  - get linear trend.

    infiname : list-like of names of consituents to be inferred
        Unclear what a more specific docstring should be. Clarify docs here.

    infirefname : list-like of names of references
        Unclear what a more specific docstring should be. Clarify docs here.

    ray : float
        The Rayleigh criteria (default: 1)

    shallownames : list-like of strings
        The names of shallow-water constituents to analyze.

    errcalc : string
        Method to use for calculation of confidence limits:
            'wboot' - Boostrapped confidence intervals based on a
                      correlated bivariate white-noise model.
            'cboot' - Boostrapped confidence intervals based on an
                      uncorrelated bivariate coloured-noise model
                      (default).'
           'linear' - Linearized error analysis that assumes an
                      uncorrelated bivariate coloured noise model.

    synth : float
        The signal-to-noise ratio of constituents to use for the
        "predicted" tide (passed to t_predic, but note that the
        default value is different).
            0 - use all selected constituents
           >0 - use only those constituents with a SNR greater than
                that given (1 or 2 are good choices, 2 is the
                default).
           <0 - return result of least-squares fit (should be the same
                as using '0', except that NaN-holes in original time
                series will remain and mean/trend are included).

    lsq : string
        'direct'  - use A\ x fit
        'normal'  - use (A'A)\(A'x) (may be necessary for very large
                  input vectors since A'A is much smaller than A)
          'best'  - automatically choose based on length of series
                  (default).

    Returns
    -------

    out : `TTideCon` instance
        This class is based on dict. The dictionary contains all of
        the relavent data of the fit. It also includes a `t_predic`
        method that can be used to create tidal predictions
        (extrapolation) based on the fit. This `t_predic` method is
        also duplicated as the `__call__` method. See the TTideCon
        docstring for more info.

    Notes
    -----

    This function is based on the Matlab T_TIDE toolbox by Rich Pawlowicz [1]_.

    `stime` and `lat` are required if nodal corrections are to be computed,
    otherwise not necessary. If they are not included then the reported
    phases are raw constituent phases at the central time.

    Currently, timeseries longer than 18.6 years are not
    supported. (The Matlab version does support this.)

    PRE-FILTERING
    .............
    If the time series has been passed through a pre-filter of some
    kind (say, to reduce the low-frequency variability), then the
    analyzed constituents will have to be corrected for this. The
    correction transfer function (1/filter transfer function) has
    (possibly complex) magnitude `corr_fac` at frequency `corr_fs`
    (cph).  Corrections of more than a factor of 100 are not applied;
    it is assumed these refer to tidal constituents that were
    intentionally filtered out, e.g., the fortnightly components.

    NaN-Values
    ..........
    Although missing data can be handled with NaN, it is wise not to
    have too many of them. If your time series has a lot of missing
    data at the beginning and/or end, then truncate the input time
    series.  The Rayleigh criterion is applied to frequency intervals
    calculated as the inverse of the input series length.

    .. [1] Pawlowicz, R., B. Beardsley, and S. Lentz, "Classical Tidal
       "Harmonic Analysis Including Error Estimates in MATLAB
       using T_TIDE", Computers and Geosciences, 28, 929-937 (2002).

    Examples
    --------

    import numpy as np
    import ttide as tt

    t = np.arange(1001)
    m2_freq = 2 * np.pi / 12.42

    ####
    # Here is an example 'real' dataset:
    elev = 5 * np.cos(m2_freq * t)

    # Compute the tidal fit:
    tfit_e = tt.t_tide(elev)

    # Construct the fitted time-series:
    elev_fit = tfit_e(t)

    # Or extrapolate the fit to other times:
    extrap_fit = tfit_e(np.arange(2000,2500))

    ####
    # Here is an example 'complex' (vector) dataset:
    vel = 0.8 * elev + 1j * 2 * np.sin(m2_freq * t)

    tfit_v = tt.t_tide(vel)

    """

    # ### This text needs to be re-included in the docstring once this
    # ### is supported.
    # If your time series is longer than 18.6 years then nodal corrections
    # are not made -instead we fit directly to all satellites (start time
    # is then just used to generate Greenwich phases).

    # ### This text came from the docstring on `infiname` and `infirefname`
    # Inference of constituents.
    #        'inference'      NAME,REFERENCE,AMPRAT,PHASE_OFFSET
    #                         where NAME is an array of the names of
    #                         constituents to be inferred, REFERENCE is an
    #                         array of the names of references, and AMPRAT
    #                         and PHASE_OFFSET are the amplitude factor and
    #                         phase offset (in degrees)from the references.
    #                         NAME and REFERENCE are Nx4 (max 4 characters
    #                         in name), and AMPRAT and PHASE_OFFSET are Nx1
    #                         (for scalar time series) and Nx2 for vector
    #                         time series (column 1 is for + frequencies and
    #                         column 2 for - frequencies).
    #                         NB - you can only infer ONE unknown constituent
    #                         per known constituent (i.e. REFERENCE must not
    #                         contain multiple instances of the same name).

    if isinstance(stime, (datetime.datetime, np.datetime64)):
        stime = tm.date2num(stime)

    corr_fs = np.array(corr_fs)
    corr_fac = np.array(corr_fac)
    infiname = np.array(infiname)
    infirefname = np.array(infirefname)
    shallownames = tu.fourpad(np.array(shallownames))
    constitnames = tu.fourpad(np.array(constitnames))

    isComplex = False

    # Check to make sure that incoming data is a vector.
    inn = xin.shape
    if len(inn) != 1:
        raise ValueError('Input time series is not a vector')
    if 'complex' in xin.dtype.name:
        isComplex = True

    # Check size of incoming data.
    nobs = max(xin.shape)

    # Set matrix method if auto-choice.
    # This could be increased for nobs>10000.
    # 100,000 uses a reasonable amount of ram for current systems.
    # Kept as is for the time being to stay true to Matlab version.
    if lsq[0:3] == 'bes':
        if nobs > 1000000:
            lsq = 'normal'
        else:
            lsq = 'direct'

    # Check to see if timeseries spans 18.6 years.
    if nobs * dt > 18.6 * 365.25 * 24:
        longseries = 1
        ltype = 'full'
    else:
        longseries = 0
        ltype = 'nodal'
    nobsu = nobs - np.remainder(nobs - 1, 2)

    # Make series odd to give a center point
    # Time vector for entire time series centered at series midpoint.
    t = (dt * (np.arange(nobs) + 1 - np.ceil(nobsu / 2)))

    if stime is not None:
        centraltime = stime + np.floor(nobsu / 2) * dt / 24.0
    else:
        centraltime = np.array([])

    # -------Get the frequencies to use in the harmonic analysis-----------
    tmptuple = tu.constituents(ray / (dt * nobsu),
                               constitnames, shallownames,
                               infiname, infirefname,
                               centraltime)
    nameu, fu, ju, namei, fi, jinf, jref = tmptuple

    mu = len(fu)
    mi = len(fi)

    # # inferred
    # Find the good data points
    # (here I assume that in a complex time series, if u is bad, so is v).
    gd = np.flatnonzero(np.isfinite(xin[0:nobsu]))
    ngood = max(gd.shape)

    # Now solve for the secular trend plus the analysis. Instead of solving
    # for + and - frequencies using exp(i*f*t), I use sines and cosines to
    # keep tc real.  If the input series is real, than this will
    # Automatically use real-only computation (faster).
    # However, for the analysis, it's handy to get the + and - frequencies
    # ('ap' and 'am'), and so that's what we do afterwards.
    # The basic code solves the matrix problem Ac=x+errors where the functions
    # to use in the fit fill up the A matrix, which is of size
    # (number points)x(number constituents). This can get very, very large
    # for long time series, and for this the more complex block processing
    # algorithm was added. It should give
    # identical results (up to roundoff error)
    if lsq[0:3] == 'dir':
        if secular[0:3] == 'lin':
            tc = np.hstack([np.ones((len(t), 1)),
                            np.cos(2 * pi * np.outer(t, fu)),
                            np.sin(2 * pi * np.outer(t, fu)),
                            t.reshape(-1, 1) * (2 / dt / nobsu)])
        else:
            tc = np.hstack([np.ones((len(t), 1)),
                            np.cos(2 * pi * np.outer(t, fu)),
                            np.sin(2 * pi * np.outer(t, fu))])

        coef = np.linalg.lstsq(tc[gd, :], xin[gd])[0].T

        # z0 a+ and a- amplitudes
        z0 = coef[0]
        ap = (coef[1:mu + 1] - 1j * coef[(1 + mu):(mu * 2) + 1]) / 2
        am = (coef[1:mu + 1] + 1j * coef[(1 + mu):(mu * 2) + 1]) / 2

        if secular[0:3] == 'lin':
            dz0 = coef[-1]
        else:
            dz0 = 0

        # Save least squares fitted prediction incase synth<=0
        xout = np.dot(tc, coef)

    else:
        # More complicated code required for long
        # timeseries when memory maybe a problem.
        # Modified from code submitted by Derek Goring (NIWA Chrischurch)
        # Basically the normal equations are formed
        # (rather than using Matlab's algorithm for least squares);
        # this can be done by adding up subblocks of data.
        # Notice how the code is messier,
        # and we have to recalculate everything to get the original fit.
        nsub = 5000
        # Block length - doesn't matter really but should be small enough to
        # get allocated quickly
        if secular[0:3] == 'lin':
            lhs = np.zeros(shape=(2 * mu + 2, 2 * mu + 2), dtype='float64')
            rhs = np.zeros(shape=(2 * mu + 2, ), dtype='float64')
            for j1 in range(1, (ngood + 1), nsub):
                j2 = np.min([j1 + nsub - 1, ngood])
                tslice = t[gd[(j1 - 1):j2] - 1]
                E = np.hstack([np.ones((j2 - j1 + 1, 1)),
                               np.cos(2 * pi * np.outer(tslice, fu)),
                               np.sin(2 * pi * np.outer(tslice, fu)),
                               tslice.reshape(-1, 1) * (2 / dt / nobsu)])
                rhs = rhs + np.dot(E.T, xin[(gd[(j1 - 1):j2] - 1)])
                lhs = lhs + np.dot(E.T, E)
        else:
            lhs = np.zeros(shape=(2 * mu + 1, 2 * mu + 1), dtype='float64')
            rhs = np.zeros(shape=(2 * mu + 1, ), dtype='float64')
            for j1 in range(1, (ngood + 1), nsub):
                j2 = np.min([j1 + nsub - 1, ngood])
                tslice = t[gd[(j1 - 1):j2] - 1]
                E = np.hstack([np.ones((j2 - j1 + 1, 1)),
                               np.cos(2 * pi * np.outer(tslice, fu)),
                               np.sin(2 * pi * np.outer(tslice, fu))])
                rhs = rhs + np.dot(E.T, xin[(gd[(j1 - 1):j2] - 1)])
                lhs = lhs + np.dot(E.T, E)
        coef = np.linalg.lstsq(lhs, rhs)[0].T

        # z0 a+ and a- amplitudes
        z0 = coef[0]
        ap = (coef[1:mu + 1] - 1j * coef[(1 + mu):(mu * 2) + 1]) / 2
        am = (coef[1:mu + 1] + 1j * coef[(1 + mu):(mu * 2) + 1]) / 2

        if secular[0:3] == 'lin':
            dz0 = coef[-1]
        else:
            dz0 = 0

        xout = xin.copy()
        # Copies over NaN
        if secular[0:3] == 'lin':
            for j1 in range(1, (nobs + 1), nsub):
                j2 = np.min([j1 + nsub - 1, nobs])
                tslice = t[(j1 - 1):j2]
                E = np.hstack([np.ones((j2 - j1 + 1, 1)),
                               np.cos(2 * pi * np.outer(tslice, fu)),
                               np.sin(2 * pi * np.outer(tslice, fu)),
                               np.dot(tslice, (2 / dt / nobsu)).reshape(-1, 1)])
                xout[(j1 - 1):j2] = np.dot(E, coef)
        else:
            for j1 in range(1, (nobs + 1), nsub):
                j2 = np.min([j1 + nsub - 1, nobs])
                tslice = t[(j1 - 1):j2]
                E = np.hstack([np.ones((j2 - j1 + 1, 1)),
                               np.cos(2 * pi * np.outer(tslice, fu)),
                               np.sin(2 * pi * np.outer(tslice, fu))])
                xout[(j1 - 1):j2] = np.dot(E, coef)

    # Check variance explained
    # (but do this with the original fit, and the residuals!)
    xres = xin - xout

    ####################################################################
    # ---------- Correct for prefiltering--------------------------------
    ####################################################################
    corrfac = spi.interpolate.interp1d(corr_fs, corr_fac)(fu)
    # To stop things blowing up!
    corrfac[corrfac > 100] = 1
    corrfac[corrfac < 0.01] = 1
    corrfac[np.isnan(corrfac)] = 1
    ap = ap * np.squeeze(corrfac)
    am = am * np.squeeze(np.conj(corrfac))

    ####################################################################
    # ---------------Nodal Corrections-----------------------------------
    # Generate nodal corrections and calculate phase relative to
    # Greenwich. Note that this is a slightly weird way to do the nodal
    # corrections, but is 'traditional'. The "right" way would be to
    # change the basis functions used in the least-squares fit above.
    ####################################################################
    if lat is not None and stime is not None:
        # Time and latitude
        # Get nodal corrections at midpoint time.
        v, u, f = t_vuf(ltype, centraltime,
                        np.hstack([ju, jinf]).astype(int), lat)
        vu = (v + u) * 360
        # total phase correction (degrees)
        nodcor = 'Greenwich phase computed with nodal\n \
                  corrections applied to amplitude\n \
                  and phase relative to center time\n'
    elif stime is not None:
        # Time only
        # Get nodal corrections at midpoint time
        v, u, f = t_vuf(ltype, centraltime,
                        np.hstack([ju, jinf]).astype(int))
        vu = (v + u) * 360
        # total phase correction (degrees)
        nodcor = 'Greenwich phase computed, no nodal corrections'
    else:
        # No time, no latitude
        nshape = (len(ju) + len(jinf), 1)
        vu = np.zeros(nshape, dtype='float64')
        f = np.ones(nshape, dtype='float64')
        nodcor = 'Phases at central time'

    ####################################################################
    # ---------------Inference Corrections------------------------------
    # Once again, the "right" way to do this
    # would be to change the basis functions.
    ####################################################################
    ii = np.flatnonzero(np.isfinite(jref))
    if ii:
        print('   Do inference corrections\\n')
        snarg = nobsu * pi * dt * (fi[(ii - 1)] - fu[(jref[(ii - 1)] - 1)])
        scarg = np.sin(snarg) / snarg
        if infamprat.shape[1] == 1:
            # For real time series
            pearg = np.dot(2 * pi,ii
                           (vu[(mu + ii - 1)] -
                            vu[(jref[(ii - 1)] - 1)] +
                            infph[(ii - 1)])) / 360
            pcfac = infamprat[(ii - 1)] * f[(mu + ii - 1)] / \
                f[(jref[(ii - 1)] - 1)] * np.exp(np.dot(ii, pearg))
            pcorr = 1 + pcfac * scarg
            mcfac = np.conj(pcfac)
            mcorr = np.conj(pcorr)
        else:
            # For complex time series
            pearg = np.dot(2 * pi,
                           (vu[(mu + ii - 1)] -
                            vu[(jref[(ii - 1)] - 1)] +
                            infph[(ii - 1), 0])) / 360
            pcfac = infamprat[(ii - 1), 0] * f[(mu + ii - 1)] / \
                f[(jref[(ii - 1)] - 1)] * np.exp(np.dot(i, pearg))
            pcorr = 1 + pcfac * scarg
            mearg = np.dot(-2 * pi,
                           (vu[(mu + ii - 1)] -
                            vu[(jref[(ii - 1)] - 1)] +
                            infph[(ii - 1), 1])) / 360
            mcfac = infamprat[(ii - 1), 1] * f[(mu + ii - 1)] / \
                f[(jref[(ii - 1)] - 1)] * np.exp(np.dot(i, mearg))
            mcorr = 1 + mcfac * scarg
        ap[(jref[(ii - 1)] - 1)] = ap[(jref[(ii - 1)] - 1)] / pcorr
        # Changes to existing constituents
        ap = np.array([ap, ap[(jref[(ii - 1)] - 1)] * pcfac]).reshape(1, -1)
        # Inferred constituents
        am[(jref[(ii - 1)] - 1)] = am[(jref[(ii - 1)] - 1)] / mcorr
        am = np.array([am, am[(jref[(ii - 1)] - 1)] * mcfac]).reshape(1, -1)
        fu = np.array([fu, fi[(ii - 1)]]).reshape(1, -1)
        nameu = np.array([nameu, namei[(ii - 1), :]]).reshape(1, -1)

    ####################################################################
    # --------------Error Bar Calculations------------------------------
    #
    # Error bar calcs involve two steps:
    # 1) Estimate the uncertainties in the analyzed amplitude
    #   for both + and - frequencies (i.e., in 'ap' and 'am').
    #   A simple way of doing this is to take the variance of the
    #   original time series and divide it into the amount appearing
    #   in the bandwidth of the analysis (approximately 1/length).
    #   A more sophisticated way is to assume "locally white"
    #   noise in the vicinity of, e.g., the diurnal consistuents.
    #   This takes into account slopes in the continuum spectrum.
    #
    # 2) Transform those uncertainties into ones suitable for ellipse
    #   parameters (axis lengths, angles). This can be done
    #   analytically for large signal-to-noise ratios. However, the
    #   transformation is non-linear at lows SNR, say, less than 10
    #   or so.
    #
    ####################################################################

    # Fill in "internal" NaNs with linearly interpolated
    # values so we can fft things.
    xr = tu.fixgaps(xres)

    nreal = 1

    if errcalc.endswith('boot'):
        # print('Using nonlinear bootstrapped error estimates.');
        ################################################################
        # "noise" matrices are created with the right covariance
        # structure to add to the analyzed components to
        # create 'nreal' REPLICATES.
        ################################################################

        nreal = 300
        # Create noise matrices
        NP, NM = tu.noise_realizations(xr[(np.isfinite(xr))],
                                       fu, dt, nreal, errcalc)
        # All replicates are then transformed (nonlinearly) into
        # ellipse parameters. The computed error bars are then
        # based on the std dev of the replicates.

        AP = np.repeat(ap, nreal).reshape(len(ap), nreal) + NP
        AM = np.repeat(am, nreal).reshape(len(am), nreal) + NM
        # Add to analysis (first column of NM,NP=0
        # so first column of AP/M holds ap/m).

        # Angle/magnitude form:
        epsp = np.angle(AP) * 180 / pi
        epsm = np.angle(AM) * 180 / pi

        ap = np.absolute(AP)
        am = np.absolute(AM)
    else:
        if errcalc == 'linear':
            print('Using linearized error estimates.')
            ############################################################
            # Uncertainties in analyzed amplitudes are computed in
            # different spectral bands. Real and imaginary parts of
            # the residual time series are treated separately
            # (no cross-covariance is assumed).
            #
            # Noise estimates are then determined from a linear analysis
            # of errors, assuming that everything is uncorrelated.
            # This is OK for scalar timeseries but can fail for vector
            # time series if the noise is not isotropic.
            ############################################################

            ercx, eicx = tu.noise_stats(xr[np.isfinite(xr)], fu, dt)

            # Note - here we assume that the error in the cos and sin
            # terms is equal, and equal to total power in the
            # encompassing frequency bin. It seems like there should be
            # a factor of 2 here somewhere but it only works this way!
            # <shrug>

            emaj, emin, einc, epha = errell(ap + am, 1j * (ap - am),
                                            ercx, ercx, eicx, eicx)
            epsp = 180 / np.pi * np.angle(ap)
            epsm = 180 / np.pi * np.angle(am)
            ap = np.absolute(ap)
            am = np.absolute(am)
        else:
            print("Unrecognized type of error analysis: " +
                  errcalc + " specified!")
    # -----Convert complex amplitudes to standard ellipse parameters----
    aap = ap / np.repeat(f, nreal).reshape(f.shape[0], nreal)
    # Apply nodal corrections and
    aam = am / np.repeat(f, nreal).reshape(f.shape[0], nreal)
    # compute ellipse parameters.
    fmaj = aap + aam
    # major axis
    fmin = aap - aam
    # minor axis

    gp = np.mod(np.repeat(vu, nreal).reshape(vu.shape[0], nreal) - epsp, 360)
    # pos. Greenwich phase in deg.
    gm = np.mod(np.repeat(vu, nreal).reshape(vu.shape[0], nreal) + epsm, 360)
    # neg. Greenwich phase in deg.
    finc = ((epsp + epsm) / 2)
    finc[:, 0] = np.mod(finc[:, 0], 180)

    # Ellipse inclination in degrees
    # (mod 180 to prevent ambiguity, i.e.,
    # we always ref. against northern
    # semi-major axis.
    finc = tu.cluster(finc, 180)
    # Cluster angles around the 'true' angle to avoid 360 degree wraps.
    pha = np.mod(gp + finc, 360)
    # Greenwich phase in degrees.
    pha = tu.cluster(pha, 360)
    # Cluster angles around the 'true' angle to avoid 360 degree wraps.

    # ----------------Generate 95% CI-----------------------------------
    # For bootstrapped errors, we now compute limits of the distribution.
    if errcalc.endswith('boot'):
        def booterrcalc(para, nreal):
            errval = np.multiply(
                np.median(
                    np.absolute(
                        para - (np.median(para, 1).reshape(-1, 1) *
                                np.ones([1, nreal]))), 1) / 0.6375, 1.96)
            return errval

        emaj = booterrcalc(fmaj, nreal)
        emin = booterrcalc(fmin, nreal)
        einc = booterrcalc(finc, nreal)
        epha = booterrcalc(pha, nreal)

    else:
        # In the linear analysis, the 95 CI are computed from the sigmas
        # by this fudge factor (infinite degrees of freedom).
        emaj = 1.96 * emaj
        emin = 1.96 * emin
        einc = 1.96 * einc
        epha = 1.96 * epha

    if isComplex:
        tidecon = np.array([fmaj[:, 0], emaj, fmin[:, 0], emin,
                            finc[:, 0], einc, pha[:, 0], epha]).T
    else:
        tidecon = np.array([fmaj[:, 0], emaj, pha[:, 0], epha]).T
    tideconout = tidecon.copy()
    # Sort results by frequency (needed if anything has been inferred
    # since these are stuck at the end of the list by code above).
    if any(np.isfinite(jref)):
        fu, I = np.sort(fu)
        nameu = nameu[(I - 1), :]
        tidecon = tidecon[(I - 1), :]
    snr = (tidecon[:, 0] / tidecon[:, 1]) ** 2
    # signal to noise ratio
    # --------Generate a 'prediction' using significant constituents----
    xoutOLD = xout
    if synth >= 0:
        if lat is not None and stime is not None:
            # This does not account for latitude,
            # functionality not added to t_predic yet.
            xout = t_predic(stime + np.array([range(nobs)]) * dt / 24.0,
                            nameu, fu, tidecon, synth=synth)
        elif stime is not None:
            xout = t_predic(stime + np.array([range(nobs)]) * dt / 24.0,
                            nameu, fu, tidecon, synth=synth)
        else:
            xout = t_predic(t / 24.0, nameu, fu, tidecon, synth=synth)

    # Check variance explained (but now do this
    # with the synthesized fit) and the residuals!
    xres = xin[:] - xout[:]

    xout = xout.reshape(inn[0], 1)

    out = TTideCon()
    out['nobs'] = nobs
    out['ngood'] = ngood
    out['dt'] = dt
    out['xin'] = xin
    out['xout'] = xout
    out['xres'] = xres
    out['xingd'] = xin[gd]
    out['xoutgd'] = xout[gd]
    out['xresgd'] = xres[gd]
    out['isComplex'] = isComplex
    out['ray'] = ray
    out['nodcor'] = nodcor
    out['z0'] = z0
    out['dz0'] = dz0

    out['fu'] = fu
    out['nameu'] = nameu
    out['tidecon'] = tideconout
    out['snr'] = snr
    out['synth'] = synth
    out['lat'] = lat
    out['ltype'] = ltype
    if stime is not None:
        out['stime'] = stime

    # -----------------Output results-----------------------------------
    if out_style or outfile:
        if out_style:
            method = out_style + '_style'
        else:
            method = 'classic_style'

        if outfile:
            getattr(out, method)(fname=outfile)
        else:
            print(getattr(out, method)(), end='')

    return out
