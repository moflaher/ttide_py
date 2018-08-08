from __future__ import division
import numpy as np
from .t_getconsts import t_getconsts
from .t_vuf import t_vuf
from . import time as tm


def t_predic(t_time, names, freq, tidecon,
             lat=None, ltype='nodal', synth=0):
    """T_PREDIC Tidal prediction from tidal consituents.

    Parameters
    ----------
    time : array_like (N)
        The array of times at which to make predictions.
    names : array_like (M)
        The constituent name strings (e.g. ['M2', 'S2', ...])
    freq : array_like (M)
        The frequencies, in cph.
    tidecon : array_like (M, P)
        The tidal constituent amplitudes, phases, and orientations. P
        is either 4 (for real output), or 8 (for complex output).
    lat : flaot
        decimal degrees (+north) (default: None)
        In the simplest case, the tidal analysis was done without nodal
        corrections, and thus neither will the prediction. If nodal
        corrections were used in the analysis, then it is likely we will
        want to use them in the prediction too and these are computed
        using the latitude, if given.
    ltype : {'nodal', 'full'}
        If the original analysis was >18.6 years satellites are not
        included and we force that here. Use, 'full' if the
        consituents were computed from data longer that >18.6
        years. 'nodal' (default) is if the analysis was <18.6 years.
    synth : float
        The tidal prediction may be restricted to only some of the
        available constituents:
            0 - use all selected constituents (default)
           >0 - use only those constituents with a SNR greater than
                that given (1 or 2 are good choices).

    Returns
    -------
    yout : array (N,)
        The predicted time series (real/scaler for P = 4,
        complex/vector for P = 8).

    """

    longseries = 0  # Currently only timeseries <18.6 years are supported.
    if t_time.dtype.name.startswith('datetime64') or t_time.dtype is np.dtype("O"):
        t_time = tm.date2num(t_time)

    t_time = t_time.reshape(-1, 1)

    # Do the synthesis.
    snr = (tidecon[:, 0] / tidecon[:, 1]) ** 2
    # signal to noise ratio
    if synth > 0:
        I = snr > synth
        if not any(I):
            print('No predictions with this SNR')
            yout = np.nan + np.zeros(shape=(t_time.shape, t_time.shape),
                                     dtype='float64')
            return yout
        tidecon = tidecon[I, :]
        names = names[I]
        freq = freq[I]
    if tidecon.shape[1] == 4:
        # Real time series
        ap = np.multiply(tidecon[:, 0] / 2.0,
                         np.exp(-1j * tidecon[:, 2] * np.pi / 180))
        am = np.conj(ap)
    else:
        ap = np.multiply((tidecon[:, 0] + tidecon[:, 2]) / 2.0,
                         np.exp(1j * np.pi / 180 * (tidecon[:, 4] - tidecon[:, 6])))

        am = np.multiply((tidecon[:, 0] - tidecon[:, 2]) / 2.0,
                         np.exp(1j * np.pi / 180 * (tidecon[:, 4] + tidecon[:, 6])))

    # Mean at central point (get rid of one point at end to
    # take mean of odd number of points if necessary).
    jdmid = np.mean(t_time[0:((2 * int((max(t_time.shape) - 1) / 2)) + 1)])
    if longseries:
        const = t_get18consts
        ju = np.zeros(shape=(freq.shape, freq.shape), dtype='float64')
        for k in range(1, (names.shape[0] + 1)):
            inam = strmatch(names[(k - 1), :], const.name)
            if max(inam.shape) == 1:
                ju[(k - 1)] = inam
            else:
                if max(inam.shape) > 1:
                    minf, iminf = np.min(abs(freq[(k - 1)] - const.freq(inam)))
                    ju[(k - 1)] = inam[(iminf - 1)]
    else:
        const, sat, cshallow = t_getconsts(np.array([]))
        ju = np.zeros((len(freq),), dtype='int32')
        # Check to make sure names and frequencies match expected values.
        for k in range(0, (names.shape[0])):
            ju[k] = np.argwhere(const['name'] == names[(k)])
        # if any(freq~=const.freq(ju)),
        # error('Frequencies do not match names in input');
        # end;
    # Get the astronical argument with or without nodal corrections.
    if lat is not None and np.absolute(jdmid) > 1:
        v, u, f = t_vuf(ltype, jdmid, ju, lat)
    else:
        if np.fabs(jdmid) > 1:
            # a real date
            v, u, f = t_vuf(ltype, jdmid, ju)
        else:
            v = np.zeros((len(ju),), dtype='float64')
            u = v
            f = np.ones((len(ju),), dtype='float64')

    ap = ap * f * np.exp(+1j * 2 * np.pi * (u + v))
    am = am * f * np.exp(-1j * 2 * np.pi * (u + v))
    t_time = t_time - jdmid

    n, m = t_time.shape
    ntime = max(t_time.shape)
    nsub = 10000
    yout = np.zeros([n * m, ], dtype='complex128')

    # longer than one year hourly.
    for j1 in np.arange(0, ntime, nsub):
        j1 = j1.astype(int)
        j2 = np.min([j1 + nsub, ntime]).astype(int)
        tap = np.repeat(ap, j2 - j1).reshape(len(ap), j2 - j1)
        tam = np.repeat(am, j2 - j1).reshape(len(am), j2 - j1)

        touter = np.outer(24 * 1j * 2 * np.pi * freq, t_time[j1:j2])
        yout[j1:j2] = np.sum(np.multiply(np.exp(touter), tap), axis=0) +\
            np.sum(np.multiply(np.exp(-touter), tam), axis=0)

    if (tidecon.shape[1] == 4):
        return np.real(yout)
    else:
        return yout
