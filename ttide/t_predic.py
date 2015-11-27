from __future__ import division
import numpy as np
import os
from .t_getconsts import t_getconsts
from .t_vuf import t_vuf


def t_predic(tim, names, freq, tidecon, **kwargs):
    """T_PREDIC Tidal prediction
     YOUT=T_PREDIC(TIM,NAMES,FREQ,TIDECON) makes a tidal prediction
     using the output of T_TIDE at the specified times TIM in decimal
     days (from DATENUM). Optional arguments can be specified using
     property/value pairs:

           YOUT=T_PREDIC(...,TIDECON,property,value,...)

     Available properties are:

        In the simplest case, the tidal analysis was done without nodal
        corrections, and thus neither will the prediction. If nodal
        corrections were used in the analysis, then it is likely we will
        want to use them in the prediction too and these are computed
        using the latitude, if given.

         'latitude'        decimal degrees (+north) (default: none)

        If the original analysis was >18.6 years satellites are
        not included and we force that here:

         'anallength'      'nodal' (default)
                           'full'  For >18.6 years.

        The tidal prediction may be restricted to only some of the
        available constituents:

         'synthesis'    0 - Use all selected constituents.  (default)
                        scalar>0 - Use only those constituents with a SNR
                                   greater than that given (1 or 2 are
                                   good choices).


      It is possible to call t_predic without using property names, in
      which case the assumed calling sequence is

        YOUT=T_PREDIC(TIM,NAMES,FREQ,TIDECON,LATITUDE,SYNTHESIS);

      T_PREDIC can be called using the tidal structure available as an
      optional output from T_TIDE

        YOUT=T_PREDIC(TIM,TIDESTRUC,...)

      This is in fact the recommended calling procedure (and required
      when the analysis results are from series>18.6 years in length)
     R. Pawlowicz 11/8/99
     Version 1.0
     8/2/03 - Added block processing to generate prediction (to
              avoid memory overflows for long time series).
     29/9/04 - small bug with undefined ltype fixed
    """

    longseries = 0
    ltype = 'nodal'
    lat = np.array([])
    synth = 0
    k = 1
    tim = tim.reshape(-1, 1)

    # Use kwargs to set values other then the defaults
    if kwargs is not None:
        for key, value in kwargs.items():
            if (key == 'ltype'):
                ltype = value
            if (key == 'synth'):
                synth = value

    # Do the synthesis.
    snr = (tidecon[:, 0] / tidecon[:, 1]) ** 2
    # signal to noise ratio
    if synth > 0:
        I = snr > synth
        if not any(I):
            print('No predictions with this SNR')
            yout = np.nan + np.zeros(shape=(tim.shape, tim.shape),
                                     dtype='float64')
            return yout
        tidecon = tidecon[I, :]
        names = names[I]
        freq = freq[I]
    if tidecon.shape[1] == 4:
        # Real time series
        ap = np.multiply(tidecon[:, 0]/2.0,
                         np.exp(-1j*tidecon[:, 2]*np.pi/180))
        am = np.conj(ap)
    else:
        ap = np.multiply((tidecon[:, 0] + tidecon[:, 2]) / 2.0,
                         np.exp(np.dot(np.dot(1j, np.pi) / 180,
                                       (tidecon[:, 4] - tidecon[:, 6]))))

        am = np.multiply((tidecon[:, 0] - tidecon[:, 2]) / 2.0,
                         np.exp(np.dot(np.dot(1j, np.pi) / 180,
                                       (tidecon[:, 4] + tidecon[:, 6]))))

    # Mean at central point (get rid of one point at end to
    # take mean of odd number of points if necessary).
    jdmid = np.mean(tim[0:np.dot(2, np.fix((max(tim.shape) - 1) / 2)) + 1])
    if longseries:
        const = t_get18consts
        ju = np.zeros(shape=(freq.shape, freq.shape), dtype='float64')
        for k in range(1, (names.shape[0]+1)):
            inam = strmatch(names[(k-1), :], const.name)
            if max(inam.shape) == 1:
                ju[(k-1)] = inam
            else:
                if max(inam.shape) > 1:
                    minf, iminf = np.min(abs(freq[(k-1)] - const.freq(inam)))
                    ju[(k-1)] = inam[(iminf-1)]
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
    if ((lat.size != 0) & (np.absolute(jdmid) > 1)):
        v, u, f = t_vuf(ltype, jdmid, ju, lat)
    else:
        if np.fabs(jdmid) > 1:
            # a real date
            v, u, f = t_vuf(ltype, jdmid, ju)
        else:
            v = np.zeros((len(ju),), dtype='float64')
            u = v
            f = np.ones((len(ju),), dtype='float64')

    ap = ap * f * np.exp(+1j*2*np.pi*(u + v))
    am = am * f * np.exp(-1j*2*np.pi*(u + v))
    tim = tim - jdmid

    n, m = tim.shape
    ntim = max(tim.shape)
    nsub = 10000
    yout = np.zeros([n*m, ], dtype='complex128')

    # longer than one year hourly.
    for j1 in np.arange(0, ntim, nsub):
        j1 = j1.astype(int)
        j2 = np.min([j1 + nsub, ntim]).astype(int)
        tap = np.repeat(ap, j2-j1).reshape(len(ap), j2-j1)
        tam = np.repeat(am, j2-j1).reshape(len(am), j2-j1)

        touter = np.outer(24*1j*2*np.pi*freq, tim[j1:j2])
        yout[j1:j2] = np.sum(np.multiply(np.exp(touter), tap), axis=0) +\
            np.sum(np.multiply(np.exp(-touter), tam), axis=0)

    if (tidecon.shape[1] == 4):
        return np.real(yout)
    else:
        return yout
