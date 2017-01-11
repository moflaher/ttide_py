from __future__ import division
import numpy as np
import scipy as sp
from .t_astron import t_astron
from .t_getconsts import t_getconsts


def t_vuf(ltype, ctime, ju, lat=None):
    """T_VUF Computes nodal modulation corrections.
     [V,U,F]=T_VUF(TYPE,DATE,JU,LAT) returns the astronomical phase V, the
     nodal phase modulation U, and the nodal amplitude correction F at
     a decimal date DATE for the components specified by index JU
     at a latitude LAT.

     TYPE is either 'full' for the 18.6 year set of constitunets, or 'nodal'
     for the 1-year set with satellite modulations.

     If LAT is not specified, then the Greenwich phase V is computed with
     U=0 and F=1.

     Note that V and U are in 'cycles', not degrees or radians (i.e.,
     multiply by 360 to get degrees).

     If LAT is set to NaN, then the nodal corrections are computed for all
     satellites that do *not* have a "latitude-dependent" correction
     factor. This is for compatibility with the ways things are done in
     the xtide package. (The latitude-dependent corrections were zeroed
     out there partly because it was convenient, but this was rationalized
     by saying that since the forcing of tides can occur at latitudes
     other than where they are observed, the idea that observations have
     the equilibrium latitude-dependence is possibly bogus anyway).
     Get all the info about constituents.
     Calculate astronomical arguments at mid-point of data time series.
    """

    astro, ader = t_astron(ctime)

    if ltype == 'full':
        const = t_get18consts(ctime)
        # Phase relative to Greenwich (in units of cycles).
        v = rem(np.dot(const.doodson, astro) + const.semi, 1)
        v = v[(ju - 1)]
        u = np.zeros(shape=(v.shape, v.shape), dtype='float64')
        f = np.ones(shape=(v.shape, v.shape), dtype='float64')
    else:
        const, sat, shallow = t_getconsts(ctime)
        # Phase relative to Greenwich (in units of cycles).
        # (This only returns values when we have doodson#s,
        # i.e., not for the shallow water components,
        # but these will be computed later.)
        v = np.fmod(np.dot(const['doodson'], astro) + const['semi'], 1)

        if lat is not None:
            # If we have a latitude, get nodal corrections.
            # Apparently the second-order terms in the tidal potential
            # go to zero at the equator, but the third-order terms
            # do not. Hence when trying to infer the third-order terms
            # from the second-order terms, the nodal correction factors
            # blow up. In order to prevent this, it is assumed that the
            # equatorial forcing is due to second-order forcing OFF the
            # equator, from about the 5 degree location. Latitudes are
            # hence (somewhat arbitrarily) forced to be no closer than
            # 5 deg to the equator, as per note in Foreman.
            if abs(lat) < 5:
                lat = np.sign(lat) * 5
            slat = np.sin(np.pi * lat / 180)
            # Satellite amplitude ratio adjustment for latitude.
            rr = sat['amprat']
            # no amplitude correction
            if np.isfinite(lat):
                j = np.flatnonzero(sat['ilatfac'] == 1)
                # latitude correction for diurnal constituents
                rr[j] = rr[j] * 0.36309 * (1.0 - 5.0 * slat * slat) / slat
                j = np.flatnonzero(sat['ilatfac'] == 2)
                # latitude correction for semi-diurnal constituents
                rr[j] = rr[j] * 2.59808 * slat
            else:
                rr[sat['ilatfac'] > 0] = 0
            # Calculate nodal amplitude and phase corrections.
            uu = np.fmod(np.dot(sat['deldood'], astro.T[
                         3:6]) + sat['phcorr'], 1)
            # uu=uudbl-round(uudbl);  <_ I think this was wrong.
            # The original
            #                         FORTRAN code is:  IUU=UUDBL
            #                                           UU=UUDBL-IUU
            #                         which is truncation.
            # Sum up all of the satellite factors for all satellites.
            nsat = np.max(sat['iconst'].shape)
            nfreq = np.max(const['isat'].shape)

            fsum = np.array(1 + sp.sparse.csr_matrix(
                (np.squeeze(rr * np.exp(1j * 2 * np.pi * uu)),
                 (np.arange(0, nsat), np.squeeze(sat['iconst'] - 1))),
                shape=(nsat, nfreq)).sum(axis=0)).flatten()

            f = np.absolute(fsum)
            u = np.angle(fsum) / (2 * np.pi)

            # Compute amplitude and phase corrections
            # for shallow water constituents.
            for k in np.flatnonzero(np.isfinite(const['ishallow'])):
                ik = ((const['ishallow'][k] - 1 +
                       np.array(range(0, const['nshallow'][k]))).astype(int))
                iname = shallow['iname'][ik] - 1
                coef = shallow['coef'][ik]
                f[k] = np.prod(np.power(f[iname], coef))
                u[k] = np.dot(u[iname], coef)
                v[k] = np.dot(v[iname], coef)

            f = f[ju]
            u = u[ju]
            v = v[ju]

        else:
            # Astronomical arguments only, no nodal corrections.
            # Compute phases for shallow water constituents.
            for k in np.flatnonzero(np.isfinite(const['ishallow'])):
                ik = ((const['ishallow'][k] - 1 +
                       np.array(range(0, const['nshallow'][k]))).astype(int))
                v[k] = np.dot(v[shallow['iname'][ik] - 1], shallow['coef'][ik])

            v = v[ju]
            f = np.ones(len(v))
            u = np.zeros(len(v))

    return v, u, f
