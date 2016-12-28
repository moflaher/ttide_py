from __future__ import division
import numpy as np
from . import time


def t_astron(jd):
    """T_ASTRON Computes astronomical Variables
     [A,ADER] = ASTRON(JD) computes the astronomical variables
                A=[tau,s,h,p,np,pp] (cycles)
      and their time derivatives
                ADER=[dtau,ds,dh,dp,dnp,dpp] (cycles/day)
      at the matlab time JD (UTC, but see code for details) where

        tau = lunar time
        s = mean longitude of the moon
        h = mean longitude of the sun
        p = mean longitude of the lunar perigee
        np = negative of the longitude of the mean ascending node
        pp = mean longitude of the perihelion (solar perigee)


        The formulae for calculating these ephemerides (other than tau)
        were taken from pages 98 and 107 of the Explanatory Supplement to
        the Astronomical Ephemeris and the American Ephemeris and Nautical
        Almanac (1961). They require EPHEMERIS TIME (ET), now TERRESTRIAL
        TIME (TT) and are based on observations made in the 1700/1800s.
        In a bizarre twist, the current definition of time is derived
        by reducing observations of planetary motions using these formulas.

        The current world master clock is INTERNATIONAL ATOMIC TIME (TAI).
        The length of the second is based on inverting the actual
        locations of the planets over the period 1956-65 into "time"
        using these formulas, and an offset added to keep the scale
        continuous with previous defns. Thus

                         TT = TAI + 32.184 seconds.

        Universal Time UT is a time scale that is 00:00 at midnight (i.e.,
        based on the earth's rotation rather than on planetary motions).
        Coordinated Universal Time (UTC) is kept by atomic clocks, the
        length of the second is the same as for TAI but leap seconds are
        inserted at intervals so that it provides UT to within 1 second.
        This is necessary because the period of the earth's rotation is
        slowly increasing (the day was exactly 86400 seconds around 1820,
        it is now about 2 ms longer). 22 leap seconds have been added in
        the last 27 years.

        As of 1/1/99,    TAI = UTC + 32 seconds.

        Thus,             TT = UTC + 62.184 seconds

        GPS time was synchronized with UTC 6/1/1980 ( = TAI - 19 secs),
        but is NOT adjusted for leap seconds. Your receiver might do this
        automatically...or it might not.

        Does any of this matter? The moon longitude is the fastest changing
        parameter at 13 deg/day. A time error of one minute implies a
        position error of less than 0.01 deg. This would almost always be
        unimportant for tidal work.

        The lunar time (tau) calculation requires UT as a base.  UTC is
        close enough - an error of 1 second, the biggest difference that
        can occur between UT and UTC, implies a Greenwich phase error of
        0.01 deg.  In Doodson's definition (Proc R. Soc. A, vol 100,
        reprinted in International Hydrographic Review, Appendix to
        Circular Letter 4-H, 1954) mean lunar time is taken to begin at
        "lunar midnight".
     Compute number of days from epoch of 12:00 UT Dec 31, 1899.
     (January 0.5 1900 ET)
    """
    # ## Matlab version info
    # B. Beardsley  12/29/98, 1/11/98
    # R. Pawlowicz  9/1/01
    # Version 1.0

    d = jd - time.date2num(time.datetime(1899, 12, 31, 12, 0, 0))
    D = d / 10000

    # Compute astronomical constants at time d1.
    args = np.array([1, d, D * D, D ** 3])

    # These are the coefficients of the formulas in the Explan. Suppl.
    sc = np.array([270.434164, 13.1763965268, - 8.5e-05, 3.9e-08])
    hc = np.array([279.696678, 0.9856473354, 2.267e-05, 0.0])
    pc = np.array([334.329556, 0.1114040803, - 0.0007739, - 2.6e-07])
    npc = np.array([- 259.183275, 0.0529539222, - 0.0001557, - 5e-08])
    #  first coeff was 281.220833 in Foreman but Expl. Suppl. has 44.
    ppc = np.array([281.220844, 4.70684e-05, 3.39e-05, 7e-08])
    coef = np.vstack([sc, hc, pc, npc, ppc])

    # Compute the parameters;
    # we only need the factional part of the cycle.
    astro = np.fmod(np.dot(coef, args) / 360.0, 1)

    # Compute lunar time tau, based on fractional part of solar day.
    # We add the hour angle to the longitude of the sun and
    # subtract the longitude of the moon.
    tau = (np.fmod(jd, 1) + astro[1] - astro[0])
    astro = np.hstack([tau, astro])

    # Compute rates of change.
    dargs = np.array([0, 1, 0.0002 * D, 0.0003 * D * D])

    ader = np.dot(coef, dargs) / 360.0
    dtau = (1.0 + ader[1] - ader[0])
    ader = np.hstack([dtau, ader])

    return astro, ader
