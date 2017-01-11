import numpy as np
import os

try:
    testdir = os.path.dirname(os.path.realpath(__file__)) + '/'
except NameError:
    testdir = 'tests/'

t = np.arange(1001)
m2_freq = 2 * np.pi / 12.42
shallow_freq = 2 * np.pi / 2.484120261

# Build an M2 tidal series with 5m elevation
ein = 5 * np.cos(m2_freq * t)
sin = ein + 1 * np.cos(shallow_freq * t)

# Build an M2 tidal series for current
# Use to 2m/s horizontal currents and .5m/s vertical currents
uin = 2 * np.cos(m2_freq * t)
vin = 0.5 * np.cos(m2_freq * t)
uvin = uin + 1j * vin


# These are (arguments, fname)
cases = [
    ########################################################################
    # Elevation Test Cases
    ########################################################################

    # Basic test case
    (dict(xin=ein), 'NoArgs.elev'),
    # No output test case
    (dict(xin=ein, out_style=None), 'NoOutput.elev'),
    # Pandas output case
    (dict(xin=ein, out_style='pandas'), 'PandasOut.elev'),
    # M2 tides only
    (dict(xin=ein, constitnames=['M2']), 'M2only.elev'),
    # 5 tidal constituents (all should basically be zero other then M2)
    (dict(xin=ein, constitnames=['M2', 'N2', 'S2', 'K1', 'O1']), '5constit.elev'),
    # M2 tides with a starttime
    (dict(xin=ein, constitnames=['M2'], stime=768000), 'M2only-Stime.elev'),
    # M2 tides with a starttime and latitude
    (dict(xin=ein, constitnames=['M2'], stime=768000, lat=45), 'M2only-Stime-lat.elev'),

    ########################################################################
    # Shallow Test Cases
    ########################################################################

    # M2 tides and shallow
    (dict(xin=sin, constitnames=['M2'], shallownames=['M10']), 'M2.shallowM10'),
    # M2 tides and shallow and start_time
    (dict(xin=sin, constitnames=['M2'], stime=768000, shallownames=['M10']), 'M2-Stime.shallowM10'),

    ########################################################################
    # Current Test Cases
    ########################################################################

    # Basic test case
    (dict(xin=uvin), 'NoArgs.vel'),
    # No output test case
    (dict(xin=uvin, out_style=None), 'NoOutput.vel'),
    # Pandas output case
    (dict(xin=uvin, out_style='pandas'), 'PandasOut.vel'),
    # M2 tides only
    (dict(xin=uvin, constitnames=['M2']), 'M2only.vel'),
    # 5 tidal constituents (all should basically be zero other then M2)
    (dict(xin=uvin, constitnames=['M2', 'N2', 'S2', 'K1', 'O1']), '5constit.vel'),
    # M2 tides with a starttime
    (dict(xin=uvin, constitnames=['M2'], stime=768000), 'M2only-Stime.vel'),
    # M2 tides with a starttime and latitude
    (dict(xin=uvin, constitnames=['M2'], stime=768000, lat=45), 'M2only-Stime-lat.vel'),
]
