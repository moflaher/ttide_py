from __future__ import division, print_function
from cStringIO import StringIO
import sys
import numpy as np
from ttide.t_tide import t_tide
import os

testdir = os.path.dirname(os.path.realpath(__file__))

with open(testdir + '/data/testdata.out') as fl:
    tdata = fl.read().replace('\r\n', '\n')

# Initialize the pRNG so that the results will always be identical.
np.random.seed(29034230)

t = np.arange(1001)
m2_freq = 2 * np.pi / 12.42
shallow_freq = 2 * np.pi / 2.484120261

# Build an M2 tidal series with 5m elevation
ein = 5 * np.cos(m2_freq * t)

# Build an M2 tidal series for current
# Use to 2m/s horizontal currents and .5m/s vertical currents
uin = 2 * np.cos(m2_freq * t)
vin = 0.5 * np.cos(m2_freq * t)
uvin = uin + 1j * vin


def run_tests():

    # Redirect stdout to a StringIO instance
    sys.stdout = stdout = StringIO()

    ########################################################################
    # Elevation Test Cases
    ########################################################################

    # Basic test case
    [nameu, freq, tidecon, xout] = t_tide(ein)

    print()
    print('*' * 80)
    print()

    # No output test case
    [nameu, freq, tidecon, xout] = t_tide(ein, out_style=None)

    print()
    print('*' * 80)
    print()

    # Pandas output case
    [nameu, freq, tidecon, xout] = t_tide(ein, out_style='pandas')

    print()
    print('*' * 80)
    print()

    # M2 tides only
    [nameu, freq, tidecon, xout] = t_tide(ein, constitnames=['M2'])

    print()
    print('*' * 80)
    print()

    # 5 tidal constituents (all should basically be zero other then M2)
    [nameu, freq, tidecon, xout] = t_tide(ein, constitnames=['M2', 'N2', 'S2', 'K1', 'O1'])

    print()
    print('*' * 80)
    print()

    # M2 tides with a starttime
    [nameu, freq, tidecon, xout] = t_tide(ein, constitnames=['M2'], stime=768000)

    print()
    print('*' * 80)
    print()

    # M2 tides with a starttime and latitude
    [nameu, freq, tidecon, xout] = t_tide(ein, constitnames=['M2'], stime=768000, lat=45)

    print()
    print('*' * 80)
    print()

    # M2 tides and shallow
    # Build an M2 tidal series with 5m elevation
    sin = ein + 1 * np.cos(shallow_freq * t)
    [nameu, freq, tidecon, xout] = t_tide(sin,
                                          constitnames=['M2'], shallownames=['M10'])

    print()
    print('*' * 80)
    print()

    [nameu, freq, tidecon, xout] = t_tide(sin, stime=768000,
                                          constitnames=['M2'], shallownames=['M10'])

    print()
    print('*' * 80)
    print()

    ########################################################################
    # Current Test Cases
    ########################################################################

    # Basic test case
    [nameu, freq, tidecon, xout] = t_tide(uvin)

    print()
    print('*' * 80)
    print()

    # No output test case
    [nameu, freq, tidecon, xout] = t_tide(uvin, out_style=None)
    print('No output test case')

    print()
    print('*' * 80)
    print()

    # Pandas output case
    [nameu, freq, tidecon, xout] = t_tide(uvin, out_style='pandas')

    print()
    print('*' * 80)
    print()

    # M2 tides only
    [nameu, freq, tidecon, xout] = t_tide(uvin, constitnames=['M2'])

    print()
    print('*' * 80)
    print()

    # 5 tidal constituents (all should basically be zero other then M2)
    [nameu, freq, tidecon, xout] = t_tide(
        uvin, constitnames=['M2', 'N2', 'S2', 'K1', 'O1'])

    print()
    print('*' * 80)
    print()

    # M2 tides with a starttime
    [nameu, freq, tidecon, xout] = t_tide(
        uvin, constitnames=['M2'], stime=768000)

    print()
    print('*' * 80)
    print()

    # M2 tides with a starttime and latitude
    [nameu, freq, tidecon, xout] = t_tide(uvin, constitnames=['M2'], stime=768000, lat=45)

    # Reassign sys.stdout to original value.
    sys.stdout = sys.__stdout__

    assert tdata == stdout.getvalue()
