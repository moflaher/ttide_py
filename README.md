ttide_py
========

A direct conversion of T_Tide to Python.

This is a work in progress. It is not done.

It is now mostly functional.
Any help with finishing the conversion is welcome.

Credit for T\_Tide goes to Rich Pawlowicz, the original creator of T\_Tide. 
It is available at http://www2.ocgy.ubc.ca/~rich/.

A description of the theoretical basis of the analysis and some
implementation details of the Matlab version can be found in:

> Pawlowicz, R., B. Beardsley, and S. Lentz, "Classical Tidal
    "Harmonic Analysis Including Error Estimates in MATLAB
    using T_TIDE", Computers and Geosciences, 28, 929-937 (2002).

Citation of this article would be appreciated if you find the toolbox
useful (either the Matlab version, or this Python one).



Installation
============

This has little to no testing. Use at your own risk. To install, run:

    python setup.py install


Example Usage
=============

Imports and define some variables:

    import ttide as tt
    import numpy as np

    t = np.arange(1001)
    m2_freq = 2 * np.pi / 12.42

Here is an example 'real' (scalar) dataset:

    elev = 5 * np.cos(m2_freq * t)

Compute the tidal fit:

    tfit_e = tt.t_tide(elev)

All other input is optional. Currently `dt`, `stime`, `lat`, `constitnames`, `output`, `errcalc`, `synth`, `out_style`, and `secular` can be specified. Take a look at the t\_tide docstring for more info on these variables. 

`tfit_e` is an instance of the TTideCon ("TTide Constituents") class. It includes a `t_predic` method that is also availabe as the special `__call__` method. This makes it possible to construct the fitted time-series by simply doing:

    elev_fit = tfit_e(t)

Or extrapolate the fit to other times:

    extrap_fit = tfit_e(np.arange(2000,2500))

And here is an example 'complex' (vector) dataset:

    vel = 0.8 * elev + 1j * 2 * np.sin(m2_freq * t)

    tfit_v = tt.t_tide(vel)

And so on...

Notes
=====

1. The code to handle timeseries longer then 18.6 years has not been converted yet.

2. The code is a little messy and they are a few hacky bits that probably will need to be fixed. The most notable is in noise_realizations. It swaps eig vectors around to match Matlab's output.
Also, the returned diagonal array would sometimes be a negative on the order of 10^-10. Values between (-0.00000000001,0) are forced to 0. 

3. ttide_py was initially converted to python with SMOP. Available at, https://github.com/victorlei/smop.git.
