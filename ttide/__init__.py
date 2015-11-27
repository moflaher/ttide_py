from __future__ import print_function
import os.path as path
import sys
import numpy as np
from scipy.io.netcdf import netcdf_file as nopen


__version__ = '0.3'


_base_dir = path.join(path.dirname(__file__), 'data')
has_const = path.exists(path.join(_base_dir, 't_constituents_const.nc'))
has_sat = path.exists(path.join(_base_dir, 't_constituents_sat.nc'))
has_shallow = path.exists(path.join(_base_dir, 't_constituents_shallow.nc'))

if (has_const and has_sat and has_shallow):
    _const = {}
    _sat = {}
    _shallow = {}

    ncid = nopen(path.join(_base_dir,
                           't_constituents_const.nc'), 'r', mmap=False)
    for key in ncid.variables.keys():
        _const[key] = ncid.variables[key].data
    ncid.close()

    ncid = nopen(path.join(_base_dir,
                           't_constituents_sat.nc'), 'r', mmap=False)
    for key in ncid.variables.keys():
        _sat[key] = ncid.variables[key].data
    ncid.close()

    ncid = nopen(path.join(_base_dir,
                           't_constituents_shallow.nc'), 'r', mmap=False)
    for key in ncid.variables.keys():
        _shallow[key] = ncid.variables[key].data
    ncid.close()

    # Correct issues with name strings
    _const['name'] = np.array([b''.join([s for s in arr])
                               for arr in _const['name']])

    _const['kmpr'] = np.array([b''.join([s for s in arr])
                               for arr in _const['kmpr']])

else:
    print('You do not have t_constituents_*.npy ' +
          'check that package installation is correct.')
    _const = {}
    _sat = {}
    _shallow = {}
