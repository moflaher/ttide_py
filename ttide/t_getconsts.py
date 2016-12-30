from __future__ import division
import numpy as np
import copy
from .t_astron import t_astron
import os.path as path
from scipy.io.netcdf import netcdf_file as nopen

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


def t_getconsts(ctime):
    """
    t_getconsts - Gets constituent data structures holding
                  information for tidal analyses
    Variables are loaded from 't_constituents_*.npy'
    on init and a copy is made now.
    When ctime is specified t_getconsts recomputes the frequencies from
    the rates-of-change of astronomical parameters at the matlab TIME given.

     :Parameters:
        ctime: a datetime, the start time of the data input into t_tide.

    Note:
        Not sure if a copy has to be made here or if they can be used directly.
        For now a copy should be much fast then a load.
    """
    const = copy.deepcopy(_const)
    sat = copy.deepcopy(_sat)
    shallow = copy.deepcopy(_shallow)

    if ctime.size != 0:
        # If no time, just take the "standard" frequencies, otherwise
        # compute them from derivatives of astro parameters. This is
        # probably a real overkill - the diffs are in the
        # 10th decimal place (9th sig fig).
        astro, ader = t_astron(ctime)

        ii = np.isfinite(const['ishallow'])
        const['freq'][~ii] = np.dot(const['doodson'][~ii, :], ader) / 24

        for k in np.flatnonzero(ii):
            ik = ((const['ishallow'][k] - 1 +
                   np.array(range(0, const['nshallow'][k]))).astype(int))
            const['freq'][k] = np.dot(const['freq'][shallow['iname'][ik] - 1],
                                      shallow['coef'][ik])

    return const, sat, shallow
