from __future__ import division,absolute_import
import numpy as np
import os, copy
from t_astron import t_astron
from __init__ import _const,_sat,_shallow

def t_getconsts(ctime):
    """
    t_getconsts - Gets constituent data structures holding information for tidal analyses
    Variables are loaded from 't_constituents_*.npy' on init and a copy is made now.
    When ctime is specified t_getconsts recomputes the frequencies from the rates-of-change of astronomical parameters at the matlab TIME given.

     :Parameters: 
        ctime: a datetime, the start time of the data input into t_tide.     


    Note:
        Not sure if a copy has to be made here or if they can be used directly. For now a copy should be much fast then a load.
    """
    const=copy.deepcopy(_const)
    sat=copy.deepcopy(_sat)
    shallow=copy.deepcopy(_shallow)

       
    if ctime.size!=0:
        # If no time, just take the "standard" frequencies, otherwise compute them from derivatives of astro parameters. 
        # This is probably a real overkill - the diffs are in the 10th decimal place (9th sig fig).
        astro, ader = t_astron(ctime)
       
        ii=np.isfinite(const['ishallow'])
        const['freq'][~ii] = np.dot(const['doodson'][~ii, :], ader)/24

        for k in np.flatnonzero(ii):
            ik = ((const['ishallow'][k]-1 + np.array( range(0,const['nshallow'][k]) ) ).astype(int))    
            const['freq'][k] = np.dot(const['freq'][shallow['iname'][ik]-1], shallow['coef'][ik])

    return const, sat, shallow
