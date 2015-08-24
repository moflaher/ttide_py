from __future__ import division
import numpy as np
import os
from t_astron import t_astron

def t_getconsts(ctime):
    """
    t_getconsts - Gets constituent data structures holding information for tidal analyses
    Variables are loaded from 't_constituents_*.npy'
    When ctime is specified t_getconsts recomputes the frequencies from the rates-of-change of astronomical parameters at the matlab TIME given.

     :Parameters: 
        ctime: a datetime, the start time of the data input into t_tide.     

    """

    base_dir = os.path.join(os.path.dirname(__file__), 'data')
    if os.path.exists(os.path.join(base_dir,'t_constituents_const.npy')) & os.path.exists(os.path.join(base_dir,'t_constituents_sat.npy')) & os.path.exists(os.path.join(base_dir,'t_constituents_shallow.npy')):
        const=np.load(os.path.join(base_dir,'t_constituents_const.npy'))
        const=const[()]
        sat=np.load(os.path.join(base_dir,'t_constituents_sat.npy'))
        sat=sat[()]
        shallow=np.load(os.path.join(base_dir,'t_constituents_shallow.npy'))
        shallow=shallow[()]
    else:
        print("You do not have t_constituents_*.npy check that package installation is correct.")
        const=[]
        sat=[]
        shallow=[]
       
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
