from __future__ import print_function
import os,sys
import numpy as np
import scipy.io.netcdf as netcdf


__version__ = '0.3'


_base_dir = os.path.join(os.path.dirname(__file__), 'data')

if os.path.exists(os.path.join(_base_dir,'t_constituents_const.nc')) & os.path.exists(os.path.join(_base_dir,'t_constituents_sat.nc')) & os.path.exists(os.path.join(_base_dir,'t_constituents_shallow.nc')):
    _const={}
    _sat={}
    _shallow={}
    
    ncid = netcdf.netcdf_file(os.path.join(_base_dir,'t_constituents_const.nc'),'r',mmap=False)
    for key in ncid.variables.keys():
        _const[key]=ncid.variables[key].data
    ncid.close()

    ncid = netcdf.netcdf_file(os.path.join(_base_dir,'t_constituents_sat.nc'),'r',mmap=False)
    for key in ncid.variables.keys():
        _sat[key]=ncid.variables[key].data
    ncid.close()
    
    ncid = netcdf.netcdf_file(os.path.join(_base_dir,'t_constituents_shallow.nc'),'r',mmap=False)
    for key in ncid.variables.keys():
        _shallow[key]=ncid.variables[key].data
    ncid.close()

    _const['name']=np.array([b''.join([s for s in arr]) for arr in _const['name']])
    _const['kmpr']=np.array([b''.join([s for s in arr]) for arr in _const['kmpr']]) 

else:
    print('You do not have t_constituents_*.npy check that package installation is correct.')
    _const={}
    _sat={}
    _shallow={}
