import os
import numpy as np

__version__ = '0.3'


_base_dir = os.path.join(os.path.dirname(__file__), 'data')
if os.path.exists(os.path.join(_base_dir,'t_constituents_const.npy')) & os.path.exists(os.path.join(_base_dir,'t_constituents_sat.npy')) & os.path.exists(os.path.join(_base_dir,'t_constituents_shallow.npy')):
    _const=np.load(os.path.join(_base_dir,'t_constituents_const.npy'))
    _const=_const[()]
    _sat=np.load(os.path.join(_base_dir,'t_constituents_sat.npy'))
    _sat=_sat[()]
    _shallow=np.load(os.path.join(_base_dir,'t_constituents_shallow.npy'))
    _shallow=_shallow[()]
else:
    print "You do not have t_constituents_*.npy check that package installation is correct."
    const={}
    sat={}
    shallow={}
