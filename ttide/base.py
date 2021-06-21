import sys
import io
from . import t_utils as tu
from .t_predic import t_predic

# file object depending on Python version
if sys.version_info[0] < 3:
    FILE_OBJ = file  # Python 2
else:
    FILE_OBJ = io.IOBase  # Python 3


class TTideCon(dict):
    """The ttide consituents class.

    This class is based on a dictionary, and has key-value pairs of
    the relavent data from the tidal fit, and for tidal prediction
    (extrapolation).  These include:

    dt : The sampling interval of the fit data.

    nameu : The names of the tidal constituents.

    fu : The frequencies of the tidal constituents.

    tidecon : The tidal constituent amplitudes.

    snr : The signal to noise ratio of the constituent fits.

    """

    def t_predic(self, time):
        return t_predic(time,
                        names=self['nameu'], freq=self['fu'],
                        tidecon=self['tidecon'], lat=self['lat'],
                        ltype=self['ltype'], synth=self['synth'])

    __call__ = t_predic

    def pandas_style(self, to_file=None, to_file_df=None):
        if to_file_df is None:
            outstr = tu.pandas_style(self)            
        else:
            outstr, df = tu.pandas_style(self, True)
            df.to_csv(to_file_df)
            
        if to_file is None:
            return outstr
        elif isinstance(to_file, FILE_OBJ):
            to_file.write(outstr)
        else:
            with open(to_file, 'w') as fl:
                fl.write(outstr)

    def classic_style(self, to_file=None):
        outstr = tu.classic_style(self)
        if to_file is None:
            return outstr
        elif isinstance(to_file, FILE_OBJ):
            to_file.write(outstr)
        else:
            with open(to_file, 'w') as fl:
                fl.write(outstr)
