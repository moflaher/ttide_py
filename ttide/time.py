"""
This modules is a reimplemntation of the MatPlotLib 'mpltime' format
(ordinal + fractional time).

This is an independent implementation of the mpltime format so that a
package that uses mpltime does not require the entire MPL library.

See this for more information:
http://matplotlib.org/api/dates_api.html

"""
from __future__ import division
import numpy as np
from datetime import datetime, timedelta


def num2date(mpltime):
    if np.ndarray in mpltime.__class__.__mro__:
        out = np.empty(len(mpltime), dtype='O')
        for idx, val in enumerate(mpltime.flat):
            out[idx] = num2date(val)
        out.shape = mpltime.shape
        return out
    return datetime.fromordinal(int(mpltime)) + timedelta(days=mpltime % 1)


def date2num(dt):
    if isinstance(dt, np.ndarray):
        if dt.dtype.name.startswith('datetime64'):
            dt = dt.astype('O')
        out = np.empty(len(dt), dtype=np.float64)
        for idx, val in enumerate(dt.flat):
            out[idx] = date2num(val)
        out.shape = dt.shape
        return out
    return (dt.toordinal() +
            (((dt.microsecond / 1e6 +
               dt.second) / 60 +
              dt.minute) / 60 +
             dt.hour) / 24)


# Not sure this is useful here...
# def mpltime2matlab_datenum(time):
#     return time.view(np.ndarray) + 366
