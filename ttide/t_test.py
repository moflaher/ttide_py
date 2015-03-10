from __future__ import division
import numpy as np
import os
from t_tide import t_tide

#Build an M2 tidal series with 5m elevation
din=5*np.cos(((2*np.pi/12.42))*np.arange(1000))


#Basic test case
[nameu, freq, tidecon, xout]=t_tide(din)

print
print '===================================================================='
print

#M2 tides only
[nameu, freq, tidecon, xout]=t_tide(din,constitnames=['M2'])


print
print '===================================================================='
print

#5 tidal constituents (all should basically be zero other then M2)
[nameu, freq, tidecon, xout]=t_tide(din,constitnames=['M2','N2','S2','K1','O1'])


print
print '===================================================================='
print

#M2 tides with a starttime
[nameu, freq, tidecon, xout]=t_tide(din,constitnames=['M2'],stime=768000)


print
print '===================================================================='
print

#M2 tides with a starttime and latitude
[nameu, freq, tidecon, xout]=t_tide(din,constitnames=['M2'],stime=768000,lat=45)




