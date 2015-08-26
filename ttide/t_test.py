from __future__ import division,print_function
import numpy as np
import os
from t_tide import t_tide

#Build an M2 tidal series with 5m elevation
ein=5*np.cos(((2*np.pi/12.42))*np.arange(1001))

#Build an M2 tidal series for current
#Use to 2m/s horizontal currents and .5m/s vertical currents
uin=2*np.cos(((2*np.pi/12.42))*np.arange(1001))
vin=0.5*np.cos(((2*np.pi/12.42))*np.arange(1001))
uvin=uin+1j*vin


###################################################################################################
#Elevation Test Cases
###################################################################################################

#Basic test case
[nameu, freq, tidecon, xout]=t_tide(ein)

print()
print('====================================================================')
print()

#No output test case
[nameu, freq, tidecon, xout]=t_tide(ein,output=False)

print()
print('====================================================================')
print()

#Pandas output case
[nameu, freq, tidecon, xout]=t_tide(ein,out_style='pandas')

print()
print('====================================================================')
print()

#M2 tides only
[nameu, freq, tidecon, xout]=t_tide(ein,constitnames=['M2'])


print()
print('====================================================================')
print()

#5 tidal constituents (all should basically be zero other then M2)
[nameu, freq, tidecon, xout]=t_tide(ein,constitnames=['M2','N2','S2','K1','O1'])


print()
print('====================================================================')
print()

#M2 tides with a starttime
[nameu, freq, tidecon, xout]=t_tide(ein,constitnames=['M2'],stime=768000)


print()
print('====================================================================')
print()

#M2 tides with a starttime and latitude
[nameu, freq, tidecon, xout]=t_tide(ein,constitnames=['M2'],stime=768000,lat=45)




###################################################################################################
#Current Test Cases
###################################################################################################

#Basic test case
[nameu, freq, tidecon, xout]=t_tide(uvin)

print()
print('====================================================================')
print()

#No output test case
[nameu, freq, tidecon, xout]=t_tide(uvin,output=False)

print()
print('====================================================================')
print()

#Pandas output case
[nameu, freq, tidecon, xout]=t_tide(uvin,out_style='pandas')

print()
print('====================================================================')
print()

#M2 tides only
[nameu, freq, tidecon, xout]=t_tide(uvin,constitnames=['M2'])


print()
print('====================================================================')
print()

#5 tidal constituents (all should basically be zero other then M2)
[nameu, freq, tidecon, xout]=t_tide(uvin,constitnames=['M2','N2','S2','K1','O1'])


print()
print('====================================================================')
print()

#M2 tides with a starttime
[nameu, freq, tidecon, xout]=t_tide(uvin,constitnames=['M2'],stime=768000)


print()
print('====================================================================')
print()

#M2 tides with a starttime and latitude
[nameu, freq, tidecon, xout]=t_tide(uvin,constitnames=['M2'],stime=768000,lat=45)



