ttide_py
========

A direct conversion of T_Tide to Python.



This is a work in progress. It is not done.

It is now mostly functional.
Any help with finishing the conversion is welcome.

All credit for T_Tide goes to Rich Pawlowicz, the original creator of T_Tide. 
It is available at http://www2.ocgy.ubc.ca/~rich/.

This is just a programming language conversion.


Installation
============

This has little to no testing. Use at your own risk. It has only been tested on my personal machine.

Run,

python setup.py install



Usage
=====

from ttide.t_tide import t_tide

[name, freq, tidecon, xout]=t_tide(xin)




All other input is variable.
Currently dt,stime,lat,constitnames,output and errcalc can be specified. To do so use key=value (ex dt=0.5).


dt -              Sampling interval (hours)   default = 1

stime -           Start time as number of days. using mpl.dates rather then datetime because it handles fractional days   default = empty

lat -             Decimal degress (+north)    default = empty

errcalc -         Method for calculation of confidence limits. (cboot,wboot,linear(not finished)) default = 'cboot'

constitnames -    Names of constituents to use. default = empty

output -          Flag to disable output. default = True



Notes
=====

1) The code to handle timeseries longer then 18.6 years has not been converted yet.

2) t_tide doesn't call t_predic at the end yet. The returned xout is incorrect. 

3) t_predic is working. Call it with [xout]=t_predic(time_in,names,freq,tidecon)

4) The code is still pretty messy and they are a few dirty hacks, that probably won't hold up and will need to be fixed. The most notable is in noise_realizations. I swap eig vectors around to match matlabs output.
Also, the returned diagonal matrix would sometimes be a negative on the order of 10^-10. Values between (-0.00000000001,0) are forced to 0. 


