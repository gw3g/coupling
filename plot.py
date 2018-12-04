import matplotlib.pyplot as P
from numpy import loadtxt
from alphas import asymp
from math import log

P.xscale('log')
P.xlabel("t")
P.yscale('log')
P.ylabel("alpha")

a,b,c = loadtxt("nf0_3loop.dat",
                delimiter="  ",unpack=True)
d,e,f = loadtxt("nf0_3loop.dat",
                delimiter="  ",unpack=True)

P.plot(a,b,'r')
P.plot(a,c,'rx')
P.plot(d,e,'purple')
P.plot(d,f,'g.')
P.legend()
P.show()
