import matplotlib.pyplot as P
from numpy import loadtxt
from alphas import asymp
from math import log

P.xscale('log')
P.xlabel("k_0/T")
#P.yscale('log')
P.ylabel("alpha/pi")

a,b,c = loadtxt("coupling_nf0_11.dat",
                delimiter="  ",unpack=True)

P.plot(a,c,'ro')
P.legend()
P.show()
