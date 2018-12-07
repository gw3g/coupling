import matplotlib.pyplot as P
from numpy import loadtxt
from alphas import asymp
from math import log

P.rc('text', usetex=True)

P.title(r"Running coupling, $n_f=0$")
P.xscale('log')
P.xlabel(r"$t=2\log(\mu/\Lambda)$")
# P.yscale('log')
P.ylabel(r"$\alpha$")
P.ylim(0,1)

a,b,c = loadtxt("nf0_1loop.dat",
                delimiter="  ",unpack=True)
d,e,f = loadtxt("nf0_3loop.dat",
                delimiter="  ",unpack=True)

P.plot(a,b,'r:',label="1-loop, UV")
P.plot(a,c,'r.',label="        RG")
P.plot(d,e,'b:',label="3-    , UV")
P.plot(d,f,'b+',label="        RG")
P.legend()
P.savefig("1.pdf")
