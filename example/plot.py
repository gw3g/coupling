import matplotlib.pyplot as P
from numpy import loadtxt
from math import log, pi
P.rc('text', usetex=True)

def plot_mu():
    P.title(r"Running coupling, $n_f=2$")
    P.xlabel(r"$\mu /\Lambda$")
    P.xlim(2.,10)
    P.ylabel(r"$\alpha (\mu)$")
    P.ylim(0,1)

    a1,b1 = loadtxt("coupling_nf2_5loop.dat",
                delimiter="  ",unpack=True)

    #P.axvline(x=.341,color='.7',linestyle=':')
    #P.text(.45,.05,r"$\Lambda_{\overline{\rm MS}}$")

    P.plot(a1,b1,'g-',label=r"5-loop (RG)")
    #P.plot(a1,c1,'c-',label="        RG")
    #P.plot(a2,b2,'b:',label="2-loop, UV")
    P.legend()
    P.show()
    # P.savefig("nf=3.png")


def plot_k0():
    T = 1.2
    #k = (3**.5)*pi/2.
    k = 1.5*pi
    tag = "1"

    P.title(r"$k=$"+str(k)+r" : 5-loop coupling, $n_f=2$")
    P.xlabel(r"$k_0/T$")
    P.ylabel(r"$\alpha\big(\sqrt{|K^2|+(\xi 2\pi T)^2}\big)$")

    fin = "coupling_nf2_{"+"k={0:.2f}".format(k)+",t="+str(T)+"}."+tag+".dat"


    a1,b1,c1,d1,e1 = loadtxt(fin,
                delimiter="  ",unpack=True)

    P.plot(a1,c1,'r:',label=r"$\xi = .5 $")
    P.plot(a1,e1,'r-.',label=r"$\xi = 2 $")
    P.fill_between(a1,c1,e1,color='y',alpha=.5)
    P.plot(a1,d1,'r-')
    # P.plot(a2,b2,'g-',label="timelike $t>0$")
    # P.plot(a2,c2,'g:',label="spacelike $t<0$")
    P.legend()
    P.show()
    #P.savefig("a.pdf")

plot_k0()
