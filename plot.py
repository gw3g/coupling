import matplotlib.pyplot as P
from numpy import loadtxt
from alphas import asymp
from math import log
P.rc('text', usetex=True)

def plot_nf0():
    P.title(r"Running coupling, $n_f=0$")
    P.xscale('log')
    P.xlabel(r"$t=2\log(\mu/\Lambda)$")
    P.xlim(1,10)
    P.ylabel(r"$\alpha$")
    P.ylim(0,1)

    a1,b1,c1 = loadtxt("nf0_1loop.dat",
                delimiter="  ",unpack=True)
    a2,b2,c2 = loadtxt("nf0_2loop.dat",
                delimiter="  ",unpack=True)
    a3,b3,c3 = loadtxt("nf0_3loop.dat",
                delimiter="  ",unpack=True)
    a4,b4,c4 = loadtxt("nf0_4loop.dat",
                delimiter="  ",unpack=True)

    P.plot(a1,b1,'c:',label="1-loop, UV")
    P.plot(a1,c1,'c-',label="        RG")
    P.plot(a2,b2,'b:',label="2-loop, UV")
    P.plot(a2,c2,'b-',label="        RG")
    P.plot(a3,b3,'g:',label="3-loop, UV")
    P.plot(a3,c3,'g-',label="        RG")
    P.plot(a4,b4,'y:',label="4-loop, UV")
    P.plot(a4,c4,'y-',label="        RG")
    P.legend()
    P.savefig("nf=0.pdf")

def plot_nf3():
    P.title(r"Running coupling, $n_f=3$")
    P.xscale('log')
    P.xlabel(r"$\mu /$GeV")
    P.xlim(.1,200)
    P.ylabel(r"$\alpha (\mu)$")
    P.ylim(0,1)

    a1,b1,c1 = loadtxt("nf3_1loop.dat",
                delimiter="  ",unpack=True)
    a2,b2,c2 = loadtxt("nf3_2loop.dat",
                delimiter="  ",unpack=True)
    a3,b3,c3 = loadtxt("nf3_3loop.dat",
                delimiter="  ",unpack=True)
    a4,b4,c4 = loadtxt("nf3_4loop.dat",
                delimiter="  ",unpack=True)
    a5,b5,c5 = loadtxt("nf3_5loop.dat",
                delimiter="  ",unpack=True)
    a,b,c = loadtxt("nf3_analytic.dat",
                delimiter="  ",unpack=True)

    P.axvline(x=.341,color='.7',linestyle=':')
    P.text(.45,.05,r"$\Lambda_{\overline{\rm MS}}$")

    P.plot(a1,b1,'g-',label=r"1-loop ($\overline{\rm MS}$)")
    #P.plot(a1,c1,'c-',label="        RG")
    #P.plot(a2,b2,'b:',label="2-loop, UV")
    P.plot(a,c,'b-',label=r"analytic extension, timelike")
    P.plot(a,b,'c-',label=r"analytic extension, spacelike")
    P.plot(a2,c2,'r:',label=r"2-loop")
    P.plot(a3,c3,'r--',label=r"3-loop")
    P.plot(a4,c4,'r-.',label=r"4-loop")
    P.plot(a5,c5,'r-',label=r"5-loop")
    P.legend()
    P.savefig("nf=3.pdf")


def plot_check():
    P.title(r"1-loop couplings, $n_f=0$")
    P.xscale('log')
    P.xlabel(r"$t=2\log(\mu/\Lambda)$")
    P.yscale('log')
    P.ylabel(r"$\alpha$")

    a1,b1,c1,d1 = loadtxt("error2.dat",
                delimiter="  ",unpack=True)
    a2,b2,c2   = loadtxt("analytic.dat",
                delimiter="  ",unpack=True)

    P.plot(a1,c1,'m-',label="$t_\infty = 10^4$")
    P.plot(a1,d1,'y-',label="$t_\infty = 10^2$")
    P.fill_between(a1,c1,d1,color='y',alpha=.5)
    P.plot(a1,b1,'r:',label="2-loop, full")
    # P.plot(a2,b2,'g-',label="timelike $t>0$")
    # P.plot(a2,c2,'g:',label="spacelike $t<0$")
    P.legend()
    P.savefig("error2.pdf")

plot_nf3()
