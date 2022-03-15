#from alphas import *
from math import pi, log, tanh
from numpy import loadtxt

#-----------------------------------------------------------------------------#
def interpolate_Alpha(tbl,mu):
    ''' tbl is a 2D array, mu is argument of coupling '''
    mu0,muF=tbl[0][0],tbl[0][-1]

    # no misbehaviour please ...
    if (mu>muF):
        print( "Out of bounds: "+str(mu)+" = mu > MAX = "+str(muF) )
        return 0
    if (mu<mu0):
        print( "Out of bounds: "+str(mu)+" = mu < MIN = "+str(mu0) )
        return 0
    dmu=tbl[0][1]-tbl[0][0]
    mu0=tbl[0][0]
    i=int((mu-mu0)/dmu)
    t =log(mu)
    tA=log(i*dmu+mu0)
    tB=log((i+1)*dmu+mu0)
    x=(t-tA)/(tB-tA) # linear param. in [0,1]
    A1,A2=tbl[1][i],tbl[1][i+1]
    return (1.-x)*A1+(x)*A2

def ReadTable(nf,l):
    tbl = loadtxt("../table_Coupling_{nf"+str(nf)+","+str(l)+"-loop}.dat",usecols=(0,2),unpack=True)
    return tbl

def alpha_mu(nf,l):
    Tb = ReadTable(nf,l)
    ''' prepare coupling, mu = [2,10] '''
    out = open("coupling_nf"+str(nf)+"_"+str(l)+"loop.dat",'w')
    out.write("# Columns: mu/Lambda, alpha\n")
    mu = 2

    while mu<10.:
        alpha1 = interpolate_Alpha(Tb,mu)
        out.write(
        "{0:.5e}  {1:.5e}\n".format(mu,alpha1) )
        mu+=1e-1

    print("Done.")
    out.close()
    return 0

def mu_opt(k0,k,T,lam,r=1,xi=1):
    ''' argument list: k0/T, k/T, T/Tc, Lambda/Tc '''
    K = ( T*T*abs(k0*k0-k*k) )**.5

    qT = xi*pi*T
    ## schemes R.1, R.2, R.3
    if (r==1): mu = (K**2 + qT**2)**.5
    if (r==2): mu = max(K,qT)
    if (r==3): mu = qT/tanh(qT/K)
    #else: return -1

    #print(mu/lam)
    return mu/lam


def alpha_k0(tbl,nf,loops,k,T,tag):
    ''' prepare coupling dependence, compatible w/ 1604.07533 '''
    r = 1 # flag for running scheme

    # output
    out = open("R."+str(r)+"/coupling_nf"+str(nf)+"_{"+"k={0:.2f}".format(k)+",t="+str(T)+"}."+tag+".dat",'w')
    out.write("# Columns: k0/T, mu_opt/Lambda, alpha xi={1,.5,2}\n")

    if (nf==0):
        lam = 1/1.24
        xi  = 1
    if (nf==2):
        lam = 1/0.56
        xi  = 2

    out.write("# (Tc = "+str(1/lam)+" Lambda, T="+str(T)+"Tc, nf="+str(nf)+", "+str(loops)+"-loop RG)\n")

    blurb = ["mu_opt=sqrt{("+str(xi)+".pi.T)^2+|K^2|}","mu_opt=max[|K|,pi.T]","pi.T/tanh(pi.T/|K|)"]
    out.write("# ("+blurb[r-1]+", vary mu=mu_opt.xi)\n")

    k0=1e-1    # initial k0 value

    while k0<=1e3:
        mu=mu_opt(k0,k,T,lam,r,xi)

        # xi = 1
        alpha1 = interpolate_Alpha(tbl,mu)
        # xi = .9
        alpha2 = interpolate_Alpha(tbl,.5*mu)
        # xi = 1.1
        alpha3 = interpolate_Alpha(tbl,2.*mu)

        out.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  {4:.5e}\n".format(k0,mu,alpha1,alpha2,alpha3) )
        k0+=1e-1

    print("Done.")
    out.close()
    return 0

loops = 5
nf    = 0
tag   = "R1(5)"
print("Reading table ...")
tbl = ReadTable(nf,loops)

if (nf==0):
    T = 1.3
    for i in [1,2,3]:
        k = i*2*pi*7/24.
        alpha_k0(tbl,nf,loops,k,T,tag)
    T = 1.1
    for i in [1,2,3]:
        k = i*2*pi/3.
        alpha_k0(tbl,nf,loops,k,T,tag)

if (nf==2):
    T = 1.2
    alpha_k0(tbl,nf,loops,pi,T,tag)
    alpha_k0(tbl,nf,loops,1.5*pi,T,tag)
    alpha_k0(tbl,nf,loops,(14.**.5)*pi/2.,T,tag)
    for i in [1,2,3]:
        k = (i**.5)*pi/2.
        alpha_k0(tbl,nf,loops,k,T,tag)

