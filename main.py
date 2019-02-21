from alphas import *
from numpy import loadtxt

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

def ReadTable():
    tbl = loadtxt("table_Coupling_{nf0, 3-loop}.dat",usecols=(0,2),unpack=True)
    return tbl

def alpha_k0():
    T = 1.1
    k = 3*2*pi/3.
    # output
    out = open("coupling2_k3_nf0_GJ.dat",'w')
    out.write("# Columns: k0/T, mu_opt/Lambda, alpha: xi=1,.5,2\n")
    out.write("# (Tc = 1.25 Lambda, T=1.1Tc, nf=0, 3-loop RG)\n")
    #out.write("# (mu_opt=max[(pi.T)^2+|K^2|], vary mu=mu_opt.xi)\n")
    out.write("# (mu_opt=sqrt{(pi.T)^2+|K^2|}, vary mu=mu_opt.xi)\n")
    k0=1e-1    # initial k0 value
    Tc = 1.25
    tbl = ReadTable()

    while k0<1e3:
        K = ( abs(k0*k0-k*k) )**.5
        #mu=Tc*T*max(K, pi)
        mu=Tc*T*(K**2 + pi**2)**.5
        t_curr=2*log(mu)
        print("k0 = ",k0)
        alpha1 = interpolate_Alpha(tbl,mu)
        alpha2 = interpolate_Alpha(tbl,.5*mu)
        alpha3 = interpolate_Alpha(tbl,2.*mu)
        out.write( "{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  {4:.5e}\n".format(k0,mu,alpha1,alpha2,alpha3) )
        k0+=1e-1

    out.close()
    return 0



alpha_k0()

#out = open("error2.dat",'w')
#out.write("# Columns: t, UV, a1(tinf=1e4), a2(tinf=1e2)\n")
#t_s,a_U=solver(.5,1e2,50,1e4)
#t_s,a_L=solver(.5,1e2,50,1e3)
#for i in range(len(t_s)):
    #out.write( "{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}\n".format(t_s[i],a2loop(t_s[i]),a_U[i],a_L[i]) )
#.out.close()
