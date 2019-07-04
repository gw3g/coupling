from alphas import *
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
    tbl = loadtxt("table_Coupling_{nf"+str(nf)+","+str(l)+"-loop}.dat",usecols=(0,2),unpack=True)
    return tbl

def alpha_k0(tbl,k,T):
    ''' prepare coupling dependence, compatible w/ 1604.07533 '''
    # output
    out = open("coupling_nf0_{"+"k={0:.2f}".format(k)+",t="+str(T)+"}.1.dat",'w')
    out.write("# Columns: k0/T, mu_opt/Lambda, alpha xi={1.,.5,2.}\n")
    out.write("# (Tc = 1.25 Lambda, T="+str(T)+"Tc, nf=0, 5-loop RG)\n")
    # out.write("# (mu_opt=max[(2.pi.T)^2+|K^2|], vary mu=mu_opt.xi)\n")
    out.write("# (mu_opt=sqrt{(xi.pi.T)^2+|K^2|}, vary xi)\n")
    k0=1e-1    # initial k0 value
    Tc = 1.25
    # print("Reading table ...")
    # tbl = ReadTable(0,3)

    while k0<2e1:
        K = ( abs(k0*k0-k*k) )**.5
        # mu=Tc*T*max(K, 2*pi)

        # xi = 1
        mu=Tc*T*(K**2 + (pi)**2)**.5
        t_curr=2*log(mu)
        alpha1 = interpolate_Alpha(tbl,mu)
        # xi = .9
        mu=Tc*T*(K**2 + (.5*pi)**2)**.5
        t_curr=2*log(mu)
        alpha2 = interpolate_Alpha(tbl,mu)

        # xi = 1.1
        mu=Tc*T*(K**2 + (2.*pi)**2)**.5
        t_curr=2*log(mu)
        alpha3 = interpolate_Alpha(tbl,mu)

        out.write(
        "{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  {4:.5e}\n".format(k0,mu,alpha1,alpha2,alpha3) )
        k0+=1e-1

    print("Done.")
    out.close()
    return 0

print("Reading table ...")
tbl = ReadTable(0,5)

T = 1.3
for i in [1,2,3]:
    k = i*2*pi*7/24.
    alpha_k0(tbl,k,T)

T = 1.2
for i in [1,2,3]:
    k = (i**.5)*pi/2.
    alpha_k0(tbl,k,T)

T = 1.1
for i in [1,2,3]:
    k = i*2*pi/3.
    alpha_k0(tbl,k,T)

