# How to put lipstick on a pig
from math import pi, log, exp, atan
from scipy.special import zeta, lambertw

b  = [0]*5

def qcd(nf):
    z3 = zeta(3)
    z4 = zeta(4)
    z5 = zeta(5)

    # arxiv:1606.08659
    global b
    b[0]= (33-2*nf)/(12*pi)
    b[1]= (153-19.*nf)/(24*pi**2)
    b[2]= (2857-(5033/9)*nf+(325/27)*nf**2)/(128*pi**3)
    b[3]= ( (149753/6+z3*3564)-(1078361/162+z3*6508/27)*nf+
            (50065/162+z3*6472/81)*nf**2+(1093/728)*nf**3 )/(256*pi**4)
    b[4]= ( (8157455/16+z3*621885/2-z4*88209/2-z5*288090)+
            (-336460813/1944-z3*4811164/81+z4*33935/6+z5*1358995/27)*nf+
            (25960913/1944+z3*698531/81-z4*10526/9-z5*381760/81)*nf**2+
            (-630559/5832-z3*48722/243+z4*1618/27+z5*460/9)*nf**3+
            (1205/2916-z3*152/81)*nf**4 )/(1024*pi**5)

qcd(0)

#-----------------------------------------------------------------------------#
# some functions

def A_analy(t):
    """ analytic running, [ space, time ]-like """
    return [ (1/t+1/(1-exp(t)))/b[0], (.5-atan(-t/pi)/pi)/b[0] ]

def A_2loop(t):
    """ implicit solution for 2-loop coupling """
    r = b[0]**2/b[1]
    z = -r*exp(-1-t*r)      # arg of product log
    return -b[0]/(b[1]*(1+lambertw(z,-1).real) )

def A_asymp(t, L):
    """ L-loop UV alpha, t=2*log(mu/lam) """
    res  = 1
    if L>1: res += -b[1]*log(t)/(b[0]**2*t)
    if L>2: res += +(b[1]**2*(log(2)**2-log(t)-1)+b[0]*b[2])/(b[0]**4*t**2)
    if L>3: res += -(b[1]**3-b[0]**2*b[3]-4*b[1]**3*log(t)+6*b[0]*b[1]*b[2]*log(t)
                     -5*b[1]**3*log(t)**2+2*b[1]**3*log(t)**3)/(2*b[0]**6*t**3)
    if L>4: res += +(7*b[1]**4-18*b[0]*b[1]**2*b[2]+10*b[0]**2*b[2]**2-b[0]**2*b[1]*b[3]
                     +2*b[0]**3*b[4]+24*b[1]**4*log(t)-18*b[0]*b[1]**2*b[2]*log(t)
                     -12*b[0]**2*b[1]*b[3]*log(t)-9*b[1]**4*log(t)**2
                     +36*b[0]*b[1]**2*b[2]*log(t)**2-26*b[1]**4*log(t)**3
                     +6*b[1]**4*log(t)**4)/(6*b[0]**8*t**4)
    if L>5: print("six loop (or higher)? C'mon ...")
    res /= b[0]*t
    return res

#-----------------------------------------------------------------------------#
# numerical RG-evolution
from scipy.integrate import ode

def F(t, a, L): # da/dt = ...
    if L==1:
        return -b[0]*a**2
    else:
        return -b[L-1]*a**(L+1)+F(t, a, L-1)

def J(t, a, L): # jacobian
    if L==1:
        return -2*b[0]*a
    else:
        return -(L+1)*b[L-1]*a**L+J(t, a, L-1)

_F = lambda t, a, L :  [F(t,a,L)]
_J = lambda t, a, L : [[J(t,a,L)]]

def solver(t_min,l,t_inf=1e3):
    """ l loop-order, RG from t_inf """
    if (t_min>t_inf/10): t_inf=1e3*t
    crank=ode(_F,_J).set_integrator('vode', method='bdf', with_jacobian=True)
    #crank=ode(_F,_J).set_integrator('dopri5', rtol=1e-3); # it seems that Dormand-Prince RK works 
                                                          # for higher accuracy // 7.12.18
    r=(t_min/t_inf)**(1e-5) # here 10^5 iterations
    crank.set_initial_value(A_asymp(t_inf,l), t_inf).set_f_params(l).set_jac_params(l)

    while crank.successful() and crank.t > t_min:
        crank.integrate(crank.t*r)
    if not crank.successful():
        print("\n  ... Ooops, the integrator failed!\n")
        return (-1)
    return crank.y[0]

#-----------------------------------------------------------------------------#

def alpha_mu(nf,l):
    qcd(nf)
    crank=ode(_F,_J).set_integrator('vode', method='bdf', with_jacobian=True)
    out = open("mu_nf"+str(nf)+"_"+str(l)+"loop.dat",'w')
    st_l= []
    t_inf=1e3  # UV boundary conditions
    mu=200
    L_msbar=.341
    t_curr=2*log(mu/L_msbar)

    crank.set_initial_value(A_asymp(t_inf,l), t_inf).set_f_params(l).set_jac_params(l)
    while crank.successful() and crank.t > t_curr:
        crank.integrate(crank.t*.99)
    while crank.successful() and (crank.y[0]<1.):
        mu*=.99
        t_curr=2*log(mu/L_msbar)
        crank.integrate(t_curr)
        st_l.insert(0,"{0:.5e}  {1:.5e}  {2:.5e}\n".format(mu,A_asymp(t_curr,l),crank.y[0]))

    # output
    out.write("# Columns: mu/GeV, UV asymptotics, alpha\n")
    out.write("# ( Lambda=341[MeV], nf="+str(nf)+", "+str(l)+"-loop )\n")
    for st in st_l: out.write(st)
    out.close()

#-----------------------------------------------------------------------------#

def alpha_T(nf,l):
    ''' temperature dependence '''
    qcd(nf)
    #crank=ode(_F,_J).set_integrator('dopri5', rtol=1e-4);
    #crank=ode(_F,_J).set_integrator('vode', method='bdf', with_jacobian=True)
    out = open("coupling_k3_nf"+str(nf)+"_"+str(l)+"loop.dat",'w')
    st_l= []
    k = 3.*2*pi/3.

    t_inf=1e4  # UV boundary conditions
    k0=10    # initial k0 value
    Tc = 1.25
    Tt = 1.1
    K = (abs(k0*k0-k*k))**.5
    mu=max(K*Tc*Tt, pi*Tc*Tt)
    t_curr=2*log(mu)

    #crank.set_initial_value(A_asymp(t_inf,l), t_inf).set_f_params(l).set_jac_params(l)
    #while crank.successful() and crank.t > t_curr:
    #    crank.integrate(crank.t*.99)
    while k0>1e-4:
        K = (abs(k0*k0-k*k))**.5
        mu=max(K*Tc*Tt, pi*Tc*Tt)
        t_curr=2*log(mu)
        print("k0 = ",k0)
        #crank.integrate(t_curr)
        res = solver(t_curr,3,1e3)
        # out.write("{0:.5e}  {1:.5e}  {2:.5e}\n".format(k0,mu,crank.y[0]/pi))
        st_l.insert(0,"{0:.5e}  {1:.5e}  {2:.5e}\n".format(k0,mu,res))
        k0-=1e-2
    # output
    out.write("# Columns: k0/T, mu/Lambda, alpha/pi\n")
    out.write("# (Tc = 1.25 Lambda, T=1.1Tc, nf=0, 3-loop )\n")
    for st in st_l: out.write(st)
    out.close()

