# How to put lipstick on a pig
from scipy.integrate import ode
from scipy.special import zeta
from math import pi, log, exp
from matplotlib.pyplot import plot, ylabel, show, xscale, yscale

nf = 0 # flavours
z3 = zeta(3)
b = [(33-2*nf)/(12*pi),
     (153-19.*nf)/(24*pi**2),
     (2857-(5033/9)*nf+(325/27)*nf**2)/(128*pi**3),
     ((149753/6+z3*3564)-(1078361/162+z3*6508/27)*nf+
      (50065/162+z3*6472/81)*nf**2+(1093/728)*nf**3)/(256*pi**4)]

def asymp(t, L):
    """ L-loop UV alpha, t=2*log(mu/lam) """
    res  = 1
    if L>0: res += -b[1]*log(t)/(b[0]**2*t)
    if L>1: res += +(b[1]**2*(log(2)**2-log(t)-1)+b[0]*b[2])/(b[0]**4*t**2)
    if L>2: res += -(b[1]**3-b[0]**2*b[3]-4*b[1]**3*log(t)+6*b[0]*b[1]*b[2]*log(t)
                     -5*b[1]**3*log(t)**2+2*b[1]**3*log(t)**3)/(2*b[0]**6*t**3)
    res /= b[0]*t
    return res

def F(t, a, L): # da/dt = ...
    if L==0:
        return -b[0]*a**2
    else:
        return -b[L]*a**(L+2)+F(t, a, L-1)

def J(t, a, L): # jacobian
    if L==0:
        return -2*b[0]*a
    else:
        return -(L+2)*b[L]*a**(L+1)+J(t, a, L-1)

_F = lambda t, a, L :  [F(t,a,L)]
_J = lambda t, a, L : [[J(t,a,L)]]

def solver(tA,tB,n):
    """ tA<tB, log scaling with n points """
    crank=ode(_F,_J).set_integrator('vode', method='bdf', with_jacobian=True)
    tinf=1000   # UV boundary conditions
    l=2         # loop-order
    r=(tB/tinf)**(1/n)
    t_val,a_val = [],[]
    crank.set_initial_value(asymp(tinf,l), tinf).set_f_params(l).set_jac_params(l)
    while crank.successful() and crank.t > tB:
        crank.integrate(crank.t*r)
    r=(tA/tB)**(1/n)
    while crank.successful() and crank.t > tA:
        crank.integrate(crank.t*r)
        t_val.append(crank.t)
        a_val.append(crank.y[0])
    return t_val, a_val

def alpha_ML():
    crank=ode(_F,_J).set_integrator('vode', method='bdf', with_jacobian=True)
    out = open("out.dat",'w')
    st_l= []
    tinf=1e2  # UV boundary conditions
    l=3         # loop-order
    k0=1e3
    crank.set_initial_value(asymp(tinf,l), tinf).set_f_params(l).set_jac_params(l)
    mu=max(k0*1.25*1.1, pi*1.1*1.25)
    t_curr=2*log(mu)
    while crank.successful() and crank.t > t_curr:
        crank.integrate(crank.t*.99)
    while crank.successful() and k0>0:
        mu=max(k0*1.25*1.1, pi*1.1*1.25)
        t_curr=2*log(mu)
        crank.integrate(t_curr)
        # out.write("{0:.5e}  {1:.5e}  {2:.5e}\n".format(k0,mu,crank.y[0]/pi))
        st_l.insert(0,"{0:.5e}  {1:.5e}  {2:.5e}\n".format(k0,mu,crank.y[0]/pi))
        k0-=.1
    for st in st_l: out.write(st)
    out.close()

alpha_ML()

# t_UV=[exp(t*.01) for t in range(1,500)]
# a_UV=[asymp(exp(t*.01),2) for t in range(1,500)] 

# t_s,a_s = solver(1.,100,50)

# out = open("out.dat",'w')

# for i in range(1,len(t_s)):
    # out.write("{0:.5e}  {1:.5e}\n".format(t_s[i],a_s[i]))

# out.close()


# xscale('log')
# yscale('log')
# plot(t_UV,a_UV)
# plot(t_s,a_s,'ro')
# show()
