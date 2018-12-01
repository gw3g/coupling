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

r = ode(_F,_J).set_integrator('vode', method='bdf', with_jacobian=True)
l=2
tinf=1000
ainf=asymp(tinf,l)
r.set_initial_value(ainf, tinf).set_f_params(l).set_jac_params(l)
tf = 1.4
tvs = []
res = []

while r.successful() and r.t > tf:
    r.integrate(r.t*.98)
    if r.t<10:
        tvs.append(r.t)
        res.append(r.y[0])
    print(r.t, r.y[0])

t_vals=[exp(t*.01) for t in range(1,200)]
t1 = [asymp(exp(t*.01),1) for t in range(1,200)]
t1 = [asymp(exp(t*.01),1) for t in range(1,200)]
t2 = [asymp(exp(t*.01),2) for t in range(1,200)]
t3 = [asymp(exp(t*.01),3) for t in range(1,200)] 

xscale('log')
# yscale('log')
plot(t_vals,t3)
plot(tvs,res,'ro')
show()
