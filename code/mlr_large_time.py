# Rigorous large time bound for the MLR property of the SU(2) heat kernel,
# p_t(theta) = sum_n n e^{-(n^2-1)t} U_{n-1}(cos theta), time in units where the Casimir is n^2-1.
# Mixed partial of log p is >= sin(theta)/p * D(t) with D below; MLR holds wherever D(t) > 0.
import numpy as np
from mpmath import mp, mpf, exp, nsum, inf
mp.dps=30
def D(t):
    t=mpf(t); e=lambda n: n*exp(-(n*n-1)*t)
    main=6*e(2)
    s1=nsum(lambda n:(n*n-1)*(n**3-n)/3*e(n),[3,inf])
    A=nsum(lambda n:(n*n-1)*n*e(n),[2,inf]); B=nsum(lambda n:(n**3-n)/3*e(n),[2,inf])
    pmin=1-nsum(lambda n:n*e(n),[2,inf])
    return main-s1-A*B/pmin if pmin>0 else mpf(-1)
lo,hi=mpf('0.1'),mpf(3)
for _ in range(60):
    mid=(lo+hi)/2
    lo,hi=(lo,mid) if D(mid)>0 else (mid,hi)
print("t* =",hi, " D(t*+0.01)=",D(hi+mpf('0.01')))
# check D>0 for all t>=t*: D(t) e^{3t} is increasing; sample
print([float(D(t)*exp(3*t)) for t in [float(hi)+0.01,0.6,1,2,5]])
# numerical region t<t*: mixed partial / sin via image sum
def logf(t,th,K=6):
    k=np.arange(-K,K+1)[:,None]; x=th[None,:]+2*np.pi*k
    return np.log(np.abs((x*np.exp(-x**2/(4*t))).sum(0)))
th=np.linspace(1e-3,np.pi-1e-3,3000); worst=np.inf
for t in np.geomspace(0.01,float(hi)+0.05,600):
    # Casimir n^2-1 time t corresponds to Poisson time t (p ~ e^{t} f / sin)
    h=1e-5*t; d=1e-4
    m=((logf(t+h,th+d)-logf(t-h,th+d))-(logf(t+h,th-d)-logf(t-h,th-d)))/(4*h*d)
    worst=min(worst,(m/np.sin(th)).min())
print("min of mixed partial / sin theta on t<t*:",worst)
