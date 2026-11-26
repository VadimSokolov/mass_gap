# Numerical companion to the initial layer of the contraction theorem (all t <= 0.1).
#  1. The two image model z_2 on eps in [t/4, pi/2]: its minimum, and the three lower bounds used in the proof.
#  2. The error terms: bounds 3.2e6 t^-5 exp(-2 pi^2/t) (near the identity) and 1.6e7 t^-6 exp(-51.8/t) (the middle), at t = 0.1.
#  3. The true z from thirteen images against the models, in 60 digit arithmetic, at t = 0.1.
import mpmath as mp, numpy as np
mp.mp.dps = 60
def z2(t, e):
    u = 2*mp.pi*e/t; tau = mp.tanh(u); D = mp.pi*tau-e; return mp.cot(e)+1/D-(u/mp.sinh(u))**2*(mp.pi*tau/D)**2*(1-(2*e*e+t)/(2*mp.pi**2))/e
def phi(t, th): return mp.log(sum((th+2*mp.pi*j)*mp.exp(-(th+2*mp.pi*j)**2/t) for j in range(-6, 7))/mp.sin(th))
z = lambda t, th: mp.diff(lambda s: s*mp.diff(lambda a: phi(s, a), th), t)
for t in (0.1, 0.03, 0.002):
    t = mp.mpf(t); es = [t/4*(mp.pi/2/(t/4))**(mp.mpf(i)/400) for i in range(401)]; v = [z2(t, e) for e in es]; i = int(np.argmin([float(x) for x in v]))
    lo = lambda e: (mp.mpf("0.517")/e-e/2 if e <= 0.05 else (mp.mpf("0.92")/e-mp.mpf("0.05") if e <= 0.1 else mp.mpf("0.29")))
    print("t = %-6s min of z_2 on [t/4, pi/2] = %s at eps = %s;  lower bounds of the proof hold on the grid: %s" % (mp.nstr(t, 3), mp.nstr(v[i], 6), mp.nstr(es[i], 4), all(v[k] >= lo(es[k]) for k in range(401))))
t = mp.mpf("0.1"); print("error bounds at t = 0.1: near the identity %s, the middle %s" % (mp.nstr(mp.mpf("3.2e6")*t**-5*mp.exp(-2*mp.pi**2/t), 3), mp.nstr(mp.mpf("1.6e7")*t**-6*mp.exp(-mp.mpf("51.8")/t), 3)))
for th in (0.01, 0.5, 1.5): print("t = 0.1, theta = %.2f: z = %s, 1/theta - cot theta = %s" % (th, mp.nstr(z(t, mp.mpf(th)), 20), mp.nstr(1/mp.mpf(th)-mp.cot(th), 20)))
for e in (t/4, 0.05, 0.1, 1.0, mp.pi/2): print("t = 0.1, eps = %s: z = %s, z_2 = %s, relative difference %s" % (mp.nstr(e, 4), mp.nstr(z(t, mp.pi-e), 12), mp.nstr(z2(t, e), 12), mp.nstr(abs(z(t, mp.pi-e)/z2(t, e)-1), 3)))
