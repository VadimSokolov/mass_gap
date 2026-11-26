# Exact weak coupling expansion of the projected ladder for SU(2), general bond moving power k and decimation power D.
# Small time heat kernel (leading image; the others are exponentially small): K_t(theta) ~ (theta/sin theta) exp(-theta^2/t).
# Bond moving: density ~ K_t^k sin^2(theta). New u_2 = E cos(theta). Decimation: u_2 -> u_2^D. Coordinate beta = 2/t, u_2 = exp(-3t/4).
import sympy as sp
tau, r, k, D, be = sp.symbols("tau r k D beta", positive=True)
th = sp.sqrt(tau)*r; m = k-2; ORD = 5
lw = m*(th**2/6+th**4/180+th**6/2835+th**8/37800+th**10/467775)               # log (theta/sin theta)^m
w = sp.series(sp.exp(lw), tau, 0, ORD+1).removeO()
cs = sp.series(sp.cos(th), tau, 0, ORD+1).removeO()
mom = lambda e: sp.integrate(sp.expand(e)*r**2*sp.exp(-r**2), (r, 0, sp.oo))   # theta^2 exp(-theta^2/tau), tau = t/k
num = sp.series(mom(sp.expand(cs*w)), tau, 0, ORD).removeO(); den = sp.series(mom(w), tau, 0, ORD).removeO()
lnu = sp.series(sp.log(num/den), tau, 0, ORD).removeO()
t = sp.symbols("t", positive=True)
tp = sp.series((-sp.Rational(4, 3)*D*lnu).subs(tau, t/k), t, 0, ORD).removeO()          # new time
bp = sp.series((2/tp).subs(t, 2/be), be, sp.oo, 4).removeO()
dec = sp.simplify(sp.expand(be-bp)); print("t' =", sp.factor(sp.expand(tp)))
for kk, DD, name in [(sp.Symbol("b")**2, sp.Symbol("b")**2, "d = 4"), (sp.Symbol("b"), sp.Symbol("b")**2, "d = 3")]:
    e = sp.expand(sp.simplify(dec.subs({k: kk, D: DD}))); print(name, ": beta - R(beta) =", sp.collect(sp.simplify(e), be))
print("d = 4, b = 2:", sp.nsimplify(sp.expand(dec.subs({k: 4, D: 4}))))
