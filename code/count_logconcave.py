# Does one coalescence preserve log concavity?  Exact rational arithmetic, random log concave starting laws.
from fractions import Fraction as F
from math import comb, exp
import random
def Q(j, a, b): return F(comb(a, j)*comb(b+2, b-j), comb(a+b+2, b))
Z = lambda j: F(comb(2*j+2, j+1), (j+2)*2**j)
def coal(pi):
    out = [F(0)]*len(pi)
    for a, pa in enumerate(pi):
        for b, pb in enumerate(pi):
            for j in range(min(a, b)+1): out[j] += pa*pb*Q(j, a, b)
    return out
def ratio(v):
    while v and v[-1] == 0: v = v[:-1]
    return min([float(v[i]**2/(v[i-1]*v[i+1])) for i in range(1, len(v)-1)] or [9.0])
if __name__ == "__main__":
    random.seed(3); worst_nu = worst_pi = 9.0
    for t in range(400):
        inc = sorted([random.gauss(0, random.choice([.3, 1, 3])) for _ in range(random.randint(2, 14))], reverse=True); nu = [F(1)]
        for d in inc: nu.append(nu[-1]*F(exp(d)).limit_denominator(1000))
        pi = [nu[j]*Z(j) for j in range(len(nu))]; s = sum(pi); out = coal([p/s for p in pi]); worst_nu = min(worst_nu, ratio([out[j]/Z(j) for j in range(len(out))]))
        s = sum(nu); worst_pi = min(worst_pi, ratio(coal([p/s for p in nu])))
    print("nu_0 log concave: smallest ratio nu(j)^2/(nu(j-1)nu(j+1)) after coalescence over 400 laws: %.4f" % worst_nu)
    print("pi   log concave: smallest ratio pi(j)^2/(pi(j-1)pi(j+1)) after coalescence over 400 laws: %.4f" % worst_pi)
