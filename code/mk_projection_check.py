import numpy as np
from mk_verify import make_R, setup
th, haar, n, CHI = setup(6000, 1500); R = make_R()
def full(beta0, steps):
    u = np.exp(-(2/beta0)*(n**2-1)/4); out = [beta0]
    for _ in range(steps):
        f = np.clip((n*u)@CHI, 1e-300, None); lf = 4*np.log(f); g = np.exp(lf-lf.max())
        c = (CHI@(g*haar))/n; c = np.clip(c/c[0], 0, None); u = c**4
        out.append(2/(-(4/3)*np.log(u[1])))
    return np.array(out)
for b0 in [3.0, 10.0]:
    F = full(b0, 8); P = [b0]
    for _ in range(8): P.append(R(P[-1]))
    print(b0, "max |full - projected| over 8 steps:", np.max(np.abs(F-np.array(P))))
    print("  full:", np.round(F, 4)); print("  proj:", np.round(P, 4))
