# Certificate that zeta(t, x) = -d/dt [ t d/dx log K_t ] > 0 on -1 <= x <= 1 at a fixed time t0 (x = cos theta, K_t the SU(2) heat kernel with u_2 = exp(-3t/4)).
# Product formula: d/dx log K = S1 - S2/S1 with S_p = sum_n a_n^p, a_n = 1/(C_n + x), C_n = cosh((2n-1)t/4).  With D_p = sum_n C_n' a_n^p,
#   zeta = t (D2 - 2 D3/S1 + D2 S2/S1^2) - S1 + S2/S1.
# Every sum is a sum of positive terms decreasing in x, so on a box [x1, x2] it lies between its values at x2 and x1. Boxes are taken in 1 + x. The enclosure of zeta
# combines these with an allowance EPS for rounding and for the tail of the series; boxes are bisected until the lower bound is positive.
import numpy as np, sys
EPS = 1e-9
def sums(t, dlt, N):                 # dlt = 1 + x; C_n + x = 2 sinh^2(y_n/2) + dlt, with no cancellation at the antipode
    n = np.arange(1, N+1); y = (2*n-1)*t/4; Cp = (2*n-1)/4*np.sinh(y); a = 1/(2*np.sinh(y/2)**2+dlt)
    return np.array([a.sum(), (a*a).sum(), (Cp*a*a).sum(), (Cp*a**3).sum()])
def lower(t, x1, x2, N):
    hi = sums(t, x1, N)*(1+EPS)+EPS; lo = sums(t, x2, N)*(1-EPS); S1h, S2h, D2h, D3h = hi; S1l, S2l, D2l, D3l = lo
    return t*(D2l-2*D3h/S1l+D2l*S2l/S1h**2)-S1h+S2l/S1h
def certify(t, depth=60):
    N = int(np.ceil((180/t+1)/2)); stack = [(0.0, 2.0)]; nbox = 0; worst = np.inf; width = 2.0      # terms up to y = 45; boxes in dlt = 1 + x
    while stack:
        x1, x2 = stack.pop(); lb = lower(t, x1, x2, N)
        if lb > 0: nbox += 1; worst = min(worst, lb); width = min(width, x2-x1)
        elif x2-x1 < 2.0**-depth: return None
        else: m = 0.5*(x1+x2); stack += [(x1, m), (m, x2)]
    return N, nbox, worst, width
if __name__ == "__main__":
    for t0 in [float(a) for a in sys.argv[1:]] or [2.0, 0.5, 0.1, 0.02]:
        r = certify(t0); print("t0 = %-6g" % t0, "FAILED" if r is None else "certified: terms %d, boxes %d, least lower bound %.1e, smallest box %.1e" % r, flush=True)
