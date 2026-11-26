# Heat flow of the Wilson law exp(kappa cos theta): sign of the mixed partial of the log density in (t, theta).
# Claim: at t = 0 it equals (kappa/4) sin(theta) (3 + 2 kappa cos(theta)), nonnegative iff kappa <= 3/2, and the
# maximum principle then keeps it nonnegative for all t.
import numpy as np
from scipy.special import ive
n = np.arange(1, 80); th = np.linspace(0.01, np.pi-0.01, 600); CHI = np.sin(np.outer(n, th))/np.sin(th)
def logf(kappa, t): return np.log((n*(ive(n, kappa)/ive(1, kappa))*np.exp(-(n**2-1)*t/4))@CHI)
def minmixed(kappa, ts):
    out = []
    for t in ts:
        h = 1e-5; dt = (logf(kappa, t+h)-logf(kappa, max(t-h, 0)))/(h+min(h, t)); out.append((np.gradient(dt, th)/np.sin(th)).min())
    return np.array(out)
ts = np.concatenate([[0], np.geomspace(1e-3, 8, 40)])
for kappa in [0.5, 1.0, 1.5, 2.0, 3.0]:
    m = minmixed(kappa, ts); print("kappa %.1f  min over t of min over theta: %.4f  (t=0 theory %.4f)" % (kappa, m.min(), kappa/4*(3-2*kappa)))
