# Extended scan of the bad cap of antipode_radius.py: maximum cap mass, its t, closure t, and the mass at beta >= 2.2 (t = 2/beta). Reuses eta() from antipode_radius.py.
import numpy as np, os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "antipode_radius.py")).read().split('print("  t')[0]; exec(src)
for k, lo, hi in ((4, 0.5, 3.2), (2, 0.5, 2.0)):
    ts = np.arange(lo, hi, 0.01); E = np.array([eta(t, k) for t in ts]); m = E[:, 1]; i = m.argmax()
    closed = ts[(E[:, 0] == 0) & (ts > ts[i])]
    print("k=%d: max cap mass %.3e at t=%.2f (eta %.4f); first t with empty cap after the maximum: %s" % (k, m[i], ts[i], E[i, 0], "%.2f" % closed[0] if closed.size else "none"))
    for b in (2.2, 2.5, 3.0, 4.0): print("   beta %.1f t %.3f eta %.4f mass %.2e" % ((b, 2/b)+eta(2/b, k)))
    for t in (1.45, 2.47): print("   t %.2f eta %.4f mass %.2e" % ((t,)+eta(t, k)))
