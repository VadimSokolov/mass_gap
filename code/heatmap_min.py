# Probe for the "normalised minimum 0.67 on [0.04, 3.41]" of the mixed partial heat map in book_figs5.py (block 1 only).
import numpy as np, os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "book_figs5.py")).read(); a = src.index("# 1. mixed partial"); b = src.index("fig, ax", a)
pre = src[:src.index("# 1.")]; exec(pre.replace("matplotlib.use", "(lambda *a: None)")) if False else None
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(src[:src.index("# 1. mixed partial")]); exec(src[a:b])
sel = (T <= 3.41)
print("raw min of m on t in [0.04,3.41]: %.4f  (at t %.3f, theta %.3f); m.max %.4f" % (m[sel].min(), T[sel][m[sel].argmin()], TH[sel][m[sel].argmin()], m[sel].max()))
print("min over theta of m/max_theta m, minimised over t: %.4f" % (m[sel].reshape(-1, 300).min(1)/m[sel].reshape(-1, 300).max(1)).min())
mt = m.min(1); print("per t min over theta, min over t<=3.41: %.4f; ratio min/ (value at theta=pi/2) min over t: %.4f" % (mt[t <= 3.41].min(), (m[:, 150]**-1*mt).min()))
print("min over theta of m per t: t=0.04 %.4f, 0.5 %.4f, 1 %.4f, 2 %.4f, 3.41 %.4f" % tuple(np.interp([0.04, .5, 1, 2, 3.41], t, mt)))
