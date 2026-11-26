# A latent variable for the SU(2) heat kernel. Product formula: K_t(theta)/K_t(pi/2) = prod_m (1 + x/gamma_m), x = cos theta, gamma_m > 1.
# Each factor is (1 - p_m) + p_m (1 + x) with p_m = 1/gamma_m, so K_t is a mixture over M = sum of Bernoulli(p_m) of the laws (1 + x)^M,
# and with y = cos^2(theta/2) the law given M is Beta(M + 3/2, 3/2). Check the product, the law of M, and an exact sampler.
import numpy as np
from scipy.optimize import brentq
from scipy.special import betaln, gammaln
def roots(t, nr):
    N = nr+40; C = np.cosh((2*np.arange(1, N+1)-1)*t/4); f = lambda g: np.sum(1/(C-g)); return np.array([brentq(f, C[m]*(1+1e-13)+1e-13, C[m+1]*(1-1e-13)-1e-13) for m in range(nr)])
def K(t, th): n = np.arange(1, 400)[:, None]; return (n*np.exp(-(n*n-1)*t/4)*np.sin(n*th)/np.sin(th)).sum(0)
rng = np.random.default_rng(2)
for t in (4.0, 2.0, 1.0, 0.5):
    nr = int(60/t)+10; g = roots(t, nr); p = 1/g; th = np.linspace(0.05, 3.09, 9); prod = np.prod(1+np.cos(th)[None, :]*p[:, None], 0); err = np.max(np.abs(prod/(K(t, th)/K(t, np.array([np.pi/2])))-1))
    pm = np.array([1.0])
    for q in p: pm = np.convolve(pm, [1-q, q])                                       # Poisson binomial law of the count
    m = np.arange(pm.size); lw = np.log(np.clip(pm, 1e-300, None))+m*np.log(2)+betaln(m+1.5, 1.5); w = np.exp(lw-lw.max()); w /= w.sum()   # weight by the mass of (1+x)^M under Haar
    M = rng.choice(m, 200000, p=w); y = rng.beta(M+1.5, 1.5); cos = 2*y-1
    print("t %.1f beta %.1f  product formula error %.1e  largest p %.4f  E[M] %.3f  sd[M] %.3f  sampled E cos %.5f  exact %.5f" % (t, 2/t, err, p[0], (m*w).sum(), np.sqrt((m*m*w).sum()-(m*w).sum()**2), cos.mean(), np.exp(-3*t/4)))
