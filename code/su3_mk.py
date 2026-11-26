"""SU(3) Migdal Kadanoff decrement (b=2, d=4) versus the formula (1-b^-2) N^2 E[S]/4 = 9/16.
Wilson weight exp((beta/3) Re tr U). Coupling coordinate: beta_eff by matching the fundamental
character coefficient u_F to the Wilson value. Polson and Sokolov, Nov 2026."""
import numpy as np
M = 192
g = (np.arange(M)+0.5)*2*np.pi/M - np.pi
T1, T2 = np.meshgrid(g, g+0.31*2*np.pi/M, indexing="ij"); T3 = -T1-T2
TH = np.stack([T1, T2, T3])
E = np.exp(1j*TH)
vd = (E[0]-E[1])*(E[0]-E[2])*(E[1]-E[2])
haar = np.abs(vd)**2; haar /= haar.sum()
den = vd  # Weyl denominator det[e^{i(3-j)theta_k}] equals Vandermonde up to a phase
def char(p, q):
    lam = [p+q+2, q+1, 0]
    A = np.array([[np.exp(1j*l*TH[k]) for k in range(3)] for l in lam])  # 3x3xMxM
    num = (A[0,0]*(A[1,1]*A[2,2]-A[1,2]*A[2,1]) - A[0,1]*(A[1,0]*A[2,2]-A[1,2]*A[2,0]) + A[0,2]*(A[1,0]*A[2,1]-A[1,1]*A[2,0]))
    Ad = np.array([[np.exp(1j*l*TH[k]) for k in range(3)] for l in [2, 1, 0]])
    dn = (Ad[0,0]*(Ad[1,1]*Ad[2,2]-Ad[1,2]*Ad[2,1]) - Ad[0,1]*(Ad[1,0]*Ad[2,2]-Ad[1,2]*Ad[2,0]) + Ad[0,2]*(Ad[1,0]*Ad[2,1]-Ad[1,1]*Ad[2,0]))
    return num/dn
chiF = char(1, 0)
print("check dim F:", np.real(np.sum(haar*np.abs(chiF)**2)), "(orthonormal: 1)")
retr = np.cos(T1)+np.cos(T2)+np.cos(T3)
def uF_of_logw(logw):
    w = np.exp(logw-logw.max())*haar
    return np.real(np.sum(w*np.conj(chiF)))/(3*np.real(np.sum(w)))
bgrid = np.linspace(0.5, 60, 600)
uW = np.array([uF_of_logw((b/3)*retr) for b in bgrid])
def beta_eff(u): return np.interp(u, uW, bgrid)
P = 45
reps = [(p, q) for p in range(P+1) for q in range(P+1-p)]
CH = {r: char(*r) for r in reps}
def C2(p, q): return (p*p+q*q+p*q+3*p+3*q)/3
def dim(p, q): return (p+1)*(q+1)*(p+q+2)/2
def R(beta0):
    # heat kernel with u_F matched to Wilson at beta0: u_lambda = exp(-t C2(lambda)), C2(F)=4/3
    uF0 = uF_of_logw((beta0/3)*retr); t = -np.log(uF0)/C2(1, 0)
    f = sum(dim(*r)*np.exp(-t*C2(*r))*CH[r] for r in reps).real
    f = np.clip(f, 1e-300, None)
    u1 = uF_of_logw(4*np.log(f))
    return beta_eff(u1**4)
for b0 in [4, 8, 12, 16, 20, 24]:
    r = R(b0); print("beta=%5.1f  R=%.4f  decrement=%.4f" % (b0, r, b0-r))
print("formula (1-1/4)*9/12 = %.4f ; one loop SU(3) 11*9*ln2/(12 pi^2) = %.4f" % (9/16, 11*9*np.log(2)/(12*np.pi**2)))
