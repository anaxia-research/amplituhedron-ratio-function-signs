"""Floating-point sanity test of the region decomposition and of each analytic bound
(not part of the proof; the proof is the ball-arithmetic certificates)."""
import sys, numpy as np
sys.path.insert(0, '../src')
from scipy.special import spence
from b2 import Bs
from bclosed import H, lam
z2 = np.pi**2/6
mu = lambda x: x*lam(x)
rng = np.random.default_rng(5); N = 3000000
pts = []
for lo, hi in [(-12, 30), (-12, 2), (-3, 1)]:
    X = np.exp(rng.uniform(lo, hi, (3, N)))
    pts.append(np.sort(X, 0))
u, v, w = np.concatenate(pts, 1)
# near-edge oversampling
t = rng.uniform(0, 1, N); vv = rng.uniform(0, 0.5, N); s = 10**rng.uniform(-5, -1.5, N)
uu = s**2*t/(4*vv*(1-vv)+1e-300)
u = np.concatenate([u, uu]); v = np.concatenate([v, vv]); w = np.concatenate([w, 1-vv+s])
S = np.sort(np.stack([u, v, w]), 0); u, v, w = S
e = u+v+w-1; Dl = e*e-4*u*v*w; R = (e > 0) & (Dl >= 0)
u, v, w, e = u[R], v[R], w[R], e[R]
B = Bs(u, v, w); Hh = H(u, v, w)
print('points in region:', len(u), ' max B:', B.max())
E = e <= 0.01
Sr = (~E) & (u <= np.exp(-7)) & (v <= np.exp(5)) & (w <= np.exp(24))
L = (v <= np.exp(5)) & (w >= np.exp(24))
V = v >= np.exp(5)
C = (~E) & (~Sr) & (~L) & (~V)
print('counts E,S,L,V,C:', E.sum(), Sr.sum(), L.sum(), V.sum(), C.sum(), ' total', E.sum()+Sr.sum()+L.sum()+V.sum()+C.sum())
print('E: max H/e =', (Hh[E]/e[E]).max(), '(certified bound -0.157)')
H0 = spence(v)+spence(w)-z2+mu(v)*mu(w)
Phi = H0 + mu(np.exp(-7))*np.maximum(mu(v)+mu(w)-1, 0)
print('S: max H-Phi =', (Hh[Sr]-Phi[Sr]).max(), ' max Phi =', Phi[Sr].max())
c = np.log(w)
print('L: max H-(-c^2/2+2m(c+1)+m^2) =', (Hh[L]-(-c[L]**2/2+2*5.0339*(c[L]+1)+5.0339**2)).max())
D = e-2*u*v*w
print('V, D>=0: max H =', Hh[V & (D >= 0)].max() if (V & (D>=0)).any() else None, '(bound -4.17)')
print('V, D<0 : max B =', B[V & (D < 0)].max(), '(bound -3.89)')
print('C: max B =', B[C].max())
