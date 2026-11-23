"""Independent consistency test of the core certificate: random points are pushed
down the same bisection tree; at each leaf we compare the certified bound with
B evaluated in floating point at the points of that leaf."""
import sys, math, numpy as np
sys.path.insert(0, '../src')
from core import box_bound, split, LU0, LV0, LW0
from b2 import Bs
rng = np.random.default_rng(11); N = 400000
a = rng.uniform(LU0, LV0, N); b = rng.uniform(LU0, LV0, N); c = rng.uniform(LU0, LW0, N)
P = np.stack([a, b, c], 1)
P = P[(P[:,0] <= P[:,1]) & (P[:,1] <= P[:,2])]
u, v, w = np.exp(P.T); e = u+v+w-1; Dl = e*e-4*u*v*w
inside = (e >= 0.01) & (Dl >= 0)
P = P[inside]; Btrue = Bs(*np.exp(P.T))
print('test points in the core region:', len(P))
bad = 0; checked = 0
stack = [((LU0, LV0, LU0, LV0, LU0, LW0), np.arange(len(P)))]
while stack:
    box, idx = stack.pop()
    kind, val = box_bound(*box)
    if kind == 'fail':
        b1, b2 = split(box)
        for bb in (b1, b2):
            m = np.all((P[idx] >= np.array(bb[0::2])) & (P[idx] <= np.array(bb[1::2])), axis=1)
            stack.append((bb, idx[m]))
        continue
    if len(idx) == 0: continue
    if kind in ('empty', 'edge', 'outside', 'Binf'):
        bad += len(idx) if kind != 'Binf' else 0
        if kind != 'Binf': print('point in discarded leaf', kind, box, len(idx))
        continue
    ub = float(val.upper()); checked += len(idx)
    if np.any(Btrue[idx] > ub + 1e-12):
        bad += 1; print('bound violated', kind, box, Btrue[idx].max(), ub)
print('points checked against leaf bounds:', checked, ' violations:', bad)
