"""Core region C: e^-7 <= u <= v <= w, v <= e^5, w <= e^24, e >= e0=0.01, Delta >= 0.
Certifies B < 0 by adaptive bisection in log coordinates (a,b,c)=(log u,log v,log w).
On a box: H <= f(u_lo)+f(v_lo)+f(w_lo)-2z2+mu(u_hi)mu(v_hi)+mu(v_hi)mu(w_hi)+mu(w_hi)mu(u_hi);
if D>0 on the box: B = H - D lam(v)lam(w)/(1-u);  if D<0: B = H - |D| lam(u)lam(v)/(w-1)."""
import math, time, sys
from flint import arb
from common import fL, muL, lamL, xL, Z2
CHECK_DELTA = True
THRESH = arb("-0.001")   # certify B < THRESH (Theorem C); B < 0 follows
E0 = arb("0.01"); LU0 = -7.0; LV0 = 5.0; LW0 = 24.0
def box_bound(a1, a2, b1, b2, c1, c2):
    # restrict to the sorted part a<=b<=c
    a2 = min(a2, b2, c2); b1 = max(b1, a1); b2 = min(b2, c2); c1 = max(c1, b1)
    if a1 > a2 or b1 > b2 or c1 > c2: return 'empty', None
    ul, uh, vl, vh, wl, wh = xL(a1), xL(a2), xL(b1), xL(b2), xL(c1), xL(c2)
    eh = uh + vh + wh - 1
    if eh < E0: return 'edge', None              # e < e0 on the whole box
    el = ul + vl + wl - 1
    if el > E0: elo = el
    else: elo = E0                               # only points with e >= e0 matter
    pl = ul*vl*wl; ph = uh*vh*wh
    if CHECK_DELTA and eh*eh - 4*pl < 0: return "outside", None  # Delta < 0 on the whole box
    mu_ = (muL(a2), muL(b2), muL(c2))
    H = fL(a1) + fL(b1) + fL(c1) - 2*Z2 + mu_[0]*mu_[1] + mu_[1]*mu_[2] + mu_[2]*mu_[0]
    if H < THRESH: return 'H', H
    Dl = elo - 2*ph; Dh = eh - 2*pl
    if Dl > 0 and a1 < 0:
        B = H - Dl*lamL(b2)*lamL(c2)/(1-ul)
        if B < THRESH: return 'B+', B
    elif Dh < 0:
        if c2 <= 0: return 'Binf', None          # D<0 and all ratios <= 1: B = -infinity
        if c1 > 0:
            B = H + Dh*lamL(a2)*lamL(b2)/(wh-1)
            if B < THRESH: return 'B-', B
    return 'fail', H
def split(box):
    a1, a2, b1, b2, c1, c2 = box
    # weight each direction by its contribution to the enclosure width
    wts = []
    for (l, h, others) in ((a1, a2, (b2, c2)), (b1, b2, (a2, c2)), (c1, c2, (a2, b2))):
        if h - l <= 0: wts.append(-1); continue
        d = (fL(l) - fL(h)) + (muL(h) - muL(l))*(muL(others[0]) + muL(others[1]))
        wts.append(float(d.upper()))
    i = max(range(3), key=lambda k: wts[k])
    l, h = box[2*i], box[2*i+1]; m = (l+h)/2
    b1_ = list(box); b2_ = list(box); b1_[2*i+1] = m; b2_[2*i] = m
    return tuple(b1_), tuple(b2_)
def run(root, log=None, maxboxes=10**9):
    t0 = time.time(); stack = [root]; stats = {}; worst = None; n = 0
    while stack:
        box = stack.pop(); n += 1
        kind, val = box_bound(*box)
        if kind == 'fail':
            if min(box[1]-box[0], box[3]-box[2], box[5]-box[4]) < 1e-9 and max(box[1]-box[0], box[3]-box[2], box[5]-box[4]) < 1e-6:
                raise RuntimeError('no convergence at %r' % (box,))
            stack.extend(split(box)); continue
        stats[kind] = stats.get(kind, 0) + 1
        if val is not None and (worst is None or val.upper() > worst[0]): worst = (val.upper(), box)
        if log and n % 200000 == 0:
            print(n, stats, len(stack), '%.0fs' % (time.time()-t0), file=log, flush=True)
        if n > maxboxes: break
    return n, stats, worst, time.time()-t0
if __name__ == '__main__':
    root = (LU0, LV0, LU0, LV0, LU0, LW0)
    with open('logs/core_progress.log', 'w') as lg:
        n, stats, worst, T = run(root, lg)
    msg = ['Core region (e^-7<=u<=v<=w, v<=e^5, w<=e^24, e>=0.01, Delta>=0): B < -0.001 certified',
           'boxes examined: %d; leaf counts: %s' % (n, stats),
           'largest certified upper bound: %s on box %r' % worst, 'time %.0fs' % T]
    print('\n'.join(msg)); open('logs/core.log', 'w').write('\n'.join(msg)+'\n')
