"""Region S: u <= U0=e^-7, v <= V0=e^5, w <= W0=e^24, e > e0.
There H <= Phi(v,w) := H0(v,w) + mu(U0)*(mu(v)+mu(w)-1)^+, and v+w >= 1+e0-U0.
We certify Phi < 0 on {0<=v<=w, v<=V0, w<=W0, v+w >= 1+sigma}, sigma=e0-U0,
by adaptive bisection in log coordinates (v may be 0)."""
import math, time
from flint import arb
from common import fL, muL, xL, f, mu, Z2
E0 = 0.01; LU0 = -7.0; LV0 = 5.0; LW0 = 24.0
U0 = arb(LU0).exp(); MU0 = muL(LU0)
SIG = arb(E0) - U0
THRESH = arb('-0.001')
def down(t):  # a float <= every point of ball t
    return math.nextafter(float(t.lower()), -math.inf)
def bound(p1, p2, q1, q2):
    # sorted part: v <= w
    p2 = min(p2, q2); q1 = max(q1, p1)
    if p1 > p2 or q1 > q2: return None   # empty
    vlo, vhi, wlo, whi = xL(p1), xL(p2), xL(q1), xL(q2)
    # constraint v+w >= 1+sigma
    if (vhi + whi < 1 + SIG): return None
    t = 1 + SIG - whi
    fv = fL(p1)
    if t > vlo:
        x = max(down(t), 0.0); fv = f(x) if x > 0 else fv
    t = 1 + SIG - vhi
    fw = fL(q1)
    if t > wlo:
        x = down(t); fw = f(x)
    mv, mw = muL(p2), muL(q2)
    ex = mv + mw - 1
    pos = MU0*ex if not (ex < 0) else arb(0)
    return fv + fw - Z2 + mv*mw + pos
def run():
    t0 = time.time()
    stack = [(-math.inf, LV0, math.log(0.5), LW0)]
    done = 0; empty = 0; worst = None; maxdepth = 0
    while stack:
        p1, p2, q1, q2 = stack.pop()
        b = bound(p1, p2, q1, q2)
        if b is None: empty += 1; continue
        if b < THRESH:
            done += 1
            if worst is None or b.upper() > worst: worst = b.upper()
            continue
        # bisect the direction with the larger contribution to the enclosure width
        dv = (fL(p1) - fL(p2)) + (muL(p2) - muL(p1))*(muL(q2) + MU0)
        dw = (fL(q1) - fL(q2)) + (muL(q2) - muL(q1))*(muL(p2) + MU0)
        if dv.upper() >= dw.upper():
            m = p2 - 8.0 if p1 == -math.inf else (p1+p2)/2
            stack += [(p1, m, q1, q2), (m, p2, q1, q2)]
        else:
            m = (q1+q2)/2; stack += [(p1, p2, q1, m), (p1, p2, m, q2)]
        if (p2-p1 < 1e-7 or p2 < -700) and q2-q1 < 1e-7: raise RuntimeError('no convergence near %s' % ((p1,p2,q1,q2),))
    return done, empty, worst, time.time()-t0
if __name__ == '__main__':
    d, e, w, T = run()
    msg = ['Region S (u<=e^-7, v<=e^5, w<=e^24, e>0.01): Phi < -0.001 certified',
           'certified boxes: %d, empty boxes: %d, max upper bound: %s, time %.1fs' % (d, e, w, T)]
    print('\n'.join(msg)); open('logs/small_u.log','w').write('\n'.join(msg)+'\n')
