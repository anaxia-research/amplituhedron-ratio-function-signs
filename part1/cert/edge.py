"""Edge lemma (region e <= e0).  Certifies
  (i)  dH0/dw (v,w') = -lam(w') + mu(v) mu'(w') <= -k0  on
       K = {0<=v<=(1+e0)/2, 1-v <= w' <= 1-v+e0},
  (ii) the final constant  C = -(k0-eta2)(1-theta) + eta1 < 0,
so that H <= C*e < 0 whenever u<=v<=w, 0<e<=e0, Delta>=0."""
import sys, math
from flint import arb
from common import f, mu, lam, dmu
E0 = 0.01; K0 = arb('0.83')
def check_slope(n=2048):
    vmax = (1+E0)/2
    worst = None; boxes = 0
    for i in range(n):
        v1 = vmax*i/n; v2 = vmax*(i+1)/n
        wlo = 1-v2; whi = 1-v1+E0
        ub = -lam(whi) + mu(v2)*dmu(wlo)
        boxes += 1
        assert ub < -K0, (v1, v2, ub)
        if worst is None or ub.upper() > worst.upper(): worst = ub
    return boxes, worst
if __name__ == '__main__':
    out = []
    nb, worst = check_slope()
    out.append('slope check: %d boxes, max upper bound of dH0/dw on K = %s < -k0 = -%s' % (nb, worst.upper(), K0))
    e0 = arb(E0)
    theta = 1/(2*(1-e0).sqrt())            # u <= theta*e
    umax = theta*e0
    eta2 = (umax*umax.log().__neg__()/(1-umax))/2   # mu(umax)/2
    eta1 = e0*(2/e0**2).log()**2/(4*(1-umax)*(1-e0))
    C = -(K0-eta2)*(1-theta) + eta1
    out.append('theta = %s' % theta)
    out.append('eta1 = %s' % eta1)
    out.append('eta2 = %s' % eta2)
    out.append('C = -(k0-eta2)(1-theta)+eta1 = %s' % C)
    assert C < 0
    out.append('EDGE LEMMA CERTIFIED: H <= C*e < 0 on {u<=v<=w, 0<e<=%g, Delta>=0}' % E0)
    print('\n'.join(out))
    open('logs/edge.log','w').write('\n'.join(out)+'\n')
