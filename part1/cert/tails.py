"""Explicit tail bounds (regions L and V), evaluated in ball arithmetic.
Region L: v <= V0=e^5, w >= W0=e^24:  H <= -c^2/2 + 2m(c+1) + m^2,  m=mu(V0).
Region V: v >= V0 (so w >= v >= V0):
  D >= 0:  B <= H <= KV_plus;   D < 0:  B <= KV_minus."""
from flint import arb, acb
from common import Z2, muL, lamL
b0 = arb(5); V0 = b0.exp(); c0 = arb(24)
m = muL(5.0)
QL = -c0**2/2 + 2*m*(c0+1) + m**2
out = []
out.append('Region L: m=mu(e^5)=%s; vertex of the quadratic at c=2m=%s < 24' % (m, 2*m))
out.append('  bound at c=24 (decreasing beyond): %s' % QL)
assert 2*m < c0 and QL < 0
lV = lamL(5.0)
gV = acb(1-1/V0).polylog(2).real
rho = 1/(1-1/(2*V0**2))
muu = rho*b0/(V0-rho)
KVp = (muu+lV)**2/2 + 2*rho*b0**2/(V0-rho) - Z2 - 2*gV + 2*lV*muu + 2*b0**2/(V0-1) + lV**2
P0 = (2*V0/(2*V0-1))*(b0+arb(2).log())*b0/(V0-1)
KVm = Z2 + lV**2/2 - 2*Z2 - 2*gV + b0**2/(V0-1) + (2*b0**2+2*b0)/(V0-1) + P0*(1+V0/(V0-1))
out.append('Region V: g(V0)=%s, lam(V0)=%s, rho=%s' % (gV, lV, rho))
out.append('  D>=0: H <= %s' % KVp)
out.append('  D<0 : B <= %s' % KVm)
assert KVp < 0 and KVm < 0
out.append('TAILS CERTIFIED')
print('\n'.join(out)); open('logs/tails.log','w').write('\n'.join(out)+'\n')
