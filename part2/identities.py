"""High-precision checks of the identities used in the paper:
(1) the decomposition of E^(2) and R6^(2) agrees with the GSVV formula;
(2) J = -delta/2 = -(1/2) log(abc);
(3) the derivative formula for L4 (Theorem C) and L4(l)+L4(-l) = l^4/192;
(4) every L4 term is nonnegative at the sample points."""
import mpmath as mp, random
from e2lib import parts, uvw, ell, L4
from r62 import R62
mp.mp.dps=40
random.seed(11); worst=[0,0,0]
for k in range(200):
    a,b,c=[mp.mpf(random.random())**random.choice([1,2,3]) for _ in range(3)]
    p=parts(a,b,c); u,v,w=uvw(a,b,c)
    worst[0]=max(worst[0],abs(p['R']-R62(u,v,w))/(1+abs(p['R'])))
    J=0
    for x,y,z in ((a,b,c),(b,c,a),(c,a,b)):
        tp=x*(1-y*z)/(1-x); tm=(1-y*z)/(y*z*(1-x)); J+=ell(1,-tp)-ell(1,-tm)
    worst[1]=max(worst[1],abs(J+mp.log(a*b*c)/2))
    assert min(p['L'])>=0
def F(l,d): return mp.re(L4(-mp.e**((l+d)/2),-mp.e**((l-d)/2)))
for k in range(50):
    l=mp.mpf(random.uniform(-8,8)); d=mp.mpf(random.uniform(-8,8))
    lhs=mp.diff(lambda t:F(t,d),l); rhs=l**3/96*(1+mp.tanh((l+d)/4)/2+mp.tanh((l-d)/4)/2)
    worst[2]=max(worst[2],abs(lhs-rhs),abs(F(l,d)+F(-l,d)-l**4/192))
out=[f'decomposition vs GSVV, 200 points, max relative difference: {mp.nstr(worst[0],5)}',
     f'J + log(abc)/2, max absolute value: {mp.nstr(worst[1],5)}',
     f'L4 derivative formula and reflection identity, 50 points, max error: {mp.nstr(worst[2],5)}',
     'all L4 terms nonnegative: True']
print('\n'.join(out)); open('logs/identities.log','w').write('\n'.join(out)+'\n')
