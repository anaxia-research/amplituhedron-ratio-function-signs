"""Rigorous ball-arithmetic enclosures (python-flint / Arb) of the one-variable
functions used in the proof:
    f(x)   = Li2(1-x)                (decreasing)
    mu(x)  = x log x /(x-1)          (increasing, concave, mu(1)=1)
    lam(x) = log x /(x-1)            (decreasing, lam(1)=1)
    dmu(x) = mu'(x) = (x-1-log x)/(x-1)^2   (decreasing, dmu(1)=1/2)
Points are given either by an exact float log-coordinate a (x=e^a) or by an
exact float x.  a=-inf means x=0.  Every returned ball contains the true value."""
import math
from flint import arb, acb, ctx
ctx.prec = 128
Z2 = arb.pi()**2/6
_cache = {}
def _f_of_x(X, logX, one_minus_X):
    """Li2(1-x) for x>0 given x, log x and 1-x as balls.
    x < 1/2: Euler reflection  Li2(1-x) = zeta2 - log(x) log(1-x) - Li2(x)."""
    if X < arb('0.5'):
        return Z2 - logX*one_minus_X.log() - acb(X).polylog(2).real
    return acb(one_minus_X).polylog(2).real
def _key(kind, t): return (kind, t)

def fL(a):
    """f(e^a)"""
    k = ('fL', a)
    if k in _cache: return _cache[k]
    if a == -math.inf: r = Z2
    elif a == 0.0: r = arb(0)
    else:
        r = _f_of_x(arb(a).exp(), arb(a), -arb(a).expm1())
    _cache[k] = r; return r
def muL(a):
    k = ('muL', a)
    if k in _cache: return _cache[k]
    if a == -math.inf: r = arb(0)
    elif a == 0.0: r = arb(1)
    else:
        A = arb(a); r = A/(-(-A).expm1())
    _cache[k] = r; return r
def lamL(a):
    k = ('lamL', a)
    if k in _cache: return _cache[k]
    if a == 0.0: r = arb(1)
    else:
        A = arb(a); r = A/A.expm1()
    _cache[k] = r; return r
def xL(a):
    k = ('xL', a)
    if k in _cache: return _cache[k]
    r = arb(0) if a == -math.inf else arb(a).exp()
    _cache[k] = r; return r
# linear-coordinate versions (exact float x >= 0)
def f(x):
    k = ('f', x)
    if k in _cache: return _cache[k]
    if x == 0.0: r = Z2
    elif x == 1.0: r = arb(0)
    else:
        X = arb(x); r = _f_of_x(X, X.log(), 1-X)
    _cache[k] = r; return r
def mu(x):
    k = ('mu', x)
    if k in _cache: return _cache[k]
    if x == 0.0: r = arb(0)
    elif x == 1.0: r = arb(1)
    else:
        X = arb(x); r = X*X.log()/(X-1)
    _cache[k] = r; return r
def lam(x):
    k = ('lam', x)
    if k in _cache: return _cache[k]
    if x == 1.0: r = arb(1)
    else:
        X = arb(x); r = X.log()/(X-1)
    _cache[k] = r; return r
def dmu(x):
    k = ('dmu', x)
    if k in _cache: return _cache[k]
    if x == 1.0: r = arb(1)/2
    else:
        X = arb(x); y = X-1
        r = (y - X.log())/(y*y)
    _cache[k] = r; return r
def g(x):
    """g(x)=Li2(1-1/x) (increasing)"""
    X = arb(x); return acb(1-1/X).polylog(2).real
