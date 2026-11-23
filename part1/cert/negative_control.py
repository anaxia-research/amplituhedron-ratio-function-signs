"""Negative control: with the Delta>=0 test switched off and threshold 0, the same
bisection must fail, because B > 0 at points with Delta < 0,
e.g. (u,v,w) = (0.2450, 0.3737, 0.4857) where B = +0.0031."""
import sys, math, numpy as np
from flint import arb
import core
core.CHECK_DELTA = False
core.THRESH = arb(0)
try:
    core.run((-1.5, -1.3, -1.0, -0.9, -0.75, -0.7))
    print('NEGATIVE CONTROL FAILED (no failure detected)')
except RuntimeError as ex:
    box = eval(str(ex).split('at ', 1)[1])
    u, v, w = (math.exp((box[0]+box[1])/2), math.exp((box[2]+box[3])/2), math.exp((box[4]+box[5])/2))
    sys.path.insert(0, '../src')
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        from b2 import Bs
    Bv = Bs(np.array([u]), np.array([v]), np.array([w]))[0]
    e = u+v+w-1
    print('negative control OK: without Delta>=0 the bisection cannot terminate;')
    print('  it stops at a box of width < 1e-6 around (u,v,w)=(%.6f, %.6f, %.6f): B=%.2e, Delta=%.4f' % (u, v, w, Bv, e*e-4*u*v*w))
