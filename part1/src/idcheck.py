import numpy as np, mpmath as mp
mp.mp.dps=30
def f(x): return mp.polylog(2,1-x)
def mu(x): return x*mp.log(x)/(x-1)
def lam(x): return mp.log(x)/(x-1)
def g(x): return mp.polylog(2,1-1/x)
z2=mp.zeta(2)
rng=np.random.default_rng(0)
for _ in range(5):
    b=mp.mpf(rng.uniform(5,9)); t=mp.mpf(rng.uniform(0,3)); c=b+t
    u=mp.mpf(np.exp(rng.uniform(-5,b))); v=mp.e**b; w=mp.e**c
    H=f(u)+f(v)+f(w)-2*z2+mu(u)*mu(v)+mu(v)*mu(w)+mu(w)*mu(u)
    lv,lw,lu=lam(v),lam(w),lam(u)
    H2=-t**2/2+t*(mu(u)+lv)+2*b*mu(u)+f(u)-2*z2-g(v)-g(w)+mu(u)*(lv+lw)+b*(lv+lw)+lv*lw
    D=u+v+w-1-2*u*v*w
    B=H+D*lu*lv/(w-1)
    B2=-t**2/2+mu(u)*t+f(u)-2*z2-g(v)-g(w)+c*lv+lw*(mu(u)+b+lv)-mu(u)*lv+lu*lv+lu*lv*(u+v-2*u*v)/(w-1)
    print(float(H-H2), float(B-B2))
