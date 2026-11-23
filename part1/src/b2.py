import numpy as np
from scipy.special import spence
from bclosed import Bc, H, lam
z2=np.pi**2/6
def Bs(u,v,w):
    # simplified: B = H + D*Lambda/l*, l*=log(min) if D>0 else log(max)
    e=u+v+w-1; D=e-2*u*v*w
    L=lam(u)*lam(v)*lam(w)
    mn=np.minimum(np.minimum(u,v),w); mx=np.maximum(np.maximum(u,v),w)
    l=np.where(D>0,np.log(mn),np.log(mx))
    return H(u,v,w)+D*L/l
rng=np.random.default_rng(1); N=1000000
u,v,w=np.exp(rng.uniform(-10,7,size=(3,N)))
e=u+v+w-1; R=(e>0)&(e*e>=4*u*v*w)
b1=Bc(u[R],v[R],w[R]); b2=Bs(u[R],v[R],w[R])
print('diff',np.max(np.abs(b1-b2)/(1+np.abs(b1))))
