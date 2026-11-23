import numpy as np
from scipy.special import spence
from lp2 import Bfun
z2=np.pi**2/6
def lam(x):
    d=np.where(np.abs(x-1)<1e-13,1.0,x-1)
    return np.where(np.abs(x-1)<1e-13,1.0,np.log(x)/d)
def H(u,v,w): 
    mu=lambda x: x*lam(x)
    return spence(u)+spence(v)+spence(w)-2*z2+mu(u)*mu(v)+mu(v)*mu(w)+mu(w)*mu(u)
def Bc(u,v,w):
    e=u+v+w-1; D=e-2*u*v*w
    psi=[lam(v)*lam(w)/np.abs(u-1), lam(w)*lam(u)/np.abs(v-1), lam(u)*lam(v)/np.abs(w-1)]
    X=[u,v,w]; best=np.full(u.shape,np.inf)
    for x,p in zip(X,psi):
        rel=np.where(D>0, x<1, x>1)
        best=np.where(rel,np.minimum(best,p),best)
    return H(u,v,w)-np.abs(D)*best
if __name__=='__main__':
    rng=np.random.default_rng(4); N=2000000
    u,v,w=np.exp(rng.uniform(-10,7,size=(3,N)))
    e=u+v+w-1; R=(e>0)&(e*e>=4*u*v*w)
    b1,_=Bfun(u[R],v[R],w[R]); b2=Bc(u[R],v[R],w[R])
    print('max |diff|', np.nanmax(np.abs(b1-b2)/(1+np.abs(b1))), 'B max', b2.max())
    # extreme regime: u small, v=w huge
    for uu in [1e-2,1e-3,1e-5]:
        for L in [20,60,200]:
            X=np.exp(L)
            uu_=np.array([uu]); print('u',uu,'v=w=e^%d'%L,'H',H(uu_,np.array([X]),np.array([X]))[0],'B',Bc(uu_,np.array([X]),np.array([X]))[0])
