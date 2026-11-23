import sympy as sp, pickle, numpy as np
from scipy.special import spence
from explore import P as Pfull, uvw
c=sp.symbols('c1:7',positive=True)
N,res=pickle.load(open('N6.pkl','rb'))
fN=sp.lambdify(c,N,'numpy'); fk={k:sp.lambdify(c,res[k],'numpy') for k in 'ABC'}
z2=np.pi**2/6
def mu(x): return np.where(np.abs(x-1)<1e-12,1.0,x*np.log(x)/np.where(np.abs(x-1)<1e-12,1,x-1))
def G(cc):
    u,v,w=uvw(cc); n=fN(*cc)
    kA,kB,kC=[fk[k](*cc)/n for k in 'ABC']
    S=spence(u)+spence(v)+spence(w)-2*z2
    return S+kA*mu(u)*mu(v)+kB*mu(v)*mu(w)+kC*mu(w)*mu(u),u,v,w
rng=np.random.default_rng(7)
M=400000
cc=np.exp(rng.uniform(-9,9,size=(6,M)))
g,u,v,w=G(cc)
F,T,_=Pfull(cc); g2=F/T
ok=np.isfinite(g2)&(np.abs(T)>0)
print('identity check (kappa form vs R-invariant form): max rel diff', np.nanmax(np.abs(g-g2)[ok]/(1+np.abs(g[ok]))))
e=u+v+w-1
print('max G',g.max(),' max of G + min(e/10,1e-3):',(g+np.minimum(e/10,1e-3)).max())
# scatter data: subsample 2500 points, plus near-boundary points
idx=rng.choice(M,2500,replace=False)
c2=np.exp(rng.uniform(-9,9,size=(6,1500))); c2[4:]*=np.exp(rng.uniform(-12,-3,size=(2,1500)))
g3,u3,v3,w3=G(c2); e3=u3+v3+w3-1
E=np.concatenate([e[idx],e3]); Y=np.concatenate([-g[idx],-g3])
m=(E>1e-7)&(E<1e7)&(Y>0)
np.savetxt('scatter.dat',np.stack([E[m],Y[m]],1),fmt='%.5e',header='e minusG',comments='')
print('scatter points',m.sum(),'min -G/e',(Y[m]/np.minimum(E[m],1e9)).min())
