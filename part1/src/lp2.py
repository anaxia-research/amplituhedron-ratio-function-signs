import numpy as np
from scipy.special import spence
z2=np.pi**2/6
def mu(x): return np.where(np.abs(x-1)<1e-12,1.0,x*np.log(x)/np.where(np.abs(x-1)<1e-12,1,x-1))
def Bfun(u,v,w):
    S=spence(u)+spence(v)+spence(w)-2*z2
    m=np.array([mu(u)*mu(v),mu(v)*mu(w),mu(w)*mu(u)]); a=np.array([u*v*(w-1),v*w*(u-1),u*w*(v-1)]); d=(u-1)*(v-1)*(w-1)
    best=np.full(u.shape,-np.inf); Vs=[]
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        ki=(d-a[j]-a[k])/a[i]
        val=S+m[j]+m[k]+m[i]*ki
        feas=ki<=1+1e-12
        best=np.where(feas,np.maximum(best,val),best); Vs.append((val,feas))
    return best,Vs
if __name__=='__main__':
    rng=np.random.default_rng(8); N=4000000
    u,v,w=np.exp(rng.uniform(-10,7,size=(3,N)))
    eps=1-u-v-w; Dl=eps**2-4*u*v*w
    B,Vs=Bfun(u,v,w)
    for nm,m in [('eps<0',eps<0),('eps<0,Delta>0',(eps<0)&(Dl>0)),('eps<0,Delta<0',(eps<0)&(Dl<0)),('eps>0',eps>0)]:
        print(nm, m.sum(), 'B max', B[m].max(), 'frac B>=0', np.mean(B[m]>=0))
    m=(eps<0)&(Dl>0)
    i=np.argmax(B[m]); print('argmax', u[m][i],v[m][i],w[m][i])
