import numpy as np
from scipy.special import spence
# Li2(1-x) = spence(x) for x>0 in scipy convention: spence(z) = int_1^z log t/(1-t) dt = Li2(1-z)
z2=np.pi**2/6
def L(x): return spence(x)
def brackets(c):
    c1,c2,c3,c4,c5,c6=c
    B={}
    B['1234']=c5+c6; B['1235']=c6-c4; B['6123']=c4+c5; B['1245']=c3+c6; B['1246']=c3-c5
    B['1256']=c3+c4; B['1345']=c6-c2; B['3461']=c2+c5; B['1356']=c4-c2; B['4561']=c2+c3
    B['2345']=c1+c6; B['2346']=c1-c5; B['2356']=c1+c4; B['2456']=c1-c3; B['3456']=c1+c2
    return B
def R(c):
    c1,c2,c3,c4,c5,c6=c
    r1=1/((c1+c6)*(c1+c2)*(c1-c3)*(c1+c4)*(c1-c5))
    r2=1/((c1+c2)*(c2+c3)*(c2-c4)*(c2+c5)*(c2-c6))
    r3=1/((c2+c3)*(c3+c4)*(c3-c5)*(c3+c6)*(c3-c1))
    r4=1/((c3+c4)*(c4+c5)*(c4-c6)*(c1+c4)*(c4-c2))
    r5=1/((c4+c5)*(c5+c6)*(c1-c5)*(c2+c5)*(c3-c5))
    r6=1/((c5+c6)*(c1+c6)*(c2-c6)*(c3+c6)*(c4-c6))
    return r1,r2,r3,r4,r5,r6
def uvw(c):
    c1,c2,c3,c4,c5,c6=c
    u=(c1+c2)*(c4+c5)/((c2+c5)*(c1+c4))
    v=(c2+c3)*(c5+c6)/((c2+c5)*(c3+c6))
    w=(c1+c6)*(c3+c4)/((c1+c4)*(c3+c6))
    return u,v,w
def V1(u,v,w):
    return 0.5*(L(u)+L(v)+L(w)+np.log(u)*np.log(v)-np.log(u)*np.log(w)+np.log(v)*np.log(w)-2*z2)
def P(c):
    r=R(c); u,v,w=uvw(c)
    tree=r[0]+r[2]+r[4]
    F=(r[0]+r[3])*V1(u,v,w)+(r[1]+r[4])*V1(v,w,u)+(r[2]+r[5])*V1(w,u,v)
    return F,tree,(r[1]+r[3]+r[5])
rng=np.random.default_rng(1)
N=2000000
c=np.exp(rng.uniform(-8,0,size=(6,N)))
F,T,T2=P(c)
print('tree check', np.max(np.abs(T-T2)/np.abs(T)), 'tree min', T.min()>0)
G=F/T
print('G max', G.max(), 'G min', G.min(), 'frac >=0', np.mean(G>=0))
i=np.argmax(G); print('argmax c', c[:,i], uvw(c[:,i]))
c=rng.uniform(0,1,size=(6,N)); F,T,_=P(c); G=F/T; print('uniform G max',G.max(), uvw(c[:,np.argmax(G)]))
