import mpmath as mp, numpy as np, sympy as sp, pickle
mp.mp.dps=60
c=sp.symbols('c1:7',positive=True)
N,res=pickle.load(open('N6.pkl','rb'))
fN=sp.lambdify(c,N,'mpmath'); fk={k:sp.lambdify(c,res[k],'mpmath') for k in 'ABC'}
def full(cc):
    c1,c2,c3,c4,c5,c6=cc
    r=[1/((c1+c6)*(c1+c2)*(c1-c3)*(c1+c4)*(c1-c5)),1/((c1+c2)*(c2+c3)*(c2-c4)*(c2+c5)*(c2-c6)),
       1/((c2+c3)*(c3+c4)*(c3-c5)*(c3+c6)*(c3-c1)),1/((c3+c4)*(c4+c5)*(c4-c6)*(c1+c4)*(c4-c2)),
       1/((c4+c5)*(c5+c6)*(c1-c5)*(c2+c5)*(c3-c5)),1/((c5+c6)*(c1+c6)*(c2-c6)*(c3+c6)*(c4-c6))]
    u=(c1+c2)*(c4+c5)/((c2+c5)*(c1+c4)); v=(c2+c3)*(c5+c6)/((c2+c5)*(c3+c6)); w=(c1+c6)*(c3+c4)/((c1+c4)*(c3+c6))
    L=lambda x: mp.polylog(2,1-x)
    V=lambda a,b,cc_: (L(a)+L(b)+L(cc_)+mp.log(a)*mp.log(b)-mp.log(a)*mp.log(cc_)+mp.log(b)*mp.log(cc_)-2*mp.zeta(2))/2
    F=(r[0]+r[3])*V(u,v,w)+(r[1]+r[4])*V(v,w,u)+(r[2]+r[5])*V(w,u,v)
    Om=r[0]+r[2]+r[4]
    n=fN(*cc); kA,kB,kC=[fk[k](*cc)/n for k in 'ABC']
    mu=lambda x: x*mp.log(x)/(x-1)
    G2=L(u)+L(v)+L(w)-2*mp.zeta(2)+kA*mu(u)*mu(v)+kB*mu(v)*mu(w)+kC*mu(w)*mu(u)
    return F/Om, G2, Om, n/((c1+c2)*(c1+c4)*(c1+c6)*(c2+c3)*(c2+c5)*(c3+c4)*(c3+c6)*(c4+c5)*(c5+c6))
rng=np.random.default_rng(3); mx=0
for i in range(300):
    cc=[mp.mpf(x) for x in np.exp(rng.uniform(-9,9,6))]
    a,b,Om,Om2=full(cc)
    mx=max(mx,abs(a-b)/(1+abs(a)), abs(Om-Om2)/abs(Om))
print('max rel diff at 60 digits:',mp.nstr(mx,5))
