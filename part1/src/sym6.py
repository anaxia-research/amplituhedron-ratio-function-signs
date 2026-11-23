import sympy as sp, pickle, time
t0=time.time()
c1,c2,c3,c4,c5,c6=sp.symbols('c1:7',positive=True)
r1=1/((c1+c6)*(c1+c2)*(c1-c3)*(c1+c4)*(c1-c5)); r2=1/((c1+c2)*(c2+c3)*(c2-c4)*(c2+c5)*(c2-c6))
r3=1/((c2+c3)*(c3+c4)*(c3-c5)*(c3+c6)*(c3-c1)); r4=1/((c3+c4)*(c4+c5)*(c4-c6)*(c1+c4)*(c4-c2))
r5=1/((c4+c5)*(c5+c6)*(c1-c5)*(c2+c5)*(c3-c5)); r6=1/((c5+c6)*(c1+c6)*(c2-c6)*(c3+c6)*(c4-c6))
u=(c1+c2)*(c4+c5)/((c2+c5)*(c1+c4)); v=(c2+c3)*(c5+c6)/((c2+c5)*(c3+c6)); w=(c1+c6)*(c3+c4)/((c1+c4)*(c3+c6))
Om=sp.factor(sp.together(r1+r3+r5))
print('Omega:',Om); print(time.time()-t0)
assert sp.simplify(r1+r3+r5-r2-r4-r6)==0
A=r1-r2+r3; B=r2-r3+r4; C=r3-r4+r5
kA=sp.factor(sp.cancel(A*(u-1)*(v-1)/(Om*u*v)))
print('kA done',time.time()-t0)
kB=sp.factor(sp.cancel(B*(v-1)*(w-1)/(Om*v*w)))
kC=sp.factor(sp.cancel(C*(w-1)*(u-1)/(Om*w*u)))
print(time.time()-t0)
pickle.dump((Om,kA,kB,kC),open('k6.pkl','wb'))
for nm,k in [('kA',kA),('kB',kB),('kC',kC)]:
    n,d=sp.fraction(sp.together(k)); print(nm,'denominator factors:',sp.factor_list(d)[1][:5] if True else '')
