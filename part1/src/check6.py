import sympy as sp, pickle
c1,c2,c3,c4,c5,c6=sp.symbols('c1:7',positive=True)
Om,kA,kB,kC=pickle.load(open('k6.pkl','rb'))
N=sp.expand(sp.fraction(sp.together(Om))[0])
res={}
for nm,k in [('A',kA),('B',kB),('C',kC)]:
    n,d=sp.fraction(sp.together(k)); n=sp.expand(n); d=sp.expand(d)
    # normalize so d == N
    q=sp.cancel(d/N); assert q.is_number, q
    n=sp.expand(n/q)
    one_minus=sp.Poly(sp.expand(N-n),c1,c2,c3,c4,c5,c6)
    coeffs=one_minus.coeffs()
    print(nm,'1-kappa numerator: #terms',len(coeffs),'min coeff',min(coeffs),'max',max(coeffs))
    print('   factor:',sp.factor(one_minus.as_expr()))
    res[nm]=n
u=(c1+c2)*(c4+c5)/((c2+c5)*(c1+c4)); v=(c2+c3)*(c5+c6)/((c2+c5)*(c3+c6)); w=(c1+c6)*(c3+c4)/((c1+c4)*(c3+c6))
eps=sp.factor(sp.together(1-u-v-w)); print('eps =',eps)
Dl=sp.factor(sp.together((1-u-v-w)**2-4*u*v*w)); print('Delta =',Dl)
# linear relation check
rel=sp.simplify(kA*u*v*(w-1)+kB*v*w*(u-1)+kC*u*w*(v-1)-(u-1)*(v-1)*(w-1))
print('relation:',rel)
pickle.dump((N,res),open('N6.pkl','wb'))
