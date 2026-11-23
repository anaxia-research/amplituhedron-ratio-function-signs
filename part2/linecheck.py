import mpmath as mp

from r62 import R62
mp.mp.dps=60
K=3000
def H(word,z):
    # word as string like '0,0,1'; trailing 1
    w=[int(x) for x in word.split(',')]
    c=[mp.mpf(0)]*(K+1)
    for k in range(1,K+1): c[k]=mp.mpf(1)/k   # H_1
    for letter in reversed(w[:-1]):
        if letter==0: c=[c[k]/k if k>0 else 0 for k in range(K+1)]
        else:
            n=[mp.mpf(0)]*(K+1); s=mp.mpf(0)
            for k in range(1,K+1):
                s+=c[k-1]; n[k]=s/k
            c=n
    return sum(c[k]*z**k for k in range(1,K+1))
z2,z4=mp.zeta(2),mp.zeta(4)
u=mp.mpf('0.3'); z=1-u
def Y(u,v,w): return sum(mp.polylog(2,1-x)+mp.log(x)**2/2 for x in (u,v,w))
def E2(u,v,w):
    y=Y(u,v,w); return R62(u,v,w)+z2*y/2+y**2/8
# fit polynomial in Lv,Lw with monomials up to total degree 4
import itertools
mons=[(i,j) for i in range(5) for j in range(5) if i+j<=4 and i>=j]
pts=[(mp.mpf(10)**-a,mp.mpf(10)**-b) for a,b in [(20,21),(22,25),(24,23),(26,30),(21,28),(30,30),(25,25),(23,27),(28,22),(32,24),(27,33),(35,35),(29,31),(33,29),(31,36),(26,26),(40,30),(38,38),(36,32),(34,40)]]
def fit(F):
    rows=[];rhs=[]
    for v,w in pts:
        Lv,Lw=mp.log(1/v),mp.log(1/w)
        rows.append([ (Lv**i*Lw**j+(Lv**j*Lw**i if i!=j else 0)) for i,j in mons]); rhs.append(F(u,v,w))
    A=mp.matrix(rows); b=mp.matrix(rhs)
    sol=mp.lu_solve(A.T*A,A.T*b)
    return dict(zip(mons,sol))
fr=fit(R62)
print('R: r11',mp.nstr(fr[(1,1)],15),'pred',mp.nstr(H('0,1',z)/4,15))
print('R: r10',mp.nstr(fr[(1,0)],15),'pred',mp.nstr((2*H('0,0,1',z)+H('1,0,1',z))/4,15))
pr00=(6*H('0,0,0,1',z)+3*H('0,1,0,1',z)+4*H('1,0,0,1',z)+2*H('1,1,0,1',z)-2*z2*(H('0,1',z)+H('1,1',z)))/4
print('R: r00',mp.nstr(fr[(0,0)],15),'pred',mp.nstr(pr00,15))
fe=fit(E2)
base=6*H('0,0,0,1',z)+2*H('0,0,1,1',z)+4*H('0,1,0,1',z)+3*H('0,1,1,1',z)+4*H('1,0,0,1',z)+2*H('1,0,1,1',z)+3*H('1,1,0,1',z)+3*H('1,1,1,1',z)
print('E: e00',mp.nstr(fe[(0,0)],15),'pred',mp.nstr((base+2*z2*(H('0,1',z)+H('1,1',z))+15*z4)/4,15))
print('E: e11',mp.nstr(fe[(1,1)],15),'pred',mp.nstr(H('0,1',z)/4,15),' e10',mp.nstr(fe[(1,0)],15),'pred',mp.nstr((2*H('0,0,1',z)+H('1,0,1',z))/4,15))
print('E: e22',mp.nstr(fe[(2,2)],10),' e40',mp.nstr(fe[(4,0)],10),' e20',mp.nstr(fe[(2,0)],12),'pred',mp.nstr((H('0,1',z)+H('1,1',z)+4*z2)/8,12))
