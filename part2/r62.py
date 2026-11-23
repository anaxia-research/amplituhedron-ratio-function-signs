import mpmath as mp
mp.mp.dps=40
def ell(n,x): return (mp.polylog(n,x)-(-1)**n*mp.polylog(n,1/x))/2
def L4(xp,xm):
    l=mp.log(xp*xm)
    s=l**4/mp.mpf(8*6*4*2)
    dfact=[1,2,8,48]
    for m in range(4):
        s+=(-1)**m/mp.mpf(dfact[m])*l**m*(ell(4-m,xp)+ell(4-m,xm))
    return s
def R62(u1,u2,u3):
    u=[mp.mpf(u1),mp.mpf(u2),mp.mpf(u3)]
    D=(sum(u)-1)**2-4*u[0]*u[1]*u[2]
    sq=mp.sqrt(D)
    xm=(sum(u)-1-sq)/(2*u[0]*u[1]*u[2]); xp=1/(u[0]*u[1]*u[2]*xm)
    s=0; J=0; S2=0
    for ui in u:
        s+=L4(ui*xp,ui*xm)-mp.polylog(4,1-1/ui)/2
        J+=ell(1,ui*xp)-ell(1,ui*xm)
        S2+=mp.polylog(2,1-1/ui)
    return mp.re(s-S2**2/8+J**4/24+mp.pi**2/12*J**2+mp.pi**4/72)
if __name__=='__main__':
    for u in [0.3,0.6]:
        # fit coefficients on v,w -> 0
        import itertools
        pts=[(mp.mpf(10)**-a,mp.mpf(10)**-b) for a,b in [(8,9),(10,9),(8,12),(12,12),(9,11),(11,8)]]
        rows=[];rhs=[]
        for v,w in pts[:4]:
            Lv,Lw=mp.log(1/v),mp.log(1/w)
            rows.append([Lv*Lw,Lv+Lw,1,(Lv**2+Lw**2)]); rhs.append(R62(u,v,w))
        sol=mp.lu_solve(mp.matrix(rows),mp.matrix(rhs))
        print(u,[mp.nstr(x,10) for x in sol],'Li2(1-u)/4=',mp.nstr(mp.polylog(2,1-u)/4,10))
    print('R(1/4,1/4,1/4)=',R62(0.25,0.25,0.25), 'R(0.1,.1,.1)',R62(.1,.1,.1))
