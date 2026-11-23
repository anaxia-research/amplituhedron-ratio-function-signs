"""Two-loop six-point MHV amplitude in the MHV positive region, in y-coordinates.
All quantities follow Goncharov-Spradlin-Vergu-Volovich (GSVV) and
Dixon-von Hippel-McLeod-Trnka (DvHMT) conventions."""
import mpmath as mp
def ell(n,x): return (mp.polylog(n,x)-(-1)**n*mp.polylog(n,1/x))/2
def L4(xp,xm):
    l=mp.log(xp*xm); s=l**4/384
    for m,df in enumerate([1,2,8,48]):
        s+=(-1)**m/mp.mpf(df)*l**m*(ell(4-m,xp)+ell(4-m,xm))
    return s
def uvw(a,b,c):
    u=a*(1-b)*(1-c)/((1-a*b)*(1-a*c)); v=b*(1-c)*(1-a)/((1-b*c)*(1-b*a)); w=c*(1-a)*(1-b)/((1-c*a)*(1-c*b))
    return u,v,w
def parts(a,b,c):
    a,b,c=map(mp.mpf,(a,b,c))
    trip=((a,b,c),(b,c,a),(c,a,b))
    Ls=[];Qs=[];Ys=[]
    for x,y,z in trip:
        tp=x*(1-y*z)/(1-x); tm=(1-y*z)/(y*z*(1-x))
        Ls.append(mp.re(L4(-tp,-tm)))
        zz=(1-x)*(1-x*y*z)/(x*(1-y)*(1-z))
        Qs.append(-mp.polylog(4,-zz)/2); Ys.append(-mp.polylog(2,-zz))
    d=mp.log(a*b*c); Y=sum(Ys); z2=mp.zeta(2); z4=mp.zeta(4)
    E2=sum(Ls)+sum(Qs)+d**4/384+z2*d**2/8+5*z4/4+z2*Y/2
    R=E2-z2*Y/2-Y**2/8
    return dict(L=Ls,Q=Qs,Y=Y,delta=d,E1=Y/2,E2=E2,R=R)
