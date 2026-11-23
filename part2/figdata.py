import mpmath as mp
from e2lib import parts, L4
mp.mp.dps=30
P='figdata/'
# F1: diagonal y=a=b=c, t=u=y/(1+y)^2
with open(P+'diag.dat','w') as f:
    f.write('t E1 E2 R\n')
    for k in range(1,200):
        y=mp.mpf(k)/200
        if k<20: y=mp.mpf(10)**(-6+ mp.mpf(k)*5/20)
        p=parts(y,y,y); t=y/(1+y)**2
        f.write(f'{float(t):.6e} {float(p["E1"]):.6e} {float(p["E2"]):.6e} {float(p["R"]):.6e}\n')
# F2: L4(l,delta)
with open(P+'l4.dat','w') as f:
    f.write('l d0 d2 d5\n')
    for k in range(-60,61):
        l=mp.mpf(k)/10
        vals=[]
        for d in (0,-2,-5):
            p=(l+d)/2; q=(l-d)/2
            vals.append(float(mp.re(L4(-mp.e**p,-mp.e**q))) if l!=0 else 0.0)
        f.write(f'{float(l):.3f} '+' '.join(f'{v:.6e}' for v in vals)+'\n')
# F3: approach to face y_u -> 1
with open(P+'face.dat','w') as f:
    f.write('eps R55 R19 Rbound\n')
    for k in range(0,41):
        e=mp.mpf(10)**(-1-mp.mpf(k)/4)
        r1=parts(1-e,0.5,0.5)['R']; r2=parts(1-e,0.1,0.9)['R']
        f.write(f'{float(e):.6e} {float(r1):.6e} {float(r2):.6e} {float(e*mp.log(1/e)**2):.6e}\n')
print(open(P+'diag.dat').read().splitlines()[-1]); 
p=parts(mp.mpf(1)-mp.mpf(10)**-12,1-mp.mpf(10)**-12,1-mp.mpf(10)**-12); print('near (1/4,1/4,1/4):',p['E2'],p['R'],p['Y'])
