"""Numerical evidence for the two-loop remainder-function conjecture and
sampling of the proven bounds. 40,000 points of the y-cube, sampled with
extra weight near every face; 25 significant digits."""
import mpmath as mp, random, sys
from multiprocessing import Pool
from e2lib import parts, uvw
mp.mp.dps=25
def sample(seed):
    rnd=random.Random(seed); out=[]
    for k in range(2500):
        y=[]
        for i in range(3):
            r=rnd.random(); e=rnd.choice([1,1,2,4,8])
            t=rnd.random()**e*0.999999+1e-7
            y.append(t if r<0.5 else 1-t)
        p=parts(*y)
        out.append((float(p['R']),float(p['E2']),float(p['Y']),float(min(p['L'])),y))
    return out
if __name__=='__main__':
    with Pool(4) as P: res=sum(P.map(sample,range(16)),[])
    n=len(res)
    minR=min(res); minE=min(res,key=lambda r:r[1]); minL=min(r[3] for r in res); minY=min(r[2] for r in res)
    lines=[f'points: {n}',
      f'min R6^(2) = {minR[0]:.6e} at y = {minR[4]}',
      f'points with R6^(2) <= 0: {sum(r[0]<=0 for r in res)}',
      f'min E^(2) = {minE[1]:.6f} at y = {minE[4]}',
      f'min L4 term = {minL:.3e} (all >= 0: {minL>=0})',
      f'min Y = {minY:.6f}']
    print('\n'.join(lines)); open('logs/evidence.log','w').write('\n'.join(lines)+'\n')
