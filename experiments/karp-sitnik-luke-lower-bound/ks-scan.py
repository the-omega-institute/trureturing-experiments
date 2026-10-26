from mpmath import mp, mpf, hyper
import itertools, random
mp.dps=25
random.seed(1)
worst=(1,None)
vals=[mpf(v) for v in ('0.1','0.3','0.5','1','2','4','8')]
for sig in [mpf(s) for s in ('0.1','0.3','0.5','0.8','1','1.5','3')]:
  for a1,a2,b1,b2 in itertools.product(vals,repeat=4):
    if (b1-a1)+(b2-a2)<=0: continue
    if b1>a1 and b2>a2: continue   # proven region (for sigma>=1); focus on mixed signs
    for x in [mpf(v) for v in ('-0.9','-0.5','0.5','2','10','100')]:
        if 1+x*a1*a2/(b1*b2)<=0: continue
        L=hyper([sig,a1,a2],[b1,b2],-x)
        if not isinstance(L,type(mpf(1))): continue
        R=1/(1+x*a1*a2/(b1*b2))**sig
        r=(L-R)/abs(R)
        if r<worst[0]: worst=(r,(float(sig),float(a1),float(a2),float(b1),float(b2),float(x)))
print('min relative (F - bound):',mp.nstr(worst[0],6),'at (sigma,a1,a2,b1,b2,x)=',worst[1])
