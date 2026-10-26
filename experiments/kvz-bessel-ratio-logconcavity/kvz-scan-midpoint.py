from mpmath import mp, mpf, rf, factorial, log
mp.dps=50
def Sser(nu,w,N=300): return sum(w**k/(factorial(k)*rf(nu+1,k)) for k in range(N))
def logR(nu,a2,t): return log(Sser(nu,t*t/4))-log(Sser(nu,a2*a2*t*t/4))
for nu in [mpf(1)/10,mpf(1)/5]:
  for a2 in [2,3,4,5,10]:
    best=None
    for x in [mpf(k)/4 for k in range(2,40)]:
      for h in [mpf(1)/2,mpf(1),mpf(3)/2,2,3]:
        if x-h<=0: continue
        gap=(logR(nu,a2,x-h)+logR(nu,a2,x+h))/2-logR(nu,a2,x)
        if best is None or gap>best[0]: best=(gap,float(x),float(h))
    w=a2*a2*(best[1]+best[2])**2/4
    print(float(nu),a2,mp.nstr(best[0],5),best[1:],"max outer arg",float(w),flush=True)
