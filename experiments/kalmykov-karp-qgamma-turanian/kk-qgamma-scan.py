from mpmath import mp, mpf, log, ceil
mp.dps=30
def qpoch_inf(a,q,N): 
    p=mpf(1)
    for j in range(N): p*=(1-a*q**j)
    return p
def rg(x,q):
    N=int(200/(-log(q)))+50
    return qpoch_inf(q**x,q,N)/(qpoch_inf(q,q,N)*(1-q)**(1-x))
def delta(m,mu,a,b,q):
    return sum(rg(mu+a+k,q)*rg(mu+b+m-k,q)-rg(mu+k,q)*rg(mu+a+b+m-k,q) for k in range(m+1))
for q in [mpf("0.5"),mpf("0.9"),mpf("0.97")]:
  best=None
  for t in [mpf(k)/20 for k in range(1,20)]:
    for mu in [-t, -t+mpf('0.05'), mpf(0), mpf('0.5')]:
      if mu < max(-t,-1): continue
      for m in range(0,5):
        d=delta(m,mu,t,t,q)
        if best is None or d<best[0]: best=(d,float(t),float(mu),m)
  print(float(q),"min delta",mp.nstr(best[0],6),best[1:],flush=True)
