from mpmath import mp, mpf, besseli, log, diff
mp.dps=30
def d2(x,n1,a1,n2,a2): return diff(lambda t: log(besseli(n1,a1*t))-log(besseli(n2,a2*t)), x, 2)
pairs=[(0.5,0.5),(1.5,0.5),(2.5,0.5),(1.5,1.5),(2.5,1.5),(0.3,0.1),(0.5,0.1),(1,0.1),(0.2,0.2),(0.1,0.05),(0.5,0.25),(1,0.5),(0.15,0.1),(0.4,0.3)]
for n1,n2 in pairs:
  best=None
  for a2 in [1.2,1.5,2,3,5,10]:
    for k in range(1,121):
      x=mpf(k)/8
      v=d2(x,n1,1,n2,a2)
      if best is None or v>best[0]: best=(v,a2,float(x))
  print(n1,n2,"max d2",mp.nstr(best[0],5),"a2",best[1],"x",best[2],flush=True)
