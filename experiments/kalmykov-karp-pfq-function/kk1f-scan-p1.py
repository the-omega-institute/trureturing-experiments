# mpmath scan of min_x Delta_f / f(3/2;x)^2 for a = (a0), b = (1/2)^q, mu = 1, alpha = beta = 1/2.
from mpmath import mp, mpf, gamma, factorial
mp.dps=80
def f(mu,a,b,q,x,N=400):
    return sum(gamma(a+mu+n)/gamma(b+mu+n)**q * x**n/factorial(n) for n in range(N))
h=mpf(1)/2
for a in [mpf(1),h,mpf(3)/2,mpf(2)]:
  for q in [8,10,12,14,16,18,20]:
    best=None
    for e in range(0,141):
        x=mpf(10)**(mpf(e)/20)
        A=f(1+h,a,h,q,x); D=A*A-f(1,a,h,q,x)*f(2,a,h,q,x)
        r=D/(A*A)
        if best is None or r<best[0]: best=(r,x)
    print(float(a),q, mp.nstr(best[0],6), mp.nstr(best[1],6), flush=True)
