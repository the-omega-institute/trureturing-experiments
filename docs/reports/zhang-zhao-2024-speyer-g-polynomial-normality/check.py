#!/usr/bin/env python3
"""Zhang–Zhao arXiv:2408.15506 Conjecture 4.8: exact checks of the identity
q(t) = g_{n,d}(t)/t = sum_j C(d-1,j) C(n-d-1,j) (1+t)^j and of the variance bound
sigma^2 >= (d-1)(n-d-1)/(4(n-2)) >= (d-1)/8 for 2 <= d <= n/2, n <= N (argv[1]),
plus real-rootedness of g_{n,d} by Sturm sequences (sympy, exact)."""
import sys
from fractions import Fraction as F
from math import comb, factorial
import sympy as sp
N=int(sys.argv[1]) if len(sys.argv)>1 else 40
t=sp.symbols('t')
def S(n,d,i): return F(factorial(n-i-1), factorial(d-i)*factorial(n-d-i)*factorial(i-1))
bad=0; checked=0; minratio=None
for n in range(4,N+1):
  for d in range(2,n//2+1):
    a,b=d-1,n-d-1
    # identity: coefficient of t^k in sum_j C(a,j)C(b,j)(1+t)^j equals S_{k+1}
    for k in range(0,a+1):
      if sum(comb(a,j)*comb(b,j)*comb(j,k) for j in range(k,a+1)) != S(n,d,k+1): bad+=1
    tot=sum(S(n,d,i) for i in range(1,d+1))
    mu=sum(i*S(n,d,i) for i in range(1,d+1))/tot
    var=sum((i-mu)**2*S(n,d,i) for i in range(1,d+1))/tot
    lb=F((d-1)*(n-d-1),4*(n-2))
    if not (var>=lb and lb>=F(d-1,8)): bad+=1
    r=var/lb; minratio=r if minratio is None or r<minratio else minratio
    if n<=24:
      g=sp.Poly(sum(sp.Rational(S(n,d,i).numerator,S(n,d,i).denominator)*t**(i-1) for i in range(1,d+1)),t)
      nreal=sp.count_roots(g,-sp.oo,-1)  # strictly below -1? count on (-oo,-1]
      at=g.eval(-1)
      if nreal!=d-1 or at==0: bad+=1
    checked+=1
print(f"pairs={checked} N={N} bad={bad} min_var_over_bound={float(minratio):.6f}")
sys.exit(1 if bad else 0)
