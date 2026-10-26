# Exact rational certificate: Karp-Vishnyakova-Zhang 2408.01755 Conjecture 1 (log-concavity clause).
# nu1 = nu2 = 1/10, a1 = 1, a2 = 2. R(t) = I_nu(t)/I_nu(2t) = 2^(-nu) * S(t^2/4)/S(t^2), S(w) = sum w^k/(k! (nu+1)_k).
# Midpoint log-convexity at t = 3/2, 9/2, 15/2:  R(3/2) R(15/2) > R(9/2)^2.
from fractions import Fraction as F
nu=F(1,10)
def term(w,k):
    r=F(1)
    for j in range(k): r*= w/((j+1)*(nu+1+j))
    return r
def lower(w,N): return sum(term(w,k) for k in range(N))
def upper(w,N):
    s=lower(w,N); tN=term(w,N); rho=w/((N+1)*(nu+1+N))
    assert rho<1
    return s+tN/(1-rho)
N=60
x,m,y=F(3,2),F(9,2),F(15,2)
# need  S(x^2/4) S(y^2/4) S(m^2)^2  >  S(m^2/4)^2 S(x^2) S(y^2)
L=lower(x*x/4,N)*lower(y*y/4,N)*lower(m*m,N)**2
U=upper(m*m/4,N)**2*upper(x*x,N)*upper(y*y,N)
print("exact comparison L > U:", L>U, " quotient approx", float(L/U), " quotient > 1.02:", L/U > F(102,100))
