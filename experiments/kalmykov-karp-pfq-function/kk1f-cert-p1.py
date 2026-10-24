# Exact rational certificate: Kalmykov-Karp Conjecture 1 function clause with a = (1), b = (1/2)^12, mu = 1, alpha = beta = 1/2.
from fractions import Fraction as F
from math import factorial
q=12; x=F(2000); N=10
def odd(k):  # (2k-1)!! for k>=0, with (-1)!!=1
    r=1
    for j in range(1,2*k,2): r*=j
    return r
# Gamma(n+1/2)=sqrt(pi)*(2n-1)!!/2^n
g_half=lambda n: F(odd(n),2**n)          # Gamma(n+1/2)/sqrt(pi)
# f(3/2): Gamma(5/2+n)/Gamma(2+n)^q = sqrt(pi)*g_half(n+2)/((n+1)!)^q
Ap=lambda n: g_half(n+2)/factorial(n+1)**q * x**n/factorial(n)
# f(1): Gamma(2+n)/Gamma(3/2+n)^q = (n+1)!/(pi^6 g_half(n+1)^q)
Rt=lambda n: F(factorial(n+1))/g_half(n+1)**q * x**n/factorial(n)
# f(2): Gamma(3+n)/Gamma(5/2+n)^q = (n+2)!/(pi^6 g_half(n+2)^q)
St=lambda n: F(factorial(n+2))/g_half(n+2)**q * x**n/factorial(n)
A=sum(Ap(n) for n in range(N)); R=sum(Rt(n) for n in range(N)); S=sum(St(n) for n in range(N))
rho=max(Ap(n+1)/Ap(n) for n in range(N-1,N+30))
Aup=A+Ap(N)/(1-rho)
piup=F(31416,10000)
lhs=piup**13*Aup**2; rhs=R*S
print("rho",float(rho),"A'",float(A),"R",float(R),"S",float(S))
print("exact comparison R*S > 3.1416^13*Aup^2:", rhs>lhs, " quotient approx", float(rhs/lhs), " quotient > 1.17:", rhs/lhs > F(117,100))
