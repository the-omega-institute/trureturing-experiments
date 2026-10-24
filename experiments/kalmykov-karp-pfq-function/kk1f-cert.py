from fractions import Fraction as F
from math import factorial
q=12; x=F(5000); N=8
def dfact(n):  # (2n+1)!!
    r=1
    for k in range(1,2*n+2,2): r*=k
    return r
# Gamma(3/2+n)=sqrt(pi)*(2n+1)!!/2^(n+1); Gamma(5/2+n)=sqrt(pi)*(2n+3)!!/2^(n+2); Gamma(2+n)=(n+1)!
A=sum(x**n/(factorial(n)*factorial(n+1)**q) for n in range(N))
R=sum(x**n/(factorial(n)*F(dfact(n),2**(n+1))**q) for n in range(N))
S=sum(x**n/(factorial(n)*F(dfact(n+1),2**(n+2))**q) for n in range(N))
# tail of A: term ratio t_{n+1}/t_n = x/((n+1)(n+2)^q) ; for n>=N-1 ratio <= x/(N*(N+1)^q)
tN=x**N/(factorial(N)*factorial(N+1)**q); rho=x/((N+1)*(N+2)**q)
Aup=A+tN/(1-rho)
piup=F(314160,100000)
lhs=piup**12*Aup**2; rhs=R*S
print("rho",float(rho),"A",float(A),"R",float(R),"S",float(S))
print("ratio RS/(pi^12 A^2) >=", float(rhs/lhs), "certificate", rhs>lhs)
