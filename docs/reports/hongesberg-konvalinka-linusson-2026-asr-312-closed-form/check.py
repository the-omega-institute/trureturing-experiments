# Höngesberg–Konvalinka–Linusson arXiv:2610.07442 Conjecture 4.5: recurrences (Prop 4.2, 4.3, Thm 4.4) vs the closed form.
from functools import lru_cache
from math import comb
def C(m,j):
    # generalized binomial (m any integer, j >= 0); 0 for j < 0
    if j<0: return 0
    num=1
    for t in range(j): num*= (m-t)
    from math import factorial
    return num//factorial(j)
def Sch(n):  # large Schroeder numbers, Sch_{-1}=1
    if n<0: return 1
    s=[1]
    for m in range(1,n+1):
        s.append(s[m-1]+sum(s[k]*s[m-1-k] for k in range(m)))
    return s[n]
def Cat(n): return comb(2*n,n)//(n+1)
@lru_cache(None)
def S(r,k,d):
    if d<0 or r<d or k<d: return 0
    if d==0: return 1
    if r==k==d: return Sch(d-1)
    if d==r and r<=k:  # Prop 4.2 (1 <= r <= k), here r<k
        return S(r,k-1,r)+S(r-1,k-1,r-1)+sum(Sch(i-2)*S(r+1-i,k+1-i,r+1-i) for i in range(2,r+1))
    if d==k and k<r:   # Prop 4.3
        return sum(S(j,k,j)*comb(r-j-1,k-j) for j in range(0,k+1))
    # Thm 4.4, d < r, k
    return S(r,k-1,d)+S(r-1,k-1,d-1)+sum(Sch(i-2)*(S(r+1-i,k+1-i,d+1-i)-C(r-i,d+1-i)) for i in range(2,d+1))
def F(r,k,d):
    t1=sum(2**(i-1)*comb(k,i)*(C(r-2,d-i)-2*C(r-2,d-i-2)) for i in range(1,d+1))
    t3=sum(Cat(2*i)*C(r-1+2*i,d-1-2*i) for i in range(0,d+1))
    t4=2*sum((-1)**i*C(r,d-i) for i in range(1,d+1))
    return t1+C(r,d)+t3+t4
bad=[];n=0
for r in range(0,20):
    for k in range(0,20):
        for d in range(0,min(r,k)+1):
            n+=1
            if S(r,k,d)!=F(r,k,d): bad.append((r,k,d,S(r,k,d),F(r,k,d)))
print("triples",n,"mismatches",len(bad),bad[:8])
print("S(5,6,4) =",S(5,6,4)," sample S(3,3,3)",S(3,3,3),"Sch(2)",Sch(2))
