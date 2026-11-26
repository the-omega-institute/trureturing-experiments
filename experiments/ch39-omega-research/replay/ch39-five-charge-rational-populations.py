"""有限几何边字在电荷五、等权臂计数与轮换保护之间的算术约束。"""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt

radius = 867
pole = (17,17,17)
outer = sorted(set(permutations((1,5,29))))
inner = sorted(set(permutations((13,13,23))))
sites = outer+inner+sorted(p for p in product(range(30),repeat=3)
                         if sum(t*t for t in p)==radius and p not in outer+inner)


def dot(x,y):
    return sum(i*j for i,j in zip(x,y))


def verify():
    assert radius==3*17**2
    count = 0
    for x,y in combinations(sites,2):
        if 2*(radius-dot(x,y))>804:
            continue
        A = radius*(radius+dot(pole,x)+dot(x,y)+dot(y,pole))
        det = sum(pole[i]*(x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3])
                  for i in range(3))
        B = -17*det
        denom = A*A+3*B*B
        assert denom>0
        re,im = Fraction(A*A-3*B*B,denom),Fraction(2*A*B,denom)
        assert re*re+3*im*im==1
        count += 1
    assert count==135
    # If the loop is separately chosen for the baryon, its fifth power still
    # has to have rational real part -a/b. The flavor only imposes this window.
    lo,hi=Fraction(2014,1000)**2,Fraction(2072,1000)**2
    candidates = []
    for q in range(1,96):
        for p in range(1,q):
            delta=q*q-p*p
            if gcd(p,q)!=1 or 3*isqrt(delta//3)**2!=delta:
                continue
            x=Fraction(p,q)
            cosine=abs(16*x**5-20*x**3+5*x)
            if cosine==0:
                continue
            a,b=cosine.numerator,cosine.denominator
            assert b==q**5//(16 if q%2==0 else 1)
            if lo < Fraction(b,a) < hi:
                candidates.append((a+b,a,b,p,q))
    result=(365689129,71672243,294016886,83,86)
    assert min(candidates)==result and 96**5//16>result[0]
    assert 86**2-83**2==3*13**2
    assert 21**2-5**2==16*26 and (21**2-5**2)%3!=0
    print("135 squared edge units in Q(sqrt(-3)); 5:21 excluded.")
    print("Minimum population in the quadratic-field relaxation",result)


if __name__=="__main__":
    verify()
