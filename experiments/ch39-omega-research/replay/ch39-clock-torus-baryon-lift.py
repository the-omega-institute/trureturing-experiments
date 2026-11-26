"""A residual ninth-root torus element has a sharp three-copy chiral-character lift."""
from collections import defaultdict
from fractions import Fraction as F

def add(*terms):
    data=defaultdict(int)
    for scale,pol in terms:
        for k,v in pol.items():
            data[k]+=scale*v
    return {k:v for k,v in data.items() if v}
def product(p,q):
    result=defaultdict(int)
    for (n,h),a in p.items():
        for (m,k),b in q.items():
            result[n+m,(h+k)%9]+=a*b
    return {k:v for k,v in result.items() if v}
one={(0,0):1}
z={(12,3):1}
x={(36,0):1}
fund={(4,1):2,(-8,7):1}
anti={(-y,(-h)%9):v for (y,h),v in fund.items()}
K=add((1,product(fund,anti)),(-5,one))
psi={(3*y,(3*h)%9):v for (y,h),v in fund.items()}
assert product(z,product(z,z))==x
first=product(x,add((1,product(K,K)),(-8,one),(-2,product(x,K))))
second=product(product(x,add((2,x),(-1,one))),psi)
assert add((1,first),(1,second))=={(12,3):3}
# Discrete scalar, left Weyl matter in standard conventions, and Delta B=1.
qnums={"Q":(F(1,3),1),"uc":(-F(1,3),-4),"dc":(-F(1,3),2),"L":(0,-3),"ec":(0,6),"nc":(0,0),"H":(0,3)}
trial={k:int(3*b+2*y)%9 for k,(b,y) in qnums.items()}
assert trial=={"Q":3,"uc":0,"dc":3,"L":3,"ec":3,"nc":0,"H":6}
assert (3*trial["Q"]+trial["L"])%9==3
assert all(n%9==0 for n in (trial["Q"]+trial["uc"]+trial["H"],
    trial["Q"]+trial["dc"]-trial["H"],trial["L"]+trial["ec"]-trial["H"],
    trial["L"]+trial["nc"]+trial["H"]))
# A chiral monomial with residual H=3 has nonzero first derivative mod 9
# unless its signed multiplicity is divisible by three.
assert [n for n in range(1,10) if 3*n%9==0]==[3,6,9]
print("three-copy isolated Laurent character",add((1,first),(1,second)))
print("Z9 weights and QQQL weight",trial,3)
