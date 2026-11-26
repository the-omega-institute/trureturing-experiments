"""Four diagonal prime valuations and a simultaneous degree bound for a two-shell lock."""
from itertools import product,combinations
from fractions import Fraction as F
from pathlib import Path
from math import isqrt,prod
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-algebraic-lock.py")))
bmin=m["lo"].denominator+m["hi"].denominator
cosines=[F(14,17),F(11,17),F(10,17),F(745,867),F(505,867),F(155,289)]
rad=[r.numerator**2-r.denominator**2 for r in cosines]
def rational_square(n):
    return n>=0 and isqrt(n)**2==n
for size in range(7):
    for subset in combinations(rad,size):
        n=prod(subset)
        assert not rational_square(-n) and not rational_square(-3*n)
        if size:
            assert not rational_square(n)
valuations=[]
for r in cosines:
    a,b=r.as_integer_ratio()
    assert a%17 and (a*a-b*b)%17==a*a%17
    v=0
    while b%17==0:
        b//=17
        v+=1
    valuations.append(v)
assert valuations==[1,1,1,2,2,2]
omega=((1,0),(2,0),(4,3),(5,3))
rows=[]
for i,j in omega:
    signs=[-1]*6
    signs[i]=1
    v=[a*b for a,b in zip(signs,valuations)]
    rows.append([2*(v[k]-v[l]) for k,l in omega])
assert rows==[[4,0,0,0],[0,4,0,0],[0,0,8,0],[0,0,0,8]]
# mu(K) has order two. For b even use its entire order; this is a bound valid for every b.
loss=prod([8,8,16,16])
assert loss==16384
for field in (1,64,2**20):
    denom=16*loss*field
    # Smallest integer d such that 16*loss*field*d^4 >= bmin^4.
    ceil=(bmin**4+denom-1)//denom
    root=isqrt(isqrt(ceil))
    lower=root+int(root**4<ceil)
    assert denom*(lower-1)**4<bmin**4<=denom*lower**4
    print("field degree, simultaneous lower bound on total Fourier degree",field,lower)
print("joint coefficient field degree needed when degree <=5",
      (bmin**4+16*loss*5**4-1)//(16*loss*5**4))
