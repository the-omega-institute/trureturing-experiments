"""CP traces on every positive height and rare-event density of their joint norm."""
from fractions import Fraction as F
from itertools import combinations,permutations,product
from flint import arb,acb,fmpz_mat

def cross(x,y):
    return [x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3)]
def dot(x,y):
    return sum(a*b for a,b in zip(x,y))
def ball(x):
    return arb(x.numerator)/x.denominator
def ordered_sites(r):
    from math import isqrt
    return sorted((p for p in product(range(isqrt(r)+1),repeat=3) if dot(p,p)==r),
                  key=lambda p:(0 not in p,p))
sites=ordered_sites(17)
digits=(1,5,13,23,29)
mapped=[tuple(digits[k] for k in v) for v in sites]
normals=[cross(s[0],s[i]) for s in (sites,mapped) for i in (7,8)]
det={k:dot(normals[k[0]],cross(normals[k[1]],normals[k[2]])) for k in combinations(range(4),3)}
d2={k:F(d*d,dot(normals[k[0]],normals[k[0]])*dot(normals[k[1]],normals[k[1]])*
        dot(normals[k[2]],normals[k[2]])) for k,d in det.items()}
assert all(det.values())
mean=sum(d2.values())/8
squares=sum(c*c for c in d2.values())
second=F(3,8)**3*squares + F(3,8)**2/4*(sum(d2.values())**2-squares)
print("oriented triple determinants",det,"normalized squares",d2)
print("mean and variance of sum of squared CP traces",mean,second-mean**2)

# For a rigorous coefficient bound, integrate out the permutation ordering and the radial
# coordinate, then use monotone endpoints on the remaining three-dimensional unit cube.
def volume(n):
    o3=-arb(2)/3
    grid=[arb(i)/n for i in range(n+1)]
    sixth=[x**6 for x in grid]
    tt=[x*x for x in grid]
    total=[arb(0),arb(0)]
    for order in permutations(range(4)):
        m=[ball(d2[tuple(sorted(order[i] for i in ids))]) for ids in
           ((0,1,2),(0,1,3),(0,2,3),(1,2,3))]
        for endpoint in (0,1):
            # Integral weight y is exact on each slab. endpoint=1 lower integrand.
            for j in range(n):
                Y=sixth[j+endpoint]
                slab=arb(2*j+1)/(2*n**4)
                for i in range(endpoint,n+endpoint):
                    base=m[1]+Y*(m[2]+sixth[i]*m[3])
                    total[endpoint]+=slab*sum((m[0]+base*tt[k])**o3 for k in range(endpoint,n+endpoint))
    factor=arb(9)/4*(2/arb.pi())**4
    return (factor*total[1]).lower().union((factor*total[0]).upper())
c=volume(16)
print("leading coefficient for CP trace norm below epsilon",
      c.lower().str(12),c.upper().str(12))
assert c>0
