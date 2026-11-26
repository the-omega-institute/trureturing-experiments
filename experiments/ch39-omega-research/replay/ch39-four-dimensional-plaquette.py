"""Different rectangular Wilson expectations on two four-cubes with identical bare weak jets."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path
import runpy

weights=runpy.run_path(str(Path(__file__).with_name('ch39-gauge-jet-counterexample.py')))
extent=(2,1,1,1)
def move(point,axis):
    return tuple(v+int(i==axis) for i,v in enumerate(point))
def edge(a,b):
    return tuple(sorted((a,b)))
def plaquette(p,i,j):
    a,b=move(p,i),move(p,j)
    c=move(a,j)
    return frozenset((edge(p,a),edge(a,c),edge(c,b),edge(b,p)))
faces=[plaquette(p,i,j) for p in product(*(range(n+1) for n in extent))
       for i,j in combinations(range(4),2) if p[i]<extent[i] and p[j]<extent[j]]
p=(0,0,0,0)
left,right=plaquette(p,0,1),plaquette(move(p,0),0,1)
boundary=left.symmetric_difference(right)
assert len(faces)==42 and len(boundary)==6
assert [set(pair) for pair in combinations(faces,2) if boundary<=pair[0]|pair[1]]==[{left,right}]
s=F(1,10**10)
tail=sum(comb(42,k)*(2*s)**k for k in range(3,43))
bounds=[]
for w in (weights['a'],weights['b']):
    # Haar convolution on the two tiles followed by the boundary character integral
    # contributes exactly (p_1/3^2)^2 in front of s^2.
    leading=s*s*w[1]**2/81
    assert leading>tail
    bounds.append((leading-tail,(leading+tail)/(1-2*s)**42))
gap=bounds[0][0]-bounds[1][1]
assert gap>0
print('rectangle normalized adjoint trace intervals',[(float(l),float(r)) for l,r in bounds])
print('positive expectation gap',float(gap))

# Condition on every link outside the two tiled plaquettes. On the infinite hypercubic
# lattice, only these incident plaquette weights remain in the local conditional density.
local_links=left|right
incident={plaquette(p,i,j) for p in product(range(-1,3),repeat=4)
          for i,j in combinations(range(4),2) if local_links&plaquette(p,i,j)}
assert len(incident)<=7*6
assert [set(pair) for pair in combinations(incident,2) if boundary<=pair[0]|pair[1]]==[{left,right}]
print('incident plaquettes in the full four-dimensional lattice',len(incident))
