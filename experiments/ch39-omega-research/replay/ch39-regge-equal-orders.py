"""Find two explicit four-spheres with identical hinge orders and different scalar weights."""
from collections import Counter
from pathlib import Path
import runpy

from flint import fmpz_mat

root = Path(__file__).resolve().parent
moves = runpy.run_path(str(root/"ch39-regge-move-readout.py"))
first = moves["third"]


def polynomial(triangulation):
    sets = sorted(triangulation, key=lambda v: tuple(sorted(v)))
    return fmpz_mat([[-5 if i==j else int(len(v&w)==4)
                      for j,w in enumerate(sets)] for i,v in enumerate(sets)]).charpoly()


pair = [
    ((0,1,2,4,6),(0,1,2,5,7)),
    ((0,1,2,4,7),(0,1,2,5,7)),
]
triangulations = []
for left,right in pair:
    one = moves["one_five"](first,frozenset(left),8)
    triangulations.append(moves["one_five"](one,frozenset(right),9))
signature = {3:33,4:21,5:3,6:7}
assert all(len(T)==24 and Counter(moves["hinges"](T).values())==signature for T in triangulations)
A,B = [polynomial(T) for T in triangulations]
gcd = A.gcd(B)
p,r = A//gcd,B//gcd
difference = p-r
from flint import fmpz_poly, arb
X=fmpz_poly([0,1])
assert all(a>0 for a in difference.coeffs()) and all(a>0 for a in r.coeffs())
assert r.degree()==22 and difference.degree()==16 and difference.leading_coefficient()==42
assert all(z>=0 for z in (42*r-X**6*difference).coeffs())
# Reuse the same certified squared masses and the finite-transfer bounds.
sibling=runpy.run_path(str(root/"ch39-regge-shell-determinant.py"))
gauge=sibling["gauge"]
evaluate=sibling["evaluate"]
beta,G=sibling["beta"],sibling["G"]
q=gauge["q"]
total=arb(0)
for v in gauge["xs"]:
    lo=(v.lower()/(1+beta*(v.lower()+10))).lower()
    hi=(v.upper()/G).upper()
    cuts=[lo+(hi-lo)*j/128 for j in range(129)]
    terms=[]
    for a,b in zip(cuts,cuts[1:]):
        x=a.union(b)
        terms.append((evaluate(difference,x)/evaluate(r,x)).log1p()/2)
    total+=min(v.lower() for v in terms).union(max(v.upper() for v in terms))
tail=21*(q(20)/9000**6+270*(3/sibling["original_d"])**6*q(12)/11)
total+=arb(0).union(tail.upper())
print("equal-orders full log weight ratio",total.lower().str(14),total.upper().str(14))
print("post-18 tail",tail.upper().str(13),"reduced difference coefficients",list(map(int,difference.coeffs())))
