"""Compare induced Gaussian weights on equal-volume, equal-total-curvature four-spheres."""
from pathlib import Path
import runpy

from flint import arb, fmpz_mat, fmpz_poly

root = Path(__file__).resolve().parent
moves = runpy.run_path(str(root / "ch39-regge-move-readout.py"))
one_five = moves["one_five"]
third = moves["third"]


def dual_polynomial(triangulation):
    facets = sorted(tuple(sorted(v)) for v in triangulation)
    n = len(facets)
    rows = [[-5 if i==j else int(len(set(v)&set(w))==4)
             for j,w in enumerate(facets)] for i,v in enumerate(facets)]
    assert all(sum(r)==0 for r in rows)
    return fmpz_mat(rows).charpoly()


def print_coeff(z):
    print([int(c) for c in z.coeffs()])


polynomials = [dual_polynomial(one_five(third,frozenset(moves[k][1]),8)) for k in ("low","high")]
gcd = polynomials[0].gcd(polynomials[1])
p, r = [v//gcd for v in polynomials]
print("gcd degree",gcd.degree(),"reduced pair")
print_coeff(p)
print_coeff(r)
print("second minus first")
print_coeff(r-p)

# Exactly bound the rational function on the whole positive mass axis.
difference = p-r
X = fmpz_poly([0,1])
assert all(a>0 for a in difference.coeffs()) and all(a>0 for a in r.coeffs())
assert difference.degree()+3==r.degree() and difference.leading_coefficient()==8
assert all(a>=0 for a in (8*r-X**3*difference).coeffs())
gauge = runpy.run_path(str(root / "ch39-quaternion-gauge.py"))
q = gauge["q"]
original_d = q(460)**2/q("1.04")
h = q(".01")
a = arb.pi()**2/original_d
G = (1-16*a*h*h/arb.pi()**2)**2-9*a*h
beta = a*h*h/6
def evaluate(poly, x):
    value = arb(0)
    for coefficient in reversed(poly.coeffs()):
        value = value*x+int(coefficient)
    return value

def one_mode(x):
    return (evaluate(difference,x)/evaluate(r,x)).log1p()/2

low_modes = arb(0)
for box in gauge["xs"]:
    lo = (box.lower()/(1+beta*(box.lower()+10))).lower()
    hi = (box.upper()/G).upper()
    cuts = [lo+(hi-lo)*j/128 for j in range(129)]
    pieces = [one_mode(x.union(y)) for x,y in zip(cuts,cuts[1:])]
    low_modes += min(v.lower() for v in pieces).union(max(v.upper() for v in pieces))
# The 19th continuous energy exceeds 10^4 by the high-site Schur test,
# and the boundary-zero subspace has codimension 38.
k = 100*arb.pi()/original_d.sqrt()
assert 959990-10000-36*k/k.sin()>0
assert 10000/(1+beta*10010)>9000
tail = 4*(q(20)/9000**3+270*(3/original_d)**3*q(6)/5)
# sum_{k>=1} 1/k^6 <= 1+integral_1^infty t^{-6}dt = 6/5.
full = low_modes+arb(0).union(tail.upper())
assert full>0
print("full log weight ratio, high q^2 over low",full.lower().str(12),full.upper().str(12))
print("post-18 absolute tail bound",tail.upper().str(12))
