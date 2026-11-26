"""Finite partial fractions for every winding sector of the cubic shell return."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_poly
from sympy import Poly, fraction

ROOT = Path(__file__).parent
GRAPH = run_path(str(ROOT/"ch39-weighted-arc-zero-phase.py"))


def evaluate(polynomial, z):
    val = acb(0)
    for c in polynomial.all_coeffs():
        val = val*z+int(c)
    return val


def modes():
    u = GRAPH["U"]
    numerator,denominator = (Poly(f,u) for f in fraction(GRAPH["response"](3)))
    assert numerator.degree()==denominator.degree()==6
    integer = fmpz_poly([int(c) for c in reversed(denominator.all_coeffs())])
    assert integer.gcd(integer.derivative()).degree()==0
    blocks = []
    for root, multiplicity in integer.complex_roots():
        assert multiplicity==1
        delta = (root**2-1).sqrt()
        outer = [v for v in (root+delta,root-delta) if abs(v)>1]
        assert len(outer)==1
        tau = 1/outer[0]
        residue = evaluate(numerator,root)/evaluate(denominator.diff(),root)
        weight = -2*residue*tau/(1-tau**2)
        blocks.append((tau,weight))
    return blocks,arb(int(numerator.LC()))/int(denominator.LC())


def verify():
    ctx.prec = 220
    blocks,constant = modes()
    def coeff(m):
        z = sum(w*t**m for t,w in blocks)
        assert z.imag.contains(0)
        return z.real
    def total(m):
        z = sum(w*t**m/(1-t) for t,w in blocks)
        assert z.imag.contains(0)
        return z.real
    b,c = (arb(3)/13)**3,(arb(3)/14)**3
    ustar = (((361-2*b*b)*(361-2*c*c)).sqrt()-19)/(2*b*c)
    outer = ustar+(ustar**2-1).sqrt()
    leading = [(t,w) for t,w in blocks if (outer*t-1).contains(0)]
    assert len(leading)==1
    tstar,wstar = leading[0]
    assert wstar.real>0 and all(abs(t)<tstar.real for t,w in blocks if not (outer*t-1).contains(0))
    # Check independently against the exact value at unit edge phase.
    u = GRAPH["U"]
    exact = GRAPH["response"](3).subs(u,1)
    value = arb(int(exact.p))/int(exact.q)
    assert (constant+coeff(0)+2*total(1)-value).contains(0)
    print("nearest exterior pole",outer.str(12),"leading residue weight",wstar.real.str(10))
    for m in (131,132):
        left = 2*sum(coeff(j) for j in range(m,m+3))
        right = 2*total(m)
        assert left>0 and right-left>0 and (right-left)/left < arb("1e-18")
        print(m,"single oriented coefficient",coeff(m).str(10),
              "compatible absolute mass",left.union(right).str(11),
              "fraction at trivial phase",(left.union(right)/value).str(9))


if __name__=="__main__":
    verify()
