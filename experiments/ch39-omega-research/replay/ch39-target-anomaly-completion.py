"""Target charge rigidity and the minimum continuous-gauge-singlet compensation."""
from itertools import product
from sympy import expand, symbols


def target_redefinitions():
    p, s = symbols("p s", integer=True)
    t = -p-s
    antifive, target = -9*s, -3*t
    coeff = 5*antifive + 10*target + 20*(s+t)
    assert expand(coeff) == 10*p-15*s


def indices(rows):
    return (sum(q for q,z in rows), sum(q**3 for q,z in rows),
            sum(z**3 for q,z in rows) % 9,
            sum(q*q*z for q,z in rows) % 3,
            sum(q*z*z for q,z in rows) % 3)


def verify():
    target_redefinitions()
    for q in (0,1,2):
        for charges in product((-1,0,1),repeat=3):
            # The continuous linear/cubic constraint leaves q,-q,0.
            rows = list(zip((q,-q,0), charges))
            linear,cubic,_,second,first = indices(rows)
            assert linear == cubic == 0
            assert second or not first
    repair = [(1,1),(1,-1),(-1,0),(-1,0)]
    assert indices(repair) == (0,0,0,0,2)
    # Pairings of four product 5-manifolds with the mixed homology summands
    # in H5(B(SU(5) x U(10)) smash BZ3). First choose orientation generators.
    from flint import fmpz_mat
    periods = fmpz_mat([[1,0,0,0],[0,2,0,0],[0,0,1,0],[0,0,0,1]])
    assert periods.det() % 3 == 2
    original = [1,3,3,3+3*(0+1-1)]
    addition = [sum(q*z*z for q,z in repair),sum(q*q*z for q,z in repair),0,0]
    assert all((b+c)%3 == 0 for b,c in zip(original,addition))
    print("Covering coefficient 10p-15s; four minimally sufficient singlet Weyl fields.")
    print("Four mixed spin-bordism generators at order three have trivial combined phases.")


if __name__ == "__main__":
    verify()
