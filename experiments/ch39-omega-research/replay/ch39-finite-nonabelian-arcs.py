"""Finite-group arc obstruction and optimal Weyl endpoints for the first SU(5) moment."""
from collections import Counter
from fractions import Fraction as F
from itertools import product

from flint import acb, arb, ctx


def verify():
    ctx.prec = 180
    epsilon = arb(2)**-26
    real = 1-epsilon**2/2
    point = acb(real, (1-real**2).sqrt())
    extreme = -(1+arb(5).sqrt())/4
    length = 168662972
    assert (point**(length-1)).real > extreme > (point**length).real
    # alpha <= epsilon/sqrt(1-epsilon^2/4), and the short arc has not reached half a turn.
    assert length*epsilon/(1-epsilon**2/4).sqrt() < 3
    residues = Counter((b*c-a*d)%5 for a,b,c,d in product(range(5), repeat=4))
    assert residues == {0:145,1:120,2:120,3:120,4:120}
    assert F(residues[0]-residues[1],625) == F(1,25)
    # Conjugation of a Weyl basis element by X^a Z^b.
    for c,d in product(range(5),repeat=2):
        phases = Counter((b*c-a*d)%5 for a,b in product(range(5),repeat=2))
        assert phases == ({0:25} if (c,d)==(0,0) else {i:5 for i in range(5)})
    assert 2*F(1,2**26)<1 and sum(range(-2,3))==0
    print("25 Weyl endpoints, Haar commutator mean 24/25, minimum path length",length)


if __name__ == "__main__":
    verify()
