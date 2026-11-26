"""Minimum residual Abelian index, symmetric mass pairings, and a one-loop threshold."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from flint import arb, ctx

original = ((1,-1),(1,-1),(-1,-1),(-1,0))


def three_tests(fields):
    return tuple(sum(n for n,q in fields if (k*n+q)%3) for k in range(3))


def verify():
    possibilities = list(product((-1,1), range(3)))
    optimal = []
    for fields in combinations_with_replacement(possibilities, 4):
        if three_tests(original+fields) == (0,0,0):
            optimal.append(Counter(fields))
    mirrors = Counter((-n,(-q)%3) for n,q in original)
    assert optimal == [mirrors]
    for ranks in product(range(3), range(2), range(2)):
        massive_pairs = sum(ranks)
        assert sum(n-r for n,r in zip((2,1,1),ranks)) == 4-massive_pairs
    fundamental_charges = (1,1,0,-2,0)
    norm2 = 2*sum(n*n for n in fundamental_charges)
    assert norm2 == 12
    integer_beta = F(2,3)*4*2**2
    assert integer_beta == F(32,3) and integer_beta/norm2 == F(8,9)
    assert 3*integer_beta/norm2 == F(8,3)
    ctx.prec = 120
    decade = 4*arb(10).log()/(9*arb.pi())
    assert arb("0.32574") < decade < arb("0.32576")
    print("Unique minimum charge spectrum and one-decade normalized inverse-coupling shift",decade.str(12))


if __name__ == "__main__":
    verify()
