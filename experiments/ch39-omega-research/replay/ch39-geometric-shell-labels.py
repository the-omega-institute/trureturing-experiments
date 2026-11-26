"""Singlet/doublet tetrahedral labels versus a transitive perfect-matching fiber."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq_mat

arc = run_path(str(Path(__file__).with_name("ch39-arc-step-count.py")))
q = arc["q"]


def orbit_counts():
    # A matching is an unordered pair of unordered pairs on four vertices.
    canonical = lambda pairs: tuple(sorted(tuple(sorted(pair)) for pair in pairs))
    original = ((0,1),(2,3))
    images = Counter(canonical([(p[a],p[b]) for a,b in original]) for p in permutations(range(4)))
    assert len(images) == 3 and set(images.values()) == {8}
    one = fmpq_mat([[1]*3]*3)/3
    two = fmpq_mat([[int(i==j) for i in range(3)] for j in range(3)])-one
    assert one*one == one and two*two == two and (one*two).rank() == 0
    assert (one.rank(),two.rank()) == (1,2)
    for perm in permutations(range(3)):
        P = fmpq_mat([[int(perm[j]==i) for j in range(3)] for i in range(3)])
        assert P*one*P.transpose() == one and P*two*P.transpose() == two
    shell = [p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867]
    row = {51:0,49:1,47:2,35:2}
    result = Counter()
    for p in shell:
        for values in ((1,4,3),(2,3,5),(2,3,5)):
            result[values[row[sum(p)]]] += 1
    counts = tuple(result[k] for k in range(1,6))
    assert counts == (1,2,21,3,30)
    return counts


def verify():
    counts = orbit_counts()
    weights = [q(F(n,57)) for n in counts]
    initial = [q(F(v,1000)) for v in (523,797,829,847,900)]
    for zero,(bounds,_) in arc["dyadic"]["DATA"].items():
        targets = [arb(q(F(2*b+1,2*10**20)),q(F(1,2*10**20))) for b in bounds]
        arc["upper_certificate"](initial,targets,2**53,weights)
        print(zero,"tetrahedron without orientation or with either orientation:",counts)


if __name__ == "__main__":
    verify()
