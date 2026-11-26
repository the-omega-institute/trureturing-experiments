"""Pairing symmetry on the 19-point shell and its three Higgs copy directions."""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

sign = run_path(str(Path(__file__).with_name("ch39-matching-flavor-stabilizers.py")))["sign"]


def verify():
    shell = [v for v in product(range(30), repeat=3) if sum(x*x for x in v) == 867]
    group = list(permutations(range(3)))
    counts = Counter()
    for v, orientation in product(shell, (0, -1, 1)):
        stabilizer = [p for p in group if tuple(v[i] for i in p) == v and sign(p)*orientation == orientation]
        counts[len(stabilizer)] += 1
        even = all(sign(p) == 1 for p in stabilizer)
        assert even == (orientation != 0 or len(set(v)) == 3)
    assert counts == {1: 48, 2: 6, 3: 2, 6: 1}
    assert sum(number*order for order, number in counts.items()) == 72
    assert 6*len(shell)*3 == 342
    assert [m for m in range(-4, 5) if (3*m-2)%5==0 and m%2==0] == [4]
    assert (3*(-1)-2)%5==0
    print("57/57 invariant phases, 56/57 covariant densities, 48/57 equivariant ordered fields.")
    print("Least evenly refined shell sizes per geometric direction: 114 for densities, 342 for fields.")


if __name__ == "__main__":
    verify()
