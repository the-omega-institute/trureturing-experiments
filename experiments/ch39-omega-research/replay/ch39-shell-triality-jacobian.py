"""Baryon triality phases of the finite magnetic overlap index table."""
from collections import Counter
from itertools import combinations, product


def index(charge):
    return (0, 50, 25, 25, 50)[charge % 5]


def sector_indices(charges):
    color, weak = charges[:3], charges[3:]
    return (sum(index(-c) for c in color), sum(index(-w) for w in weak),
            sum(index(a+b) for a,b in combinations(color,2)),
            sum(index(c+w) for c in color for w in weak), index(sum(weak)))


def exponent(indices):
    down, lepton, up, quark, electron = indices
    return down - lepton - up + 2*electron


def verify():
    flux = (1, 1, 1, -3, 0)
    assert sector_indices(flux) == (150, 25, 75, 225, 25)
    assert exponent(sector_indices(flux)) == 100
    counts = Counter()
    with_fixed_weak = Counter()
    for prefix in product(range(5), repeat=4):
        q = prefix+((-sum(prefix)) % 5,)
        phase = exponent(sector_indices(q)) % 3
        counts[phase] += 1
        if 0 in q[3:]:
            with_fixed_weak[phase] += 1
            if 0 in q[:3] or q[3:] == (0, 0):
                assert phase == 0
    assert counts == {0:289, 1:168, 2:168}
    assert with_fixed_weak == {0:225-56, 1:28, 2:28}
    # For one color plus one weak parallel vector the residual group is U(2).
    # q (fundamental) and d (determinant character) denote their independent indices.
    for q,d in product((-3,0,7), repeat=2):
        assert (q-d-(q+d)+2*d) == 0
    print("Magnetic triality exponents", counts, "and one weak fixed direction", with_fixed_weak)


if __name__ == "__main__":
    verify()
