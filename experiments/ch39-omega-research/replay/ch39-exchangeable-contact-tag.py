"""Use the same rotor records to extract the one-bypass, three-comparison label."""
from fractions import Fraction
from itertools import permutations, product
from math import prod
from pathlib import Path
from runpy import run_path

coupling = run_path(str(Path(__file__).with_name("ch39-schur-winner-gauge-yukawa.py")))
perms = list(permutations(range(4)))
axes = [(0,*p) for p in permutations(range(1,4))]


def pointer(row):
    if row.count(0)==1:
        return row.index(0)
    single = [x for x in (1,2,3) if row.count(x)==1]
    return row.index(single[0]) if len(single)==1 else None


def verify():
    rotor = (1,0,0,0)
    for k in coupling["word"]:
        rotor = coupling["words"]["multiply"](rotor,coupling["words"]["GENERATORS"][int(k)])
    counts = [v*v for v in rotor]
    assert sum(counts) == 5**128
    tags = [0]*4
    marked = 0
    for row in product(range(4), repeat=4):
        tag = pointer(row)
        fixed = set(range(4))
        for labels in axes:
            for perm in perms:
                transformed = tuple(labels[row[i]] for i in perm)
                assert pointer(transformed) == (None if tag is None else perm.index(tag))
                if transformed == row:
                    fixed &= {i for i in range(4) if perm[i]==i}
        assert bool(fixed) == (tag is not None)
        if tag is not None:
            marked += 1
            tags[tag] += prod(counts[i] for i in row)
    assert marked == 144 and len(set(tags)) == 1
    success = Fraction(sum(tags),sum(counts)**4)
    assert Fraction(5666906,10**7) < success < Fraction(5666907,10**7)
    assert (1-success)**20 < Fraction(1,2**24)
    print("four equiprobable tags; optimal equivariant success",float(success),
          "20 independent blocks miss with probability below 2^-24")


if __name__ == "__main__":
    verify()
