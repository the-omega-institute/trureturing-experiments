"""Check left quotients of short rotations and recover the radial assignment."""
from collections import Counter
from itertools import combinations_with_replacement, product
from math import gcd
from pathlib import Path
from runpy import run_path


def primitive_line(q):
    divisor = gcd(*q)
    if next(x for x in q if x) < 0:
        divisor = -divisor
    return tuple(x // divisor for x in q)


def quotient_counts():
    words = run_path(str(Path(__file__).with_name("ch39-hurwitz-gauge-words.py")))
    endpoints = set()
    counts = []
    for depth in range(18):
        hits = words["hits"](depth)
        counts.append(len(hits))
        for (a, b, c, d), _ in hits:
            for i, j, k in product((-1, 1), repeat=3):
                signed = (a, i*b, j*c, k*d)
                # Check that each sign choice is actually reached at this depth.
                words["factor"](signed, depth)
                endpoints.add(primitive_line(signed))
    assert counts == [0]*16 + [4, 39]
    assert len(endpoints) == 344
    multiplicities = Counter()
    for a, b, c, d in endpoints:
        conjugate = (a, -b, -c, -d)
        for other in endpoints:
            quotient = primitive_line(words["multiply"](other, conjugate))
            multiplicities[quotient] += 1
    assert multiplicities.pop((1, 0, 0, 0)) == 344
    assert len(multiplicities) == 117992 and set(multiplicities.values()) == {1}
    print("344 rotations, 117992 distinct nonidentity ordered quotients")


def recover_assignment():
    sizes = (1, 9, 15)
    pairs = tuple(combinations_with_replacement(range(3), 2))
    counts = [sizes[a]*sizes[b]*(1 if a == b else 2) for a, b in pairs]
    observed = (1, 18, 81, 255, 270)
    solutions = []
    for labels in product(range(len(observed)), repeat=len(pairs)):
        totals = [0]*len(observed)
        for label, count in zip(labels, counts, strict=True):
            totals[label] += count
        if tuple(totals) == observed:
            solutions.append(tuple(observed[label] for label in labels))
    assert solutions == [(1, 18, 255, 81, 270, 255)]
    print("unique assignment by observed multiplicity:", solutions[0])


if __name__ == "__main__":
    quotient_counts()
    recover_assignment()
