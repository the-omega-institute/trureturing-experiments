"""Three labelled mass thresholds, with their fourth powers in geometric progression."""

from collections import Counter
from fractions import Fraction as Q
from itertools import product


def layer_differences():
    # Twelve signed occurrences in all three layers. A nonzero pole difference
    # needs at least two contributing layers, so six bin labels cover them all.
    profiles = set()
    for a, b, c, d in product(range(6), repeat=4):
        row = [0]*6
        for location, weight in ((a, 3), (b, 1), (c, -3), (d, -1)):
            row[location] += weight
        profiles.add(tuple(row))
    return profiles


def rational_profiles(ratio, profiles):
    p, q = ratio.numerator, ratio.denominator
    counts = Counter()
    for first in profiles:
        if any(x % p for x in first):
            continue
        for last in profiles:
            if any(x % q for x in last) or not (any(first) or any(last)):
                continue
            numerators = [-q*q*x-p*p*z for x, z in zip(first, last)]
            if any(v % (p*q) for v in numerators):
                continue
            middle = tuple(v//(p*q) for v in numerators)
            if middle not in profiles:
                continue
            signed = [p*q*y+2*p*p*z for y, z in zip(middle, last)]
            counts[sum(v != 0 for v in signed)] += 1
    return counts


def witness(r, slopes, rows):
    profiles = layer_differences()
    assert len(set(slopes)) == len(slopes) and min(slopes) > 0
    assert all(tuple(row)+(0,)*(6-len(row)) in profiles for row in rows)
    assert all(sum(r**i*row[j] for i, row in enumerate(rows)) == 0 for j in range(len(slopes)))
    log_weights = [sum(i*r**i*row[j] for i, row in enumerate(rows)) for j in range(len(slopes))]
    assert all(sum(w*x**k for x, w in zip(slopes, log_weights)) == 0 for k in (0, 1, 2))


def main():
    profiles = layer_differences()
    assert len(profiles) == 961
    ratios = sorted({Q(p,q) for p in range(2, 5) for q in range(1, p)})
    for ratio in ratios:
        counts = rational_profiles(ratio, profiles)
        candidates = sum(n for k, n in counts.items() if k == 0 or k >= 4)
        assert candidates == {Q(2): 30, Q(3): 360}.get(ratio, 0)
        print(ratio, "residue patterns:", sum(counts.values()), "passing the moment support test:", candidates)

    witness(2, (1, 2), ((4, -4), (-4, 4), (1, -1)))
    witness(3, (1, 5, 3, 7), ((-3, 0, 0, 3), (1, -3, 3, -1), (0, 1, -1, 0)))
    # At the golden ratio the residue relation is 1+r-r^2=0. The
    # signed distribution common to all three layers kills all moments
    # through degree two even before multiplying by the log mass factor.
    slopes, signed = (3, 7, 5, 1), (3, 1, -3, -1)
    assert signed+(0, 0) in profiles
    assert all(sum(d*x**k for d, x in zip(signed, slopes)) == 0 for k in (0, 1, 2))


if __name__ == "__main__":
    main()
