#!/usr/bin/env python3
"""Exact checks of Gerber--Norton, arXiv:2609.03738v1, Table 2.

Partitions are tuples of positive decreasing parts. A beta set is represented
by all its elements >= LOW, with every integer < LOW implicitly present.
All operations are integer operations. This is finite evidence, not a proof
of universal positivity or of the three-square theorem.
"""
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from math import isqrt

LOW = -96
CHARGE = (0, 0, 1, 2)


@lru_cache(None)
def partitions(n, bound=None):
    if n == 0:
        return ((),)
    if bound is None:
        bound = n
    return tuple((first,) + tail for first in range(min(n, bound), 0, -1)
                 for tail in partitions(n - first, first))


def beta(part, charge):
    # k is the source's positive row index. After the last part, beta=a-k+1.
    stop = max(len(part), charge - LOW + 1)
    return frozenset((part[k - 1] if k <= len(part) else 0) - k + 1 + charge
                     for k in range(1, stop + 1)
                     if (part[k - 1] if k <= len(part) else 0) - k + 1 + charge >= LOW)


def nested(parts, charges, e=4):
    rows = [beta(p, a) for p, a in zip(parts, charges)]
    return all(rows[j] <= rows[j + 1] for j in range(len(rows) - 1)) and (
        rows[-1] <= beta(parts[0], charges[0] + e))


def multi_of_size(n, ell=4):
    if ell == 1:
        yield from ((p,) for p in partitions(n))
    else:
        for a in range(n + 1):
            for p in partitions(a):
                for rest in multi_of_size(n - a, ell - 1):
                    yield (p,) + rest


def direct_counts(charges, limit):
    return [sum(nested(parts, charges) for parts in multi_of_size(n, len(charges)))
            for n in range(limit + 1)]


def reconstruct(beads):
    charge = sum(z > 0 for z in beads) - sum(z not in beads for z in range(LOW, 1))
    seq = sorted(beads, reverse=True)
    parts = [z + k - 1 - charge for k, z in enumerate(seq, 1)]
    assert parts[-1] == 0, (charge, parts[-5:])
    assert all(x >= 0 for x in parts)
    assert all(a >= b for a, b in zip(parts, parts[1:]))
    part = tuple(x for x in parts if x)
    assert beta(part, charge) == beads
    signed = sum(z for z in beads if z > 0) - sum(z for z in range(LOW, 1) if z not in beads)
    assert sum(part) == signed - charge * (charge + 1) // 2
    return part, charge


def bead_construction(t, r, s):
    assert sum(t) == 0 and 1 <= r <= 4 and 1 <= s <= 4 and r != s
    assert min(t) > LOW // 4 + 3
    b = tuple(i + 4 * t[i - 1] for i in range(1, 5))
    B = frozenset(z for z in range(LOW, max(4, max(b)) + 1)
                  if (z - ((z - 1) % 4 + 1)) // 4 < t[(z - 1) % 4])
    assert all(x not in B for x in b)
    B4 = frozenset(z for z in range(LOW, max(8, max(b) + 4) + 1)
                   if z - 4 < LOW or z - 4 in B)
    assert B4 == B | set(b)
    rows = (B, B, B | {b[r - 1]}, B | {b[r - 1], b[s - 1]})
    pc = [reconstruct(row) for row in rows]
    parts, charges = tuple(p for p, _ in pc), tuple(a for _, a in pc)
    assert charges == CHARGE and nested(parts, charges)
    M = 2 * sum(x * x for x in t) + sum(i * t[i - 1] for i in range(1, 5))
    assert sum(parts[0]) == M and parts[0] == parts[1]
    assert sum(parts[2]) == M + b[r - 1] - 1
    assert sum(parts[3]) == M + b[r - 1] + b[s - 1] - 3
    N = 4 * M + 2 * b[r - 1] + b[s - 1] - 4
    assert N == sum(map(sum, parts))
    c = tuple(i - 3 + 2 * (i == r) + (i == s) for i in range(1, 5))
    y = tuple(4 * a + z for a, z in zip(t, c))
    assert sum(c) == 1 and sum(z * z for z in c) == 4 * r + 2 * s - 7
    assert sum(y) == 1 and sum(z * z for z in y) == 2 * N + 1
    return parts, N


RESIDUE_PAIRS = ((1, 2), (1, 4), (1, 3), (2, 4), (2, 3), (4, 1), (3, 1), (3, 4))


def residue_check():
    possible = {a for a in combinations_with_replacement(range(4), 4) if sum(a) % 4 == 1}
    realized = set()
    for r, s in RESIDUE_PAIRS:
        c = tuple(i - 3 + 2 * (i == r) + (i == s) for i in range(1, 5))
        realized.add(tuple(sorted(x % 4 for x in c)))
    assert possible == realized
    print('Residue multisets:', [''.join(map(str, a)) for a in sorted(possible)])


def three_odd_squares(n):
    target = 8 * n + 3
    for a in range(1, isqrt(target) + 1, 2):
        for b in range(1, isqrt(target - a * a) + 1, 2):
            rem = target - a * a - b * b
            c = isqrt(rem)
            if c % 2 == 1 and c * c == rem:
                return tuple(x if x % 4 == 1 else -x for x in (a, b, c))
    raise AssertionError(('no three odd squares in finite test', n))


def from_n(n):
    a, b, c = three_odd_squares(n)
    numerators = (1 + a + b + c, 1 + a - b - c, 1 - a + b - c, 1 - a - b + c)
    assert all(z % 4 == 0 for z in numerators)
    y = tuple(z // 4 for z in numerators)
    assert sum(y) == 1 and sum(z * z for z in y) == 2 * n + 1
    for r, s in RESIDUE_PAIRS:
        cv = tuple(i - 3 + 2 * (i == r) + (i == s) for i in range(1, 5))
        for yp in permutations(y):
            if all((z - v) % 4 == 0 for z, v in zip(yp, cv)):
                t = tuple((z - v) // 4 for z, v in zip(yp, cv))
                parts, N = bead_construction(t, r, s)
                assert N == n
                return parts
    raise AssertionError(('unmatched residues', n, y))


def main():
    # Section 1.1: the source's 3-partition of total size 16.
    assert sum(map(sum, ((3, 2, 1), (4, 1), (1, 1, 1, 1, 1)))) == 16
    # Example immediately following Lemma 2.11 (the original-core criterion).
    assert nested(((1, 1), (3, 1, 1)), (1, 2))
    print('Fidelity: size-16 example and Lemma 2.11 bipartition example OK')
    # Table 2: equal charge cores are repeated ordinary 4-cores; the displayed
    # series is psi(q^4)*psi(q^8)^2. Count its coefficients independently.
    equal = direct_counts((0, 0, 0, 0), 12)
    series = [0] * 13
    for a, b, c in product(range(4), repeat=3):
        degree = 4 * (a * (a + 1) // 2) + 8 * (b * (b + 1) // 2 + c * (c + 1) // 2)
        if degree <= 12:
            series[degree] += 1
    assert equal == series
    print('Fidelity: Table 2 charge 0000 counts n=0..12:', equal)
    other = direct_counts((0, 0, 0, 1), 8)
    assert other[4] == 0  # the first named zero in the source table
    print('Fidelity: Table 2 charge 0001 counts n=0..8:', other)
    target = direct_counts(CHARGE, 12)
    assert all(x > 0 for x in target)
    print('Direct target counts n=0..12:', target)
    # Equation (5.5) is not an input to the construction. Its printed constant
    # coefficient is inconsistent with the source's counting definition.
    alt0 = direct_counts((0, 1, 1, 1), 0)[0]
    assert target[0] == alt0 == 1
    print('Equation (5.5) constant-term discrepancy: LHS=2, RHS=1; not used')
    residue_check()
    checked = 0
    for a, b, c in product(range(-3, 4), repeat=3):
        d = -a - b - c
        if -3 <= d <= 3:
            for r, s in permutations(range(1, 5), 2):
                bead_construction((a, b, c, d), r, s)
                checked += 1
    print('Signed beta-set/charge/size/identity cases:', checked)
    for n in range(513):
        parts = from_n(n)
        if n == 0:
            assert parts == ((), (), (), ())
        if n <= 12:
            assert parts in tuple(multi_of_size(n))
    print('Three-square-to-partition construction: every n=0..512 OK, including n=0')
    assert bead_construction((0, 0, 0, 0), 1, 2) == (((), (), (), ()), 0)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
