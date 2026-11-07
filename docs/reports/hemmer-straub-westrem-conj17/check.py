#!/usr/bin/env python3
"""Exact finite checks for Hemmer--Straub--Westrem, journal Conjecture 17.

Source: arXiv:2509.02796v1, Sections 1, 2.1 (Table 1), and 6 (6.1).
Only integers and fractions.Fraction participate in the calculations.
This is a finite sanity/fidelity check, not a proof for arbitrary n.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import permutations
from math import comb, factorial, prod


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def partitions(n, ceiling=None):
    """All partitions, exactly once, as decreasing tuples of positive parts."""
    if n == 0:
        yield ()
        return
    for a in range(min(n, n if ceiling is None else ceiling), 0, -1):
        for tail in partitions(n - a, a):
            yield (a,) + tail


def ev_types(lam):
    """Retain every binary replacement, including repeated resulting types."""
    types = [()]
    for a in lam:
        types = [tuple(sorted(t + replacement, reverse=True))
                 for t in types for replacement in [(2 * a,), (a, a)]]
    return Counter(types)


@cache
def r3(n):
    return tuple(mu for mu in partitions(n)
                 if len(mu) <= 3 and all(a % 2 == 0 for a in mu))


def z_factor(lam):
    return prod(a ** m * factorial(m) for a, m in Counter(lam).items())


@cache
def power_sum_coefficients(rho):
    """Coefficients of product_a (x^a+y^a+z^a), computed without Schur sums."""
    coefficients = {(0, 0, 0): 1}
    for a in rho:
        following = defaultdict(int)
        for exponent, value in coefficients.items():
            for i in range(3):
                target = list(exponent)
                target[i] += a
                following[tuple(target)] += value
        coefficients = dict(following)
    return coefficients


DELTA = (2, 1, 0)
ALTERNANT = tuple(
    ((-1) ** sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3)),
     tuple(DELTA[p[i]] for i in range(3)))
    for p in permutations(range(3)))


@cache
def character(mu, rho):
    """Three-variable Frobenius coefficient, matching the imported D5 owner."""
    require(sum(mu) == sum(rho) and len(mu) <= 3, "character domain")
    mu = mu + (0,) * (3 - len(mu))
    target = tuple(mu[i] + DELTA[i] for i in range(3))
    coefficients = power_sum_coefficients(rho)
    return sum(sign * coefficients.get(tuple(target[i] - shift[i]
                                            for i in range(3)), 0)
               for sign, shift in ALTERNANT)


def two_row_character(n, rho):
    """Independent two-variable coefficient difference for shape (n,n)."""
    coefficients = [1]
    for a in rho:
        following = coefficients + [0] * a
        for j, value in enumerate(coefficients):
            following[j + a] += value
        coefficients = following
    return coefficients[n] - (coefficients[n - 1] if n else 0)


@cache
def rim_hooks(mu, size):
    """Enumerate removable connected skew diagrams with no 2x2 square."""
    answer = []
    for nu in partitions(sum(mu) - size):
        if len(nu) > len(mu):
            continue
        padded = nu + (0,) * (len(mu) - len(nu))
        if any(padded[i] > mu[i] for i in range(len(mu))):
            continue
        cells = {(i, j) for i, a in enumerate(mu)
                 for j in range(padded[i], a)}
        if not cells or len(cells) != size:
            continue
        if any({(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)} <= cells
               for i, j in cells):
            continue
        reached = {next(iter(cells))}
        frontier = list(reached)
        while frontier:
            i, j = frontier.pop()
            for adjacent in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                if adjacent in cells and adjacent not in reached:
                    reached.add(adjacent)
                    frontier.append(adjacent)
        if reached == cells:
            answer.append((nu, (-1) ** (len({i for i, _ in cells}) - 1)))
    return tuple(answer)


@cache
def character_murnaghan_nakayama(mu, rho):
    """Separate rim-hook recursion, with no power-sum coefficient calculation."""
    if not rho:
        return int(not mu)
    return sum(sign * character_murnaghan_nakayama(nu, rho[1:])
               for nu, sign in rim_hooks(mu, rho[0]))


def main():
    expected_ev = Counter({(6, 4, 4): 1, (6, 4, 2, 2): 2,
                           (6, 2, 2, 2, 2): 1, (4, 4, 3, 3): 1,
                           (4, 3, 3, 2, 2): 2, (3, 3, 2, 2, 2, 2): 1})
    require(ev_types((3, 2, 2)) == expected_ev, "source equation (1.2)")
    require(r3(10) == ((10,), (8, 2), (6, 4), (6, 2, 2), (4, 4, 2)),
            "source equation (1.4)")
    print("FIDELITY equations (1.2) and (1.4): exact Ev multiplicities and R3(10) match")

    cycles = ((1,) * 8, (2,) + (1,) * 6, (2, 2) + (1,) * 4,
              (2, 2, 2, 1, 1), (2, 2, 2, 2))
    source_table = {(8,): (1, 1, 1, 1, 1), (6, 2): (20, 10, 4, 2, 4),
                    (4, 4): (14, 4, 2, 0, 6), (4, 2, 2): (56, 4, 0, 4, 8)}
    for mu, expected in source_table.items():
        actual = tuple(character(mu, rho) for rho in cycles)
        independent = tuple(character_murnaghan_nakayama(mu, rho) for rho in cycles)
        require(actual == independent == expected, f"source Table 1 row {mu}")
        print(f"FIDELITY Table 1 mu={mu}: {actual}")
    ev = ev_types((1, 1, 1, 1))
    left = sum(m * (-1) ** len(rho) * sum(character(mu, rho) for mu in r3(8))
               for rho, m in ev.items())
    right = sum(m * character((4, 4), rho) for rho, m in ev.items())
    require(left == right == 48, "source Section 2.1, both totals are 48")
    print(f"FIDELITY Section 2.1 lambda=(1,1,1,1): L={left}, R={right}")

    checked_pairs = set()
    partition_count = 0
    for n in range(13):
        left_average = Fraction(0)
        right_average = Fraction(0)
        normalization = Fraction(0)
        count = 0
        for lam in partitions(n):
            count += 1
            weight = Fraction(1, z_factor(lam))
            normalization += weight
            left = right = 0
            ev = ev_types(lam)
            require(sum(ev.values()) == 2 ** len(lam), "Ev multiplicity total")
            for rho, multiplicity in ev.items():
                ch = character((n, n), rho)
                require(ch == two_row_character(n, rho), "two-row character cross-check")
                right += multiplicity * ch
                for mu in r3(2 * n):
                    value = character(mu, rho)
                    left += multiplicity * (-1) ** len(rho) * value
                    if n <= 7:
                        pair = (mu, rho)
                        if pair not in checked_pairs:
                            require(value == character_murnaghan_nakayama(mu, rho),
                                    f"rim-hook cross-check {pair}")
                            checked_pairs.add(pair)
            require(left == right, f"pointwise sibling n={n}, lambda={lam}")
            left_average += weight * left
            right_average += weight * right
        require(normalization == 1, f"partition cycle-index normalization n={n}")
        expected = Fraction(comb(n // 2 + 2, 2) if n % 2 == 0 else 0)
        require(left_average == right_average == expected, f"weighted identity n={n}")
        partition_count += count
        print(f"n={n:2d} partitions={count:2d} left={left_average} "
              f"right={right_average} expected={expected}")
    print(f"EXACT weighted checks: all {partition_count} partitions for 0<=n<=12")
    print(f"INDEPENDENT rim-hook character comparisons: {len(checked_pairs)} distinct pairs for 0<=n<=7")
    print("INDEPENDENT two-variable character check: every Ev type for 0<=n<=12")
    print("BOUNDARY: n=0 is a generating-series convention; the conjecture quantifies n>=1")
    print("BOUNDARY: finite computations establish no universal theorem and no Lean proof")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
