#!/usr/bin/env python3
"""Exact finite checks for Takuto Nakai's OEIS A400474 conjecture.

Integer arithmetic only. These checks do not establish the geometric bijection
or the theorem for all n. Source: https://oeis.org/A400474, revision 12.
"""
from fractions import Fraction
from pathlib import Path

CONSTANTS = (6, -3, 4, 3, 0, -5)


def closed_form(n):
    answer = Fraction(n * n - 4 * n + CONSTANTS[n % 6], 3)
    assert answer.denominator == 1 and answer >= 0
    return answer.numerator


def source_a(n):
    """Literal Python program in OEIS revision 12, apart from the function name."""
    def cnt(p, q, r):
        t = 0
        for x in range(p + 1, n):
            for y in range(q + 1, n):
                z = 2*n - x - y
                if r < z < n and x != y and y != z and x != z:
                    t += 1
        return t
    best = 0
    for s in (n - 2, n - 1, n):
        m = 0 if s < n else 1
        for p in range(m, n):
            for q in range(m, n):
                r = s - p - q
                if r >= m:
                    best = max(best, cnt(p, q, r))
    return best


def distinct_triples(n):
    return [(u, v, n-u-v) for u in range(1, n)
            for v in range(1, n-u) if len({u, v, n-u-v}) == 3]


def floor_vectors(n):
    for total in (n-2, n-1, n):
        lower = int(total == n)
        for p in range(lower, n):
            for q in range(lower, n):
                r = total-p-q
                if lower <= r < n:
                    yield p, q, r


def F(r):
    # Integer subtraction, including r=0; there is no natural truncation.
    return (r-1)**2 // 2


def H(n, r):
    # Remove u=v and u=w as well as v=w. Triple equality is already absent in F.
    return (F(r) - 2*max(0, r-n//2)
            + 2*int(n % 3 == 0 and 3*r >= 2*n))


def total_formula(n):
    return (n-1)*(n-2)//2 - 3*((n-1)//2) + 2*int(n % 3 == 0)


def balanced(n):
    b, t = divmod(n-2, 3)
    return (b, b + int(t == 2), b + int(t >= 1))


def residue_minimum(n):
    q, r = divmod(n, 6)
    return (2*(q-1)*(3*q-2), 2*(q-1)*(3*q-1), 6*q*(q-1),
            2*q*(3*q-2), 2*q*(3*q-1), 6*q*q)[r]


def main():
    # Fidelity to the published source program and examples.
    values = [source_a(n) for n in range(3, 31)]
    assert values[:12] == [0, 0, 0, 6, 6, 12, 16, 20, 24, 34, 38, 48]
    assert values == [closed_form(n) for n in range(3, 31)]
    assert source_a(6) == 6 and source_a(9) == 16
    print('FIDELITY: literal source program n=3..30 agrees; a(6)=6; a(9)=16')
    print('SOURCE_PREFIX n=3..14:', ','.join(map(str, values[:12])))

    # Stored published table: deterministic and usable without a network connection.
    table_path = Path(__file__).with_name('oeis-b400474.txt')
    table = [tuple(map(int, line.split())) for line in table_path.read_text().splitlines()
             if line.strip() and not line.startswith('#')]
    assert [n for n, _ in table] == list(range(3, 1001))
    assert all(closed_form(n) == v for n, v in table)
    print('FIDELITY: all 998 published b-file entries n=3..1000 equal the closed form')

    # Independent triple enumeration, exact removal counts, and optimization over
    # all realizable floor vectors (not merely the proposed balanced ones).
    floor_count = 0
    for n in range(3, 200):
        triples = distinct_triples(n)
        assert len(triples) == total_formula(n)
        histogram = [0]*n
        for u, _, _ in triples:
            histogram[u] += 1
        counts = [0]*n
        running = 0
        for r in range(1, n):
            running += histogram[n-r]
            counts[r] = running
        assert counts == [H(n, r) for r in range(n)]
        assert all(counts[r] == F(r) for r in range(n//2+1))
        best = -1
        for rs in floor_vectors(n):
            floor_count += 1
            # These rational angles have exactly the indicated floors and sum n.
            fractional = Fraction(n-sum(rs), 3)
            angles = [r+fractional for r in rs]
            assert all(0 < alpha < n and alpha.numerator//alpha.denominator == r
                       for alpha, r in zip(angles, rs))
            assert sum(angles) == n
            via_removals = len(triples)-sum(counts[r] for r in rs)
            assert via_removals >= 0
            if n <= 22:
                direct = sum(all(u < n-r for u, r in zip(uvw, rs)) for uvw in triples)
                assert direct == via_removals
                assert all(sum(u >= n-r for u, r in zip(uvw, rs)) <= 1
                           for uvw in triples)
            best = max(best, via_removals)
        assert best == closed_form(n), (n, best, closed_form(n))
        rs = balanced(n)
        assert sum(rs) == n-2 and max(rs) <= n//2
        assert sum(F(r) for r in rs) == residue_minimum(n)
        assert len(triples)-sum(F(r) for r in rs) == best
        # Donor-row lower bound and receiver-row upper bound used by the
        # general exchange argument; n<6 has no distinct positive triple.
        if n >= 6:
            for r in range(n//2+1, n-1):
                s = (n-2-r)//2
                assert counts[r]-counts[r-1] >= r-4 >= s
                assert counts[s+1]-counts[s] <= s
    print('IDENTITY: independent distinct-triple optimization n=3..199 agrees')
    print('FLOOR_VECTORS_CHECKED:', floor_count)
    print('REMOVALS: corrected H(n,r) for every 0<=r<n, n=3..199 agrees')
    print('DIRECT_SETS: disjointness and subtraction n=3..22 agree')
    print('EXCHANGE: donor/receiver inequalities n=6..199 agree')
    assert H(6, 5) == 6 and F(5) == 8
    print('THRESHOLD_REGRESSION: H(6,5)=6 whereas F(5)=8; F requires r<=floor(n/2)')

    for r in range(1000):
        assert F(r+1)-F(r) == 2*(r//2)
    assert closed_form(180) == 10562
    # Source's angle examples, in units of one degree for n=180.
    assert total_formula(180)-3*H(180, 60) == 10446
    alphas = (Fraction(597, 10), Fraction(597, 10), Fraction(606, 10))
    assert sum(alphas) == 180
    floors = [x.numerator//x.denominator for x in alphas]
    assert total_formula(180)-sum(H(180, r) for r in floors) == 10562
    print('DEGREE_CASE: equilateral=10446; (59.7,59.7,60.6) degrees=10562')
    print('GEOMETRY: circle-arc/angle-triple bijection is an unproved formalization obligation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
