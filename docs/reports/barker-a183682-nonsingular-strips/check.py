#!/usr/bin/env python3
"""Exact consistency checks for Hardin–Barker OEIS A183682; no Lean proof."""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ROWS = tuple(product((0, 1), repeat=3))


def compatible(r, s):
    """Both consecutive 2 by 2 determinants are nonzero in the integers."""
    return all(r[j] * s[j + 1] - r[j + 1] * s[j] != 0 for j in range(2))


def literal_count(n):
    """Enumerate all labelled arrays directly, independently of row transitions."""
    count = 0
    for entries in product((0, 1), repeat=3 * (n + 1)):
        if all(entries[3*i+j] * entries[3*(i+1)+j+1]
               - entries[3*i+j+1] * entries[3*(i+1)+j] != 0
               for i in range(n) for j in range(2)):
            count += 1
    return count


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def main():
    bfile = (ROOT / 'b183682.txt').read_bytes()
    terms = dict(tuple(map(int, line.split()))
                 for line in bfile.decode().splitlines()
                 if line.strip() and not line.lstrip().startswith('#'))
    assert tuple(terms) == tuple(range(1, 201))
    print('source: https://oeis.org/A183682/b183682.txt')
    print('source_bfile_sha256:', sha256(bfile).hexdigest())
    neighbors = {''.join(map(str, r)): tuple(''.join(map(str, s))
                 for s in ROWS if compatible(r, s)) for r in ROWS}
    expected = {
        '000': (), '001': (), '010': ('101', '111'),
        '011': ('101', '110'), '100': (),
        '101': ('010', '011', '110', '111'),
        '110': ('011', '101'), '111': ('010', '101'),
    }
    assert neighbors == expected
    print('all 64 ordered row pairs: exact neighbor table OK')
    for r, ns in neighbors.items():
        print(r + ': ' + ','.join(ns))
    # This uses ad-bc literally, rather than assuming equivalence to a permanent.
    for a, b, c, d in product((0, 1), repeat=4):
        assert (a*d-b*c != 0) == (a*d+b*c == 1)
    print('16 binary squares: determinant nonzero iff permanent one OK')
    active_a = {'010', '011', '110', '111'}
    active_b = {'101'}
    for r in active_a:
        assert len(set(neighbors[r]) & active_a) == 1
        assert len(set(neighbors[r]) & active_b) == 1
    assert set(neighbors['101']) == active_a
    print('class degrees: A -> (1 A, 1 B); B -> (4 A, 0 B) OK')
    q = ((1, 1), (4, 0))
    assert matmul(q, q) == ((5, 1), (4, 4))
    assert matmul(q, q) == tuple(tuple(q[i][j] + 4 * int(i == j)
                                      for j in range(2)) for i in range(2))
    print('quotient matrix: Q^2 = Q + 4 I OK')
    counts = [8]  # Literal one-row count; it is not an OEIS n=0 term.
    ways = {r: 1 for r in ROWS}
    u, v = 1, 1  # Auxiliary zero-length walks on the FIVE active vertices.
    for n in range(1, 201):
        ways = {r: sum(ways[s] for s in ROWS if compatible(r, s)) for r in ROWS}
        counts.append(sum(ways.values()))
        u, v = u + v, 4*u
        assert len({ways[r] for r in ROWS if ''.join(map(str, r)) in active_a}) == 1
        assert counts[n] == 4*u + v == terms[n]
    print('8-state propagation and 2-class quotient match all 200 source terms OK')
    assert counts[1:7] == [12, 32, 80, 208, 528, 1360]
    for n in range(1, 6):
        assert literal_count(n) == counts[n] == terms[n]
    print('direct labelled-array enumeration, n=1..5:', counts[1:6])
    assert counts[0] == 8 and 4*1 + 1 == 5
    assert counts[2] != counts[1] + 4*counts[0]
    for n in range(3, 201):
        assert counts[n] == counts[n-1] + 4*counts[n-2]
    print('source recurrence, n=3..200: OK; n=2 with literal count(0)=8 fails as expected')
    # Expand the displayed rational function by independent coefficient division.
    numerator = [Fraction(0), Fraction(12), Fraction(20)]
    denominator = [Fraction(1), Fraction(-1), Fraction(-4)]
    gf = []
    for n in range(201):
        rhs = numerator[n] if n < len(numerator) else Fraction(0)
        gf.append((rhs - sum(denominator[k]*gf[n-k]
                             for k in range(1, min(n, 2)+1))) / denominator[0])
    assert gf[0] == 0
    assert all(gf[n] == terms[n] for n in range(1, 201))
    for n in range(201):
        assert sum(denominator[k]*gf[n-k] for k in range(min(n, 2)+1)) == (
            numerator[n] if n < 3 else 0)
    print('4*x*(3+5*x)/(1-x-4*x^2): exact coefficients x^0..x^200 OK')
    print('finite consistency evidence only; all-height induction remains for Lean implementation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
