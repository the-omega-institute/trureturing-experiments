#!/usr/bin/env python3
"""Fixed n=3 certificate replay; no hypothesis search, imports, or b enumeration.

Source: Chau et al., The stable Tamari lattice, arXiv:2610.39257v1,
Definitions 2.22/2.30, Proposition 2.24, Theorem 4.17; CC BY 4.0
(source license). This replay program is original.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
CERTIFICATE_DATA = ROOT / 'context_exchange_data.json'

def exponent(left, right, colored=0):
    return tuple(int(left <= j <= right) for j in (1, 2, 3)) + (colored,)

def weight(h):
    intervals = tuple((h[r-1], r) for r in (1, 2, 3))
    assert all(1 <= l <= r for l, r in intervals)
    for (l, r), (s, v) in product(intervals, repeat=2):
        assert r < s or v < l or l <= s <= v <= r or s <= l <= r <= v
    a = (0,) + tuple(sum(l <= i <= r for l, r in intervals)-1 for i in (1, 2, 3))
    U = []
    for j in (1, 2, 3):
        if h[j-1] == j:
            for gamma in range(a[j-1], a[j]+1):
                k = max(v for v in range(j, 4)
                        if all(a[r] >= gamma for r in range(j, v+1)))
                U.append((j, gamma, k))
    factors = [exponent(j, k) for j, gamma, k in U]
    factors += [exponent(l, r, 1) for l, r in intervals]
    return {'h': h, 'a': a[1:], 'U': U, 'numerator': a[1:] + (0,), 'factors': factors}

def swapped(e):
    return (e[1], e[0], e[2], e[3])

def expanded(numerator, factors, sign=1):
    # Direct subset expansion of the whole product.
    ans = Counter()
    for selected in product((0, 1), repeat=len(factors)):
        e = tuple(numerator[j] + sum(use*f[j] for use, f in zip(selected, factors))
                  for j in range(4))
        ans[e] += sign * (-1)**sum(selected)
    return {e: c for e, c in ans.items() if c}

def value(e, values):
    ans = Fraction(1)
    for v, power in zip(values, e):
        ans *= v**power
    return ans

def main():
    w0, w1 = weight((1, 2, 3)), weight((1, 1, 3))
    assert w0['a'] == (0, 0, 0) and w1['a'] == (1, 0, 0)
    assert w0['U'] == [(1, 0, 3), (2, 0, 3), (3, 0, 3)]
    assert w1['U'] == [(1, 0, 3), (1, 1, 1), (3, 0, 3)]
    # The single untouched h interval is [3,3], disjoint from [1,2].
    assert (w0['h'][2], 3) == (3, 3) == (w1['h'][2], 3)
    terms = [(1, w['numerator'], w['factors']) for w in (w0, w1)]
    terms += [(-1, swapped(w['numerator']), list(map(swapped, w['factors'])))
              for w in (w0, w1)]
    # Independently declare the displayed Q before comparing with stored data.
    Q = [exponent(1, 3), exponent(2, 3), exponent(3, 3),
         exponent(1, 1, 1), exponent(2, 2, 1), exponent(3, 3, 1),
         exponent(1, 1), exponent(1, 2, 1), (1, 0, 1, 0), exponent(2, 2)]
    assert Q[8] == (1, 0, 1, 0)  # 1-xz
    cleared = Counter()
    values = (Fraction(1, 5), Fraction(1, 7), Fraction(1, 11), Fraction(2))
    numeric_terms = []
    for sign, numerator, factors in terms:
        residual = list(Q)
        for f in factors:
            residual.remove(f)
        for e, c in expanded(numerator, residual, sign).items():
            cleared[e] += c
        v = value(numerator, values)
        for f in factors:
            denominator = 1-value(f, values)
            assert denominator
            v /= denominator
        numeric_terms.append(sign*v)
    cleared = {e: c for e, c in cleared.items() if c}
    cancel = [exponent(1, 3), exponent(3, 3), exponent(1, 2, 1)]
    expected = Counter(expanded((1, 0, 0, 0), cancel))
    for e, c in expanded((0, 1, 0, 0), cancel, -1).items():
        expected[e] += c
    expected = {e: c for e, c in expected.items() if c}
    assert cleared == expected and cleared[(1, 0, 0, 0)] == 1
    reduced = list(Q)
    for f in cancel:
        reduced.remove(f)
    reduced_value = values[0]-values[1]
    for f in reduced:
        reduced_value /= 1-value(f, values)
    lhs = numeric_terms[0]+numeric_terms[1]
    rhs = -numeric_terms[2]-numeric_terms[3]
    assert lhs-rhs == reduced_value == Fraction(326095, 1329696)
    # Compare only after the independently constructed identity is established.
    stored = json.loads(CERTIFICATE_DATA.read_text())['witness']
    stored_poly = {tuple(v['exponents']): v['coefficient'] for v in stored['cleared_numerator']}
    assert cleared == stored_poly
    assert str(lhs) == stored['numeric_check']['lhs']
    assert str(rhs) == stored['numeric_check']['rhs']
    report = {
      'source': {'arxiv': '2610.39257v1', 'sections': ['2.22', '2.24', '2.30', '4.17'],
                 'license': 'CC BY 4.0'},
      'weights': [w0, w1],
      'Q_factors': Q,
      'cleared_polynomial': [{'exponents': e, 'coefficient': c} for e, c in sorted(cleared.items())],
      'certificate': 'Q D = (x-y)(1-xyz)(1-z)(1-txy)',
      'nonzero_coefficient': {'monomial': 'x', 'coefficient': 1},
      'reduced_D': '(x-y)/((1-tx)(1-ty)(1-tz)(1-x)(1-y)(1-xz)(1-yz))',
      'numeric': {'lhs': str(lhs), 'rhs': str(rhs), 'D': str(lhs-rhs)},
      'checks': {'full_polynomial_identity': True, 'formal_series_units': True,
                 'stored_certificate_agrees': True, 'kernel_verification': False},
    }
    dest = ROOT / 'independent_review_data.json'
    dest.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'certificate_verified': True, 'numeric': report['numeric']}))

if __name__ == '__main__':
    main()
