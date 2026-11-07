#!/usr/bin/env python3
"""Exact checks of literal arrays and Hardin's A204612 recurrence certificate."""
from collections import Counter
from itertools import product
import json
from pathlib import Path

SOURCE = json.loads(Path(__file__).with_name('source-data.json').read_text())
B = list(product((0, 1), repeat=5))


def delta(a, b):
    return tuple(a[j] * b[j + 1] - a[j + 1] * b[j] for j in range(4))


def horizontal(d):
    return all(d[j] <= d[j + 1] for j in range(3))


V = Counter((b, delta(a, b)) for a in B for b in B if horizontal(delta(a, b)))
S = sorted(V)
I = {s: i for i, s in enumerate(S)}
EDGES = [[I[(c, delta(b, c))] for c in B
          if horizontal(delta(b, c)) and all(d[j] <= delta(b, c)[j] for j in range(4))]
         for b, d in S]


def transfer(v):
    out = [0] * len(S)
    for source, destinations in enumerate(EDGES):
        for target in destinations:
            out[target] += v[source]
    return out


# Independent representation: integer-coded rows and literal integer determinants.
# Enumerate every actual matrix, including invalid ones; do not call delta or T.
def direct_count(n):
    def bit(r, j):
        return (r >> j) & 1
    blocks = {}
    for a in range(32):
        for b in range(32):
            d = tuple(bit(a, j) * bit(b, j + 1) - bit(a, j + 1) * bit(b, j)
                      for j in range(4))
            blocks[a, b] = d if d[0] <= d[1] <= d[2] <= d[3] else None
    count = 0
    for rows in product(range(32), repeat=n + 1):
        previous = None
        for i in range(n):
            d = blocks[rows[i], rows[i + 1]]
            if d is None or (previous is not None and
                             any(previous[j] > d[j] for j in range(4))):
                break
            previous = d
        else:
            count += 1
    return count


initial = [V[s] for s in S]
powers = [initial]
for _ in range(499):
    powers.append(transfer(powers[-1]))
counts = [sum(v) for v in powers]
assert len(S) == 68 and sum(initial) == 335 and sum(map(len, EDGES)) == 405
assert all(len(js) == len(set(js)) for js in EDGES)
print('STATES 68; INITIAL_MULTIPLICITY 335; TRANSITIONS 405')
for n in range(1, 4):
    value = direct_count(n)
    assert value == counts[n - 1] == SOURCE['displayed_terms'][n - 1]
    print(f'DIRECT n={n}; matrices={32 ** (n + 1)}; valid={value}')
assert counts[:16] == SOURCE['displayed_terms']
terms = {int(k): v for k, v in SOURCE['bfile_terms'].items()}
assert sorted(terms) == list(range(1, 211))
assert all(counts[n - 1] == v for n, v in terms.items())
print('SOURCE_FIDELITY displayed_terms=16; bfile_terms=210; all equal')

# P(x), with coefficients in increasing degree order.
p = [-24576, 9984, 40064, -20992, -18224, 12612, 1216, -2060, 256, 73, -18, 1]
w = [sum(p[j] * powers[4 + j][i] for j in range(12)) for i in range(len(S))]
table = [('01000', '0000', 8832), ('01001', '0000', 5124),
         ('01001', '0001', 4442), ('01010', '0000', 5150),
         ('01011', '0000', 5150), ('01011', '0001', 1854),
         ('01100', '0000', 96), ('01101', '0000', 96),
         ('01110', '0000', -13), ('01111', '0000', -13),
         ('01111', '0001', -13)]
z = [0] * len(S)
for bs, ds, q in table:
    b, d = tuple(map(int, bs)), tuple(map(int, ds))
    z[I[(b, d)]] -= 64 * q
    z[I[((1,) + b[1:], d)]] += 64 * q
assert w == z and any(w)
print('CERTIFICATE P(T)T^4v equals the specified 11-pair vector; nonzero=True')
pairs = [(i, I[((1,) + b[1:], d)]) for i, (b, d) in enumerate(S)
         if b[0] == 0 and b[1] == 1 and d[0] == 0]
assert len(pairs) == 11
covered = {i for pair in pairs for i in pair}
for lo, hi in pairs:
    g = [0] * len(S)
    g[hi], g[lo] = 1, -1
    y = transfer(g)
    assert sum(g) == sum(y) == 0
    assert all(y[i] == 0 for i in range(len(S)) if i not in covered)
    assert all(y[a] + y[b] == 0 for a, b in pairs)
print('INVARIANT_SUBSPACE all 11 generators mapped into paired zero-sum span')
residual_powers = [w]
for _ in range(68):
    residual_powers.append(transfer(residual_powers[-1]))
assert all(sum(y) == 0 for y in residual_powers)
assert residual_powers[4] == [4*w[i] - 2*residual_powers[1][i]
                            - 6*residual_powers[2][i] + 5*residual_powers[3][i]
                            for i in range(len(S))]
print('CORROBORATION sum(T^k w)=0 for k=0..68; auxiliary degree-4 vector identity holds')
for n in range(16, 501):
    assert sum(p[j] * counts[n - 12 + j] for j in range(12)) == 0
print('RECURRENCE all n=16..500; exact integers; original threshold retained')
print('ALL_CHECKS_OK')
