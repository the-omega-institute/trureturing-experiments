#!/usr/bin/env python3
"""Exact finite checks for the twin/separator reduction; Python standard library only."""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from hashlib import sha256
from pathlib import Path


def factor(n):
    ans = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            ans[p] = ans.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        ans[n] = ans.get(n, 0) + 1
    return ans


def formal_log(q):
    """Coefficients of log_2(prime); exact multiplicative rational identities."""
    assert q > 0
    out = defaultdict(F)
    for p, e in factor(q.numerator).items():
        out[p] += e
    for p, e in factor(q.denominator).items():
        out[p] -= e
    return dict(out)


def linear(*terms):
    out = defaultdict(F)
    for c, v in terms:
        for p, a in v.items():
            out[p] += c * a
    return {p: a for p, a in out.items() if a}


def entropy_formal(p):
    return linear(*[(-x, formal_log(x)) for x in p if x])


def fidelity():
    count = 0
    for den in range(3, 13):
        for a in range(1, den - 1):
            for b in range(1, den - a):
                p = [F(a, den), F(b, den), F(den - a - b, den)]
                for i in range(3):
                    for j in range(i + 1, 3):
                        x, y = p[i], p[j]
                        q = [p[k] for k in range(3) if k not in (i, j)] + [x + y]
                        loss = linear((1, entropy_formal(p)), (-1, entropy_formal(q)))
                        source_formula = linear((x, formal_log(1 / x)),
                                                (y, formal_log(1 / y)),
                                                (-(x + y), formal_log(1 / (x + y))))
                        assert loss == source_formula
                        count += 1
    assert entropy_formal([F(1, 4), F(1, 4), F(1, 2)]) == {2: F(3, 2)}
    assert entropy_formal([F(1, 2), F(1, 2)]) == {2: F(1)}
    print(f'FIDELITY: source Lemma 1 merge-loss identity, {count} exact rational cases')


def ln_unit(u, m):
    assert 1 <= u <= 2
    t = (u - 1) / (u + 1)
    z, s = t, F(0)
    for j in range(m):
        s += 2 * z / (2 * j + 1)
        z *= t * t
    tail = 2 * z / ((2 * m + 1) * (1 - t * t))
    return s, s + tail


def entropy_bounds(x, scale, m, ln2):
    """Integer bounds: lo/scale <= -x log_2 x <= hi/scale."""
    if x == 0:
        return 0, 0
    assert x > 0
    u, e = x, 0
    while u < 1:
        u *= 2
        e -= 1
    while u >= 2:
        u /= 2
        e += 1
    a, b = ln_unit(u, m)
    l2, u2 = ln2
    lo, hi = -x * (e + b / l2), -x * (e + a / u2)
    return (lo * scale).__floor__(), (hi * scale).__ceil__()


def gadget(a, target):
    n, S = len(a), sum(a)
    assert n > 0 and all(x > 0 for x in a) and 1 <= target <= S
    N, M = 2 * n, S + 1
    weights = [M + x for x in a] + [M] * n
    W, B, K = sum(weights), n * M + target, 8 * N
    c = F(1, M * 2**K)
    xs = [c * w / 2 for w in weights]
    b = (1 - c * W) / (N - 1)
    p = []
    for i, x in enumerate(xs):
        if i:
            p.append(b)
        p.extend([x, x])
    assert sum(p) == 1 and all(x > 0 for x in p)
    assert B < W and c * W < F(1, 2)
    assert 2 * (N - 1) <= 2**N
    for x in xs:
        assert b / x >= 2**(7 * N)
        assert 7 * N * x > c * W
    return weights, W, B, K, c, p


def sums(a):
    return {sum(x for x, use in zip(a, mask) if use)
            for mask in product([False, True], repeat=len(a))}


def forced_cardinality():
    count = 0
    for n in range(1, 4):
        for a in product(range(1, 5), repeat=n):
            source = sums(a)
            for target in range(1, sum(a) + 1):
                weights, W, B, K, c, p = gadget(a, target)
                attainable = sums(weights)
                assert (target in source) == (B in attainable)
                optimum = min(u for u in attainable if u >= B)
                assert (optimum == B) == (target in source)
                count += 1
    print(f'REDUCTION: every positive target for all a_i in 1..4, n in 1..3: {count} instances')


def exhaustive(a, target, expected):
    weights, W, B, K, c, p = gadget(a, target)
    bits = K + (sum(a) + 1).bit_length() + 40
    scale = 2**bits
    m = (bits + 2) // 3 + 12
    ln2 = ln_unit(F(2), m)
    assert ln2[0] > F(1, 2)
    bounds = {}
    L = len(p)
    for i in range(L):
        mass = F(0)
        for j in range(i + 1, L + 1):
            mass += p[j - 1]
            bounds[i, j] = entropy_bounds(mass, scale, m, ln2)
    hp_lo = sum(bounds[i, i + 1][0] for i in range(L))
    hp_hi = sum(bounds[i, i + 1][1] for i in range(L))
    h = F(hp_lo + hp_hi, 2 * scale)
    assert F(hp_hi - hp_lo, 2 * scale) <= c / 16
    R = h - c * B + c / 8
    lower_threshold = h - c * B - c / 8
    assert 0 < lower_threshold < R
    best_u = min(u for u in sums(weights) if u >= B)
    assert best_u == expected
    best_entropy_lo = F(hp_lo) - c * best_u * scale
    feasible_count = good_count = crossing_count = 0
    optimum_masks = []
    window_answer = False
    for mask in range(1 << (L - 1)):
        cuts = [0] + [i + 1 for i in range(L - 1) if mask >> i & 1] + [L]
        blocks = list(zip(cuts, cuts[1:]))
        lo = sum(bounds[i, j][0] for i, j in blocks)
        hi = sum(bounds[i, j][1] for i, j in blocks)
        assert hi <= R * scale or lo > R * scale, 'undecided feasibility interval'
        feasible = hi <= R * scale
        feasible_count += feasible
        good = all(j == i + 1 or (i % 3 == 0 and j == i + 2) for i, j in blocks)
        if good:
            good_count += 1
            U = sum(weights[i // 3] for i, j in blocks if j == i + 2)
            exact_lo, exact_hi = F(hp_lo) - c * U * scale, F(hp_hi) - c * U * scale
            assert lo <= exact_hi and exact_lo <= hi
            assert feasible == (U >= B)
            if feasible and U == best_u:
                optimum_masks.append(mask)
            assert (lower_threshold * scale <= exact_lo and exact_hi <= R * scale) == (U == B)
            window_answer |= (U == B)
        else:
            crossing_count += 1
            assert hi < F(hp_lo) - c * W * scale, 'separator crossing not dominated'
            assert hi < best_entropy_lo
    assert good_count == 2**len(weights)
    assert optimum_masks and window_answer == (target in sums(a))
    print(f'EXHAUSTIVE: a={list(a)}, T={target}, B={B}, optimum_U={best_u}, '
          f'partitions={1 << (L - 1)}, good={good_count}, crossing={crossing_count}, '
          f'optimizers={len(optimum_masks)}, feasible={feasible_count}, window_yes={window_answer}')


def main():
    fidelity()
    forced_cardinality()
    for a, t, expected in [([1], 1, 3), ([1, 3], 2, 13), ([1, 3], 3, 13),
                           ([2, 3, 7], 5, 44), ([2, 4, 8], 5, 51)]:
        exhaustive(a, t, expected)
    print('Arithmetic: integers, Fraction, and exact prime-log coefficient identities; no floating point')
    print('Scope: finite checks of the proposed reduction; no complexity theorem or Lean proof')
    print('check.py sha256=' + sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
