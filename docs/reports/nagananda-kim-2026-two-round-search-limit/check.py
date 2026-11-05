#!/usr/bin/env python3
"""Exact checks for arXiv:2605.15612v1, Section 4(i); no floating point.
Finite evidence for a general argument, not a computation offered as its proof.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import comb, isqrt


def weights(n, p):
    return tuple(p**e.bit_count() * (1-p)**(n-e.bit_count()) for e in range(1 << n))


def cells(states, batch):
    out = {}
    for e in states:
        out.setdefault(tuple(bool(e & a) for a in batch), []).append(e)
    return tuple(tuple(v) for v in out.values())


def best_success(n, w, groups):
    # Reporting any index is permitted; select the maximum-posterior one.
    return sum((max(sum((w[e] for e in g if e & (1 << i)), Q(0))
                    for i in range(n)) for g in groups), Q(0))


def frontier(points):
    by_success = {}
    for s, c in points:
        by_success[s] = min(c, by_success.get(s, c))
    out, best = [], None
    for s in sorted(by_success, reverse=True):
        c = by_success[s]
        if best is None or c < best:
            out.append((s, c))
            best = c
    return tuple(reversed(out))


def two_round_frontier(n, p):
    """Every reduced first batch, then DP over independent second-batch branches.
    Removing repeated/empty tests loses no information. Order within a batch
    is immaterial. Hence enumerate all subsets of the 2^n-1 nonempty tests.
    At fixed tests, a maximum-posterior decoder dominates every other decoder.
    Pareto pruning preserves the minimum cost subject to a success lower bound.
    """
    w = weights(n, p)
    batches = tuple(tuple(a+1 for a in range((1 << n)-1) if b & (1 << a))
                    for b in range(1 << ((1 << n)-1)))
    states, cache = tuple(range(1 << n)), {}

    def branch(g):
        if g not in cache:
            mass = sum((w[e] for e in g), Q(0))
            cache[g] = frontier((best_success(n, w, cells(g, b)), len(b)*mass)
                                for b in batches)
        return cache[g]

    out = []
    for first in batches:
        options = ((Q(0), Q(len(first))),)
        for g in cells(states, first):
            options = frontier((s+t, c+d) for s, c in options for t, d in branch(g))
        out.extend(options)
    ans = frontier(out)
    assert ans[0] == (p, 0)  # Untested guessing is allowed.
    assert ans[-1][0] == 1-(1-p)**n
    return ans, len(batches), len(cache)


def optimum(points, gamma, randomized=False):
    costs = [c for s, c in points if s >= gamma]
    if randomized:
        # The constrained convex-envelope optimum uses at most two vertices.
        for (s, c), (t, d) in combinations(points, 2):
            if s < gamma < t:
                theta = (gamma-s)/(t-s)
                costs.append((1-theta)*c + theta*d)
    return min(costs)


def lower_bound(p, gamma):
    ell = 0
    while (1 << (ell+1))*p <= 1:
        ell += 1
    # Sum gamma - 2^t p for t=0,...,floor(log_2(1/p))-1.
    return ell*gamma - p*((1 << ell)-1)


def construction(n, lam, gamma):
    p, L = lam/n, (n-1).bit_length()
    m = max(1, isqrt(L) + (isqrt(L)**2 < L))
    target = gamma + lam*lam/m
    h = next((h for h in range(n+1) if 1-(1-p)**h >= target), None)
    if h is None:
        return None  # This is an asymptotic construction; no finite claim here.
    q = 1-(1-p)**h
    assert target <= q < target+p
    a, r = divmod(h, m)
    sizes = [a + (j < r) for j in range(m)]
    blocks, start = [], 0
    for size in sizes:
        blocks.append(((1 << size)-1) << start)
        start += size
    assert start == h
    collision = p*p*sum(comb(s, 2) for s in sizes)
    assert collision <= lam*lam/(2*m)
    success = cost = singleton_success = Q(0)
    for e, mass in enumerate(weights(n, p)):
        positive = [b for b in blocks if e & b]
        cost += mass*(m + (L if positive else 0))
        if not positive:
            continue
        b, decoded = positive[0], 0
        # Bit-query answers are the bitwise OR of the excellent labels in b.
        for i in range(n):
            if e & b & (1 << i):
                decoded |= i
        good = decoded < n and bool(e & (1 << decoded))
        if good:
            success += mass
        if (e & b).bit_count() == 1:
            singleton_success += mass
            assert good
    assert cost == m+q*L
    assert success >= singleton_success >= q-collision >= gamma
    assert cost <= gamma*L + m + (lam*lam/m+p)*L
    return h, m, success, cost


def main():
    fidelity = 0
    for n in range(1, 9):
        for lam in (Q(1, 4), Q(1, 2), Q(1)):
            p = lam/n
            w = weights(n, p)
            assert sum(w) == 1 and w[0] == (1-p)**n
            for k in range(n+1):
                assert sum((w[e] for e in range(1 << n) if e.bit_count() == k), Q(0)) \
                    == comb(n, k)*p**k*(1-p)**(n-k)
            q = 1-(1-p)**n
            assert sum(w[1:]) == q
            # Source construction identity, with singleton tests as the second batch.
            assert sum((w[e]*(1+(n if e else 0)) for e in range(1 << n)), Q(0)) == 1+q*n
            fidelity += 1
    print('FIDELITY: source binomial law, finite-n feasibility, expected-cost identity:', fidelity, 'pairs OK')
    checks = 0
    for n in range(1, 4):
        for p in (Q(1, 8), Q(1, 4), Q(1, 2)):
            points, batches, states = two_round_frontier(n, p)
            q = 1-(1-p)**n
            print(f'FRONTIER n={n} p={p} batches={batches} branch-states={states}:',
                  '; '.join(f'({s},{c})' for s, c in points))
            for gamma in sorted({p, q/2, 3*q/4, q}):
                lo, det = lower_bound(p, gamma), optimum(points, gamma)
                rnd = optimum(points, gamma, True)
                assert lo <= rnd <= det
                print(f'OPT n={n} p={p} gamma={gamma}: lower={lo} deterministic={det} randomized={rnd}')
                checks += 1
    print('CONVERSE: ALL two-round deterministic procedures and randomized convex envelopes:', checks, 'cases OK')
    done = 0
    for n in range(2, 15):
        for lam, gamma in ((Q(1, 4), Q(1, 16)), (Q(1, 2), Q(1, 8))):
            ans = construction(n, lam, gamma)
            if ans is None:
                continue
            h, m, success, cost = ans
            print(f'CONSTRUCTION n={n} lambda={lam} gamma={gamma} h={h} m={m}: success={success} cost={cost}')
            done += 1
    assert done > 0
    print('CONSTRUCTION: exhaustive configurations, singleton decoding, unconditional success, collision and upper bounds:', done, 'cases OK')
    balanced = 0
    for h in range(65):
        for m in range(1, 65):
            a, r = divmod(h, m)
            twice_pairs = r*(a+1)*a + (m-r)*a*(a-1)
            assert twice_pairs*m <= h*h
            balanced += 1
    print('BALANCED BLOCKS: exact integer pair-count bound:', balanced, 'cases OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
