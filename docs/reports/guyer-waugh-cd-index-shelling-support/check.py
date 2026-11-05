#!/usr/bin/env python3
"""Exact finite-poset fidelity and support checks for arXiv:2610.25491v1, Q3.16.
Only Python integers are used. No claimed support pattern enters the constructors.
"""
from collections import defaultdict
from functools import lru_cache
from math import comb


def add(*polys):
    out = defaultdict(int)
    for p in polys:
        for w, a in p.items():
            out[w] += a
    return {w: a for w, a in out.items() if a}


def scale(p, a):
    return {w: a * b for w, b in p.items() if a * b}


def append(p, letter):
    return {w + letter: a for w, a in p.items()}


@lru_cache(None)
def words(n):
    if n == 0:
        return ('',)
    return tuple('c' + w for w in words(n - 1)) + (
        tuple('d' + w for w in words(n - 2)) if n >= 2 else ())


@lru_cache(None)
def expand_word(w):
    p = {'': 1}
    for x in w:
        q = defaultdict(int)
        for u, a in p.items():
            for v in (('a', 'b') if x == 'c' else ('ab', 'ba')):
                q[u + v] += a
        p = dict(q)
    return p


def decode_cd(p, n):
    """Triangular exact substitution inverse; also verifies zero residual."""
    residual = defaultdict(int, p)
    result = {}
    for w in sorted(words(n), key=lambda w: min(expand_word(w))):
        expansion = expand_word(w)
        pivot = min(expansion)
        assert expansion[pivot] == 1
        a = residual[pivot]
        if a:
            result[w] = a
            for v, b in expansion.items():
                residual[v] -= a * b
    assert all(a == 0 for a in residual.values()), (n, 'not a cd-index', residual)
    return result


def cumulative_from_poset(n, k):
    """cd-index of semisuspension of the first k facets of the n-simplex.

    A facet omits one of n+1 vertices. Selected omitted vertices are 0,...,k-1.
    Old nonempty faces miss a selected vertex. None is the new cap facet.
    An old face lies below the cap iff it also misses an unselected vertex.
    Empty face and greatest element are excluded from flag chains.
    """
    assert 1 <= k <= n
    all_vertices = (1 << (n + 1)) - 1
    selected = (1 << k) - 1
    unselected = all_vertices ^ selected
    levels = {r: [] for r in range(1, n + 1)}
    for s in range(1, all_vertices):
        if selected & (all_vertices ^ s):
            levels[s.bit_count()].append(s)
    levels[n].append(None)

    def below(s, t):
        if t is None:
            return s is None or bool((all_vertices ^ s) & unselected)
        return s is not None and (s & t) == s

    flag_f = []
    for mask in range(1 << n):
        ranks = [r + 1 for r in range(n) if mask & (1 << r)]
        if not ranks:
            flag_f.append(1)
            continue
        dp = {s: 1 for s in levels[ranks[0]]}
        for r in ranks[1:]:
            dp = {t: sum(a for s, a in dp.items() if below(s, t)) for t in levels[r]}
        flag_f.append(sum(dp.values()))
    flag_h = flag_f.copy()
    for r in range(n):
        for mask in range(1 << n):
            if mask & (1 << r):
                flag_h[mask] -= flag_h[mask ^ (1 << r)]
    ab = {''.join('b' if mask & (1 << r) else 'a' for r in range(n)): a
          for mask, a in enumerate(flag_h) if a}
    return decode_cd(ab, n)


@lru_cache(None)
def recurrence(n, i):
    """Corollary 3.10, independently of the finite-poset constructor."""
    assert 0 <= i <= n
    if n == 0:
        return {'': 1}
    if i == n:
        return {}
    cpart = append(add(*(recurrence(n - 1, j) for j in range(i, n))), 'c')
    if i == 0:
        return cpart
    h = [sum(1 for a in range(i) for b in range(n - i) if a + b == j)
         for j in range(n - 1)]
    dpart = append(add(*(scale(recurrence(n - 2, j), h[j])
                         for j in range(n - 1))), 'd')
    return add(cpart, dpart)


def expected_support(n, i, w):
    if i == 0:
        return w.endswith('c')
    if i == n or 'd' not in w:
        return False
    t = len(w) - len(w.rstrip('c'))
    return i + t + 1 <= n


def main():
    fidelity = 0
    for n in range(1, 8):
        previous = {}
        for k in range(1, n + 1):
            cumulative = cumulative_from_poset(n, k)
            contribution = add(cumulative, scale(previous, -1))
            assert contribution == recurrence(n, k - 1), (n, k, contribution)
            previous = cumulative
            fidelity += 1
    print(f'FIDELITY: {fidelity} semisuspended initial segments (1 <= n <= 7) match Corollary 3.10')
    print('FIDELITY: H(2,0)=cc, H(2,1)=d, H(2,2)=0')
    assert [recurrence(2, i) for i in range(3)] == [{'cc': 1}, {'d': 1}, {}]
    assert recurrence(0, 0) == {'': 1}
    support_checks = identity_checks = known_cases = 0
    for n in range(1, 15):
        for i in range(n + 1):
            p = recurrence(n, i)
            assert all(a > 0 for a in p.values())
            assert all(sum(1 if x == 'c' else 2 for x in w) == n for w in p)
            for w in words(n):
                assert (p.get(w, 0) > 0) == expected_support(n, i, w), (n, i, w, p)
                support_checks += 1
                if i == 1 and n >= 2:
                    assert (w in p) == (w != 'c' * n)
                    known_cases += 1
                if i == 2 and n >= 3:
                    assert (w in p) == (w not in ('c' * n, 'd' + 'c' * (n - 2)))
                    known_cases += 1
                if 0 < i < n and 'd' in w:
                    t = len(w) - len(w.rstrip('c'))
                    if t == 0:
                        assert p.get(w, 0) > 0
                        continue
                    r = n - t
                    ud = w[:-t]
                    rhs = sum(comb(j - i + t - 1, t - 1) * recurrence(r, j).get(ud, 0)
                              for j in range(i, r))
                    assert p.get(w, 0) == rhs, (n, i, w, p.get(w, 0), rhs)
                    identity_checks += 1
    print(f'FIDELITY: {known_cases} word tests reproduce Corollaries 3.17 and 3.18')
    print(f'SUPPORT: {support_checks} (n,i,word) cases, every 1 <= n <= 14 and 0 <= i <= n')
    print(f'IDENTITY: {identity_checks} trailing-c coefficient identities, every applicable case through degree 14')
    print('ARITHMETIC: Python integers only; empty sums are 0; H(0,0)=1')
    print('SCOPE: finite evidence, not a proof for unbounded n')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
