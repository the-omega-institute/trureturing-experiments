#!/usr/bin/env python3
"""Exact checks for Wani, arXiv:2607.17677v1, Question 6.11.

All arithmetic is integer arithmetic.  The quotient vertices are the nonempty
proper supports; a support clique has capacity product(a_i) except that a
co-atom has capacity one.
"""
from itertools import combinations, combinations_with_replacement, product


def supports(t):
    return [frozenset(A) for r in range(1, t) for A in combinations(range(t), r)]


def maximal_cliques(t):
    """Bron--Kerbosch maximal cliques of the support intersection graph."""
    vs = supports(t)
    n = len(vs)
    nbr = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if vs[i] & vs[j]:
                nbr[i] |= 1 << j
                nbr[j] |= 1 << i

    out = []

    def bk(rmask, pmask, xmask):
        if not pmask and not xmask:
            out.append(rmask)
            return
        # Tomita pivot: branch only on P \ N(u).
        union = pmask | xmask
        pivot = (union & -union).bit_length() - 1
        best = (pmask & ~nbr[pivot]).bit_count()
        u = pivot
        bits = union & ~(1 << pivot)
        while bits:
            z = bits & -bits
            k = z.bit_length() - 1
            d = (pmask & ~nbr[k]).bit_count()
            if d < best:
                best, u = d, k
            bits ^= z
        candidates = pmask & ~nbr[u]
        while candidates:
            z = candidates & -candidates
            k = z.bit_length() - 1
            bk(rmask | z, pmask & nbr[k], xmask & nbr[k])
            pmask ^= z
            candidates ^= z
            xmask |= z

    bk(0, (1 << n) - 1, 0)
    return vs, out


_CLIQUE_CACHE = {}


def capacity(a, A, t):
    if len(A) == t - 1:
        return 1
    z = 1
    for i in A:
        z *= a[i]
    return z


def formula(a):
    t = len(a)
    U = frozenset(range(t))
    w = lambda A: capacity(a, A, t)
    M = t
    seen = set()
    for r in range(2, t - 1):
        for A in map(frozenset, combinations(range(t), r)):
            B = U - A
            key = tuple(sorted((tuple(sorted(A)), tuple(sorted(B)))))
            if key in seen:
                continue
            seen.add(key)
            M += max(w(A), w(B))
    S = []
    for i in range(t):
        S.append(t - 1 + sum(
            w(A)
            for r in range(1, t - 1)
            for A in map(frozenset, combinations(range(t), r))
            if i in A
        ))
    return max([M] + S), M, S


def exact_support_clique_number(a):
    t = len(a)
    if t not in _CLIQUE_CACHE:
        vs, cls = maximal_cliques(t)
        _CLIQUE_CACHE[t] = ([tuple(capacity([1] * t, A, t) for A in vs)], vs, cls)
    _, vs, cls = _CLIQUE_CACHE[t]
    weights = [capacity(a, A, t) for A in vs]
    return max(sum(weights[i] for i in range(len(vs)) if mask >> i & 1) for mask in cls)


def ring_elements(qs):
    return list(product(*(range(q) for q in qs)))


def mul(x, y):
    return tuple(a * b for a, b in zip(x, y))


def literal_difference_graph(qs):
    elems = ring_elements(qs)
    zero = tuple(0 for _ in qs)
    zstar = [x for x in elems if x != zero and any(v == 0 for v in x)]
    nonzero = [x for x in elems if x != zero]
    edges = set()
    for i, x in enumerate(zstar):
        for j in range(i + 1, len(zstar)):
            y = zstar[j]
            zero_product = mul(x, y) == zero
            weak = False
            annx = [r for r in nonzero if mul(r, x) == zero]
            anny = [s for s in nonzero if mul(s, y) == zero]
            for r in annx:
                if any(mul(r, s) == zero for s in anny):
                    weak = True
                    break
            if weak and not zero_product:
                edges.add((i, j))
    return zstar, edges


def support_graph(qs):
    elems = ring_elements(qs)
    zero = tuple(0 for _ in qs)
    zstar = [x for x in elems if x != zero and any(v == 0 for v in x)]
    edges = set()
    for i, x in enumerate(zstar):
        A = {k for k, v in enumerate(x) if v}
        for j in range(i + 1, len(zstar)):
            y = zstar[j]
            B = {k for k, v in enumerate(y) if v}
            if A & B and not (len(A) == len(B) == len(qs) - 1 and A == B):
                edges.add((i, j))
            elif A == B and len(A) <= len(qs) - 2:
                edges.add((i, j))
    return zstar, edges


def check_ring_fidelity():
    cases = ([2, 2, 2], [2, 3, 3], [2, 2, 2, 2], [2, 2, 2, 3], [2, 2, 2, 2, 2])
    counts = []
    for qs in cases:
        z1, e1 = literal_difference_graph(qs)
        z2, e2 = support_graph(qs)
        assert len(z1) == len(z2), (qs, len(z1), len(z2))
        assert e1 == e2, qs
        counts.append(len(z1))
    assert counts == [6, 13, 14, 21, 30], counts


def check_published_values():
    # The paper's t=3 and t=4 closed forms, and its three t=5 examples.
    assert formula((1, 2, 4))[0] == max(2, 3, 5) + 1
    assert formula((10, 10, 10, 1))[0] == 304
    assert formula((10, 10, 10, 10, 1))[0] == 4605
    assert formula((10, 10, 10, 1, 1))[0] == 1905


def source_closed_form(a):
    t = len(a)
    if t == 3:
        return max(a) + 2  # max(q_i) + 1
    if t == 4:
        stars = []
        tris = []
        for k in range(4):
            others = [j for j in range(4) if j != k]
            stars.append(a[k] * (1 + sum(a[j] for j in others)) + 3)
            tris.append(4 + sum(a[x] * a[y] for x, y in combinations(others, 2)))
        return max(stars + tris)
    raise ValueError(t)


def check_small_exhaustive():
    alphabet = (1, 2, 3, 4, 6, 7, 8, 10)
    checked = 0
    maximal_counts = {}
    for t in range(3, 7):
        if t not in _CLIQUE_CACHE:
            vs, cls = maximal_cliques(t)
            _CLIQUE_CACHE[t] = (None, vs, cls)
        maximal_counts[t] = len(_CLIQUE_CACHE[t][2])
        for a in combinations_with_replacement(alphabet, t):
            assert formula(a)[0] == exact_support_clique_number(a), (t, a)
            if t in (3, 4):
                assert formula(a)[0] == source_closed_form(a), (t, a, formula(a)[0], source_closed_form(a))
            checked += 1
    assert checked == 2958, checked
    assert maximal_counts == {3: 4, 4: 12, 5: 81, 6: 2646}, maximal_counts
    return checked, maximal_counts


def main():
    check_ring_fidelity()
    check_published_values()
    checked, counts = check_small_exhaustive()
    print("ring_cases=5 vertices=[6,13,14,21,30]")
    print("published_values=304,4605,1905")
    print(f"exhaustive_vectors={checked} maximal_support_cliques={counts}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
