#!/usr/bin/env python3
"""Exact finite-group checks for Legarreta--Sica--Tota--Gavioli Question 3.3."""
from __future__ import annotations
from itertools import product
from math import gcd


def semidirect_group(mod_i: int, mod_j: int, action: int):
    els = tuple((i, j) for i in range(mod_i) for j in range(mod_j))
    def mul(u, v):
        i, j = u; k, l = v
        return ((i + pow(action, j, mod_i) * k) % mod_i, (j + l) % mod_j)
    e = (0, 0)
    def inv(u):
        for v in els:
            if mul(u, v) == e and mul(v, u) == e:
                return v
        raise AssertionError("no inverse")
    return els, mul, inv, e


def direct_product(G, H):
    els1, mul1, inv1, e1 = G
    els2, mul2, inv2, e2 = H
    els = tuple((a, b) for a in els1 for b in els2)
    def mul(u, v): return (mul1(u[0], v[0]), mul2(u[1], v[1]))
    def inv(u): return (inv1(u[0]), inv2(u[1]))
    return els, mul, inv, (e1, e2)


def subgroup_generated(G, seed):
    els, mul, inv, e = G
    H = {e}
    queue = list(seed)
    while queue:
        g = queue.pop()
        if g in H:
            continue
        old = tuple(H)
        additions = [g, inv(g)]
        H.update(additions)
        # Saturate by products with the current frontier.  Every new element is
        # multiplied on both sides by the elements already present.
        frontier = list(set(additions))
        while frontier:
            x = frontier.pop()
            current = tuple(H)
            for y in current:
                for z in (mul(x, y), mul(y, x)):
                    if z not in H:
                        H.add(z)
                        frontier.append(z)
        _ = old
    return frozenset(H)


def all_subgroups(G):
    els, mul, inv, e = G
    subs = {frozenset({e})}
    queue = [frozenset({e})]
    while queue:
        H = queue.pop()
        for g in els:
            if g not in H:
                K = subgroup_generated(G, tuple(H) + (g,))
                if K not in subs:
                    subs.add(K); queue.append(K)
    return tuple(sorted(subs, key=lambda H: (len(H), sorted(H))))


def conjugate_set(G, H, g):
    els, mul, inv, e = G
    gi = inv(g)
    return frozenset(mul(mul(g, h), gi) for h in H)


def normalizer(G, H):
    return frozenset(g for g in G[0] if conjugate_set(G, H, g) == H)


def normalizer_data(G, subs=None):
    if subs is None:
        subs = all_subgroups(G)
    norms = {normalizer(G, H) for H in subs}
    classes = []
    unseen = set(norms)
    while unseen:
        N = next(iter(unseen))
        orb = {conjugate_set(G, N, g) for g in G[0]}
        classes.append(orb)
        unseen -= orb
    return subs, norms, classes


def centre(G):
    els, mul, inv, e = G
    return frozenset(g for g in els if all(mul(g, h) == mul(h, g) for h in els))


def assert_product_subgroups(G, H, prodG, prod_subs, A_subs, B_subs):
    expected = {frozenset((a, b) for a in A for b in B)
                for A in A_subs for B in B_subs}
    assert set(prod_subs) == expected
    # Coprime direct-product splitting gives all subgroups; each Cartesian
    # product is closed because both coordinates are subgroups.
    assert len(prod_subs) == len(A_subs) * len(B_subs)


def report_product(label, G, A_subs, B_subs, A, B):
    subs = tuple(frozenset((a, b) for a in H1 for b in H2)
                 for H1 in A_subs for H2 in B_subs)
    norms = set()
    for H1 in A_subs:
        for H2 in B_subs:
            norms.add(frozenset((a, b) for a in normalizer(A, H1) for b in normalizer(B, H2)))
    classes = []
    unseen = set(norms)
    while unseen:
        N = next(iter(unseen))
        orb = {conjugate_set(G, N, g) for g in G[0]}
        classes.append(orb)
        unseen -= orb
    Z = centre(G)
    print(f"{label}: order={len(G[0])} subgroups={len(subs)} normalizers={len(norms)} omega={len(classes)} centre={len(Z)} index={len(G[0])//len(Z)}")
    return subs, norms, classes


def report_with_subgroups(label, G, subs):
    _, norms, classes = normalizer_data(G, subs)
    Z = centre(G)
    print(f"{label}: order={len(G[0])} subgroups={len(subs)} normalizers={len(norms)} omega={len(classes)} centre={len(Z)} index={len(G[0])//len(Z)}")
    return subs, norms, classes


def report(label, G):
    return report_with_subgroups(label, G, all_subgroups(G))


def main():
    S = semidirect_group(8, 2, 3)
    P3 = semidirect_group(9, 3, 4)       # 1+p modulo p^2
    P5 = semidirect_group(25, 5, 6)
    assert all(len({S[1](g, h) for h in S[0]}) == len(S[0]) for g in S[0])
    s_subs, s_norms, s_classes = report("S", S)
    p3_subs, p3_norms, p3_classes = report("P_3", P3)
    p5_subs, p5_norms, p5_classes = report("P_5", P5)
    for p, P in [(3, P3), (5, P5)]:
        assert len(P[0]) == p**3
        assert len(centre(P)) == p
        assert len(p3_norms if p == 3 else p5_norms) == 2
        assert len(p3_classes if p == 3 else p5_classes) == 2
    G3 = direct_product(S, P3)
    g3_subs, g3_norms, g3_classes = report_product("S x P_3", G3, s_subs, p3_subs, S, P3)
    assert_product_subgroups(S, P3, G3, g3_subs, s_subs, p3_subs)
    assert len(g3_norms) == 10 and len(g3_classes) == 8
    assert len(centre(G3)) == 2 * 3
    assert len(G3[0]) // len(centre(G3)) == 8 * 3**2
    # The normalizers in S are exactly the five subgroups listed in the paper.
    # Their sizes are a compact exact check of the complete list.
    assert sorted(map(len, s_norms)) == [4, 4, 8, 8, 16]
    print("listed S normalizer orders:", sorted(map(len, s_norms)))
    print("ALL_CHECKS_OK")

if __name__ == "__main__":
    main()
