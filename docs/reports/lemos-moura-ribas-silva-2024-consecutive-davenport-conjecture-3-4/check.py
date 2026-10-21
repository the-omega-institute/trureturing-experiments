"""Lemos-Moura-Ribas-Silva, arXiv:2404.11312v1, Conjecture 3.4: C(G) = |G| for every finite group G.

C(G) is the least positive l such that every ordered sequence over G of length >= l has a nonempty block
g_i g_(i+1) ... g_j equal to 1. For the nonabelian groups S_3, D_4, Q_8 (and the abelian C_2 x C_2, C_6 for
comparison) the script computes C(G) by exhaustive search over all sequences of length |G| - 1 and |G|, and checks
that the prefix-quotient word g_i = h_(i-1)^(-1) h_i, built from an enumeration h_0 = 1, h_1, ... of G, is
product-one free of length |G| - 1. It also counts the product-one free words of length |G| - 1, which are
exactly the words whose prefix products enumerate G, so there are (|G| - 1)! of them. For A_4, S_4 and the nonabelian group of order 21 it checks the construction only.
"""
import itertools
import math
import sys


def perm_group(gens, n):
    idt = tuple(range(n))
    G = {idt}
    frontier = [idt]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = tuple(g[s[i]] for i in range(n))
                if h not in G:
                    G.add(h)
                    new.append(h)
        frontier = new
    mul = lambda a, b: tuple(a[b[i]] for i in range(n))
    return sorted(G), mul, idt


def quaternion():
    # elements as (sign, unit) with unit in {1, i, j, k}
    table = {('1', '1'): (1, '1'), ('1', 'i'): (1, 'i'), ('1', 'j'): (1, 'j'), ('1', 'k'): (1, 'k'),
             ('i', '1'): (1, 'i'), ('i', 'i'): (-1, '1'), ('i', 'j'): (1, 'k'), ('i', 'k'): (-1, 'j'),
             ('j', '1'): (1, 'j'), ('j', 'i'): (-1, 'k'), ('j', 'j'): (-1, '1'), ('j', 'k'): (1, 'i'),
             ('k', '1'): (1, 'k'), ('k', 'i'): (1, 'j'), ('k', 'j'): (-1, 'i'), ('k', 'k'): (-1, '1')}
    G = [(s, u) for s in (1, -1) for u in '1ijk']

    def mul(a, b):
        s, u = table[(a[1], b[1])]
        return (a[0] * b[0] * s, u)
    return G, mul, (1, '1')


def free(word, mul, e):
    for i in range(len(word)):
        p = e
        for j in range(i, len(word)):
            p = mul(p, word[j])
            if p == e:
                return False
    return True


def inverse(g, G, mul, e):
    return next(h for h in G if mul(g, h) == e)


def construction(G, mul, e):
    order = [e] + [g for g in G if g != e]
    return [mul(inverse(order[i - 1], G, mul, e), order[i]) for i in range(1, len(order))]


def C_exhaustive(G, mul, e):
    """Returns (C(G) if it equals |G|, number of product-one free words of length |G| - 1)."""
    n = len(G)
    n_free = sum(1 for w in itertools.product(G, repeat=n - 1) if free(w, mul, e))
    all_hit = all(not free(w, mul, e) for w in itertools.product(G, repeat=n))
    return (n if (n_free > 0 and all_hit) else None), n_free


groups = {
    "C_2 x C_2": perm_group([(1, 0, 3, 2), (2, 3, 0, 1)], 4),
    "C_6": perm_group([(1, 2, 3, 4, 5, 0)], 6),
    "S_3": perm_group([(1, 0, 2), (1, 2, 0)], 3),
    "D_4": perm_group([(1, 2, 3, 0), (3, 2, 1, 0)], 4),
    "Q_8": quaternion(),
}
big = {
    "A_4": perm_group([(1, 2, 0, 3), (1, 0, 3, 2)], 4),
    "S_4": perm_group([(1, 0, 2, 3), (1, 2, 3, 0)], 4),
    "C_7 x| C_3": perm_group([(1, 2, 3, 4, 5, 6, 0), (0, 2, 4, 6, 1, 3, 5)], 7),
}
ok = True
for name, (G, mul, e) in groups.items():
    c, n_free = C_exhaustive(G, mul, e)
    w = construction(G, mul, e)
    good = c == len(G) and len(w) == len(G) - 1 and free(w, mul, e) and n_free == math.factorial(len(G) - 1)
    ok &= good
    print(f"{name}: |G|={len(G)} C(G) by exhaustive search={c} free words of length |G|-1: {n_free} "
          f"(|G|-1)! = {math.factorial(len(G) - 1)} construction free={free(w, mul, e)} ok={good}")
for name, (G, mul, e) in big.items():
    w = construction(G, mul, e)
    good = len(w) == len(G) - 1 and free(w, mul, e)
    ok &= good
    print(f"{name}: |G|={len(G)} construction of length {len(w)} product-one free={free(w, mul, e)} ok={good}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
