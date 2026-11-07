#!/usr/bin/env python3
"""Exact checks for Elek, arXiv:2610.22730v1, Conjecture 3.3.

All arithmetic is integer arithmetic.  The first check reproduces the paper's
S_4 magic-word example.  The second checks the B_2 counterexample and the
rank-two I_2(m) family by exact signed-matrix/dihedral normal-form arithmetic.
"""
from itertools import product


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for i in range(len(a)))


def ident(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def word_prod(gens, word):
    out = ident(len(gens[0]))
    for i in word:
        out = mm(out, gens[i])
    return out


def all_words_products(gens, max_len):
    out = {ident(len(gens[0])): {()}}
    for n in range(1, max_len + 1):
        for w in product(range(len(gens)), repeat=n):
            out.setdefault(word_prod(gens, w), set()).add(w)
    return out


def lengths_and_words(gens, max_len):
    allp = all_words_products(gens, max_len)
    lengths = {g: min(len(w) for w in ws) for g, ws in allp.items()}
    reduced = {g: [w for w in ws if len(w) == lengths[g]] for g, ws in allp.items()}
    return lengths, reduced


def subword_product(gens, word, mask):
    out = ident(len(gens[0]))
    for j, i in enumerate(word):
        if mask & (1 << j):
            out = mm(out, gens[i])
    return out


def contains_subword(gens, word, target):
    return any(subword_product(gens, word, mask) == target for mask in range(1 << len(word)))


def bruhat_leq(gens, reduced, u, v):
    return any(contains_subword(gens, q, u) for q in reduced[v])


def check_s4_fidelity():
    # S_4 simple transpositions, represented as maps and multiplied in display order.
    def trans(i):
        p = list(range(4))
        p[i], p[i + 1] = p[i + 1], p[i]
        return tuple(p)

    def pmul(p, q):
        return tuple(p[q[i]] for i in range(4))

    def pword(gens, word):
        out = tuple(range(4))
        for i in word:
            out = pmul(out, gens[i])
        return out

    gens = [trans(0), trans(1), trans(2)]
    q = (0, 2, 1, 0, 2)  # (s_1,s_3,s_2,s_1,s_3), one-based paper notation
    v = gens[1]
    w = pword(gens, q)
    lengths = {}
    for n in range(7):
        for z in product(range(3), repeat=n):
            g = pword(gens, z)
            lengths[g] = min(lengths.get(g, n), n)
    assert lengths[w] == len(q), (lengths[w], q)
    hits = [mask for mask in range(1 << len(q)) if pword(gens, tuple(q[j] for j in range(len(q)) if mask & (1 << j))) == v]
    assert len(hits) == 1, hits
    print("S4 fidelity: Elek's Example magic word has one positional subword for s_2")


def check_b2():
    e = ((1, 0), (0, 1))
    s = ((-1, 0), (0, 1))
    t = ((0, 1), (1, 0))
    gens = [s, t]
    lengths, reduced = lengths_and_words(gens, 8)
    w = mm(mm(s, t), s)
    assert lengths[w] == 3
    assert reduced[w] == [(0, 1, 0)], reduced[w]
    # The displayed reduced expression has exactly these lower-interval elements.
    sub = {subword_product(gens, (0, 1, 0), m) for m in range(8)}
    assert sub == {e, s, t, mm(s, t), mm(t, s), w}
    lower = sub
    interval = {u for u in lower if bruhat_leq(gens, reduced, s, u)}
    st, ts = mm(s, t), mm(t, s)
    assert interval == {s, st, ts, w}
    # The four-element interval is Boolean: s is bottom, st and ts are incomparable, w top.
    assert bruhat_leq(gens, reduced, s, st) and bruhat_leq(gens, reduced, s, ts)
    assert not bruhat_leq(gens, reduced, st, ts) and not bruhat_leq(gens, reduced, ts, st)
    assert bruhat_leq(gens, reduced, st, w) and bruhat_leq(gens, reduced, ts, w)
    hits = [m for m in range(8) if subword_product(gens, (0, 1, 0), m) == s]
    assert hits == [1, 4], hits  # 001 and 100, positional masks
    print("B2 counterexample: unique reduced word, Boolean interval of size 4, magic count = 2")


def dihedral_mul(x, y, m):
    # x,y are (kind, index), where kind 0 is r_i and kind 1 is s r_i.
    k, i = x
    l, j = y
    if k == 0 and l == 0:
        return (0, (i + j) % m)
    if k == 0 and l == 1:
        return (1, (j - i) % m)
    if k == 1 and l == 0:
        return (1, (i + j) % m)
    return (0, (j - i) % m)


def dihedral_word(word, m):
    out = (0, 0)
    gens = [(1, 0), (1, 1)]  # s=sr_0, t=sr_1
    for i in word:
        out = dihedral_mul(out, gens[i], m)
    return out


def check_family():
    for m in range(4, 101):
        s, t = (1, 0), (1, 1)
        w = dihedral_word((0, 1, 0), m)
        words3 = {dihedral_word(z, m): z for z in product((0, 1), repeat=3)}
        assert words3[w] == (0, 1, 0), (m, w, words3[w])
        short = {dihedral_word(z, m) for n in range(3) for z in product((0, 1), repeat=n)}
        assert w not in short, (m, w)
        hits = [mask for mask in range(8) if dihedral_word(tuple((0, 1, 0)[j] for j in range(3) if mask & (1 << j)), m) == s]
        assert hits == [1, 4], (m, hits)
    print("I_2(m) family: exact normal-form obstruction verified for every m = 4,...,100")


if __name__ == "__main__":
    check_s4_fidelity()
    check_b2()
    check_family()
    print("ALL_CHECKS_OK")
