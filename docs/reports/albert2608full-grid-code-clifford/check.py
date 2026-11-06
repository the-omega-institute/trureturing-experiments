#!/usr/bin/env python3
"""Exact binary checks for Albert's bipartite-grid fullness conjecture.

All vectors and matrices are represented by integers over F_2. No third-party
program source or numerical group package is used. The all-parameter assertion
is not proved by this finite experiment.
"""
from array import array
from itertools import combinations
from math import prod


def parity(x):
    return x.bit_count() & 1


class BinaryBasis:
    def __init__(self):
        self.rows = {}

    def reduce(self, x):
        for p in sorted(self.rows, reverse=True):
            if (x >> p) & 1:
                x ^= self.rows[p]
        return x

    def add(self, x):
        x = self.reduce(x)
        if x:
            self.rows[x.bit_length() - 1] = x
            return True
        return False

    def __len__(self):
        return len(self.rows)


def span_rank(vs):
    b = BinaryBasis()
    for v in vs:
        b.add(v)
    return len(b)


def orbit(seeds, dim, start=None):
    """Orbit under transvections, each specified by (center, pairing mask)."""
    seen = bytearray(1 << dim)
    q = array('I', [seeds[0][0] if start is None else start])
    seen[q[0]] = 1
    for v in q:
        for c, f in seeds:
            if parity(v & f):
                w = v ^ c
                if not seen[w]:
                    seen[w] = 1
                    q.append(w)
    return len(q)


def symplectic_order(k):
    return (1 << (k * k)) * prod((1 << (2 * i)) - 1 for i in range(1, k + 1))


def check_grid(a, b, exhaustive=False):
    n = a * b
    k = (a - 2) * (b - 2)
    bit = lambda i, j: 1 << (i * b + j)
    rows = [sum(bit(i, j) for j in range(b)) for i in range(a)]
    cols = [sum(bit(i, j) for i in range(a)) for j in range(b)]
    code = BinaryBasis()
    for r in rows:
        for c in cols:
            code.add(r ^ c)
    cb = list(code.rows.values())
    assert len(code) == a + b - 2
    assert n - 2 * len(code) == k
    assert all(parity(x & y) == 0 for x in cb for y in cb)
    if (a + b) % 4 == 2:
        assert all(x.bit_count() % 4 == 0 for x in cb)
    rectangles = []
    for i, h in combinations(range(a), 2):
        for j, l in combinations(range(b), 2):
            rectangles.append(bit(i, j) ^ bit(i, l) ^ bit(h, j) ^ bit(h, l))
    assert all(parity(w & c) == 0 and code.reduce(w) for w in rectangles for c in cb)

    # The chosen tensor basis of E_a/<1> tensor E_b/<1>, in physical coordinates.
    logical = [bit(i, j) ^ bit(i, b-1) ^ bit(a-1, j) ^ bit(a-1, b-1)
               for i in range(a-2) for j in range(b-2)]
    reduced = {}
    for t, v in enumerate(logical):
        v = code.reduce(v)
        label = 1 << t
        for p in sorted(reduced, reverse=True):
            if (v >> p) & 1:
                v ^= reduced[p][0]
                label ^= reduced[p][1]
        assert v
        reduced[v.bit_length()-1] = (v, label)
    assert len(reduced) == k

    def coord(v):
        v = code.reduce(v)
        label = 0
        for p in sorted(reduced, reverse=True):
            if (v >> p) & 1:
                v ^= reduced[p][0]
                label ^= reduced[p][1]
        assert v == 0
        return label

    gram = [sum(parity(x & y) << j for j, y in enumerate(logical)) for x in logical]
    assert span_rank(gram) == k
    assert all(parity(x & x) == 0 for x in logical)
    assert span_rank(coord(w) for w in rectangles) == k

    # No undetectable vector of weight 1, 2 or 3. A rectangle witnesses d=4.
    syndromes = [sum(parity((1 << q) & c) << j for j, c in enumerate(cb))
                 for q in range(n)]
    assert all(syndromes)
    assert all(syndromes[i] ^ syndromes[j] for i, j in combinations(range(n), 2))
    assert all(syndromes[i] ^ syndromes[j] ^ syndromes[l]
               for i, j, l in combinations(range(n), 3))

    def permute(v, f):
        out = 0
        while v:
            z = v & -v
            q = z.bit_length()-1
            out ^= 1 << f(q)
            v ^= z
        return out

    # Every single row or column transposition preserves the stabilizer code.
    for i, h in combinations(range(a), 2):
        f = lambda q: ((h if q//b == i else i if q//b == h else q//b) * b + q % b)
        assert all(code.reduce(permute(c, f)) == 0 for c in cb)
    for j, l in combinations(range(b), 2):
        f = lambda q: q//b*b + (l if q % b == j else j if q % b == l else q % b)
        assert all(code.reduce(permute(c, f)) == 0 for c in cb)

    # For every pair of transpositions, check I+R+Q+RQ = w w^T on every basis vector.
    identities = 0
    for i, h in combinations(range(a), 2):
        for j, l in combinations(range(b), 2):
            w = bit(i, j) ^ bit(i, l) ^ bit(h, j) ^ bit(h, l)
            for q in range(n):
                r, c = divmod(q, b)
                rr = h if r == i else i if r == h else r
                cc = l if c == j else j if c == l else c
                image = bit(r, c) ^ bit(rr, c) ^ bit(r, cc) ^ bit(rr, cc)
                assert image == (w if ((w >> q) & 1) else 0)
            identities += 1

    # The rectangle pairing-one graph is connected. Its bipartite double is
    # connected because there is a triangle: fix row pair and use 3 column pairs.
    cs = sorted({coord(w) for w in rectangles})

    def pairing_mask(c):
        out = 0
        while c:
            z = c & -c
            out ^= gram[z.bit_length()-1]
            c ^= z
        return out

    seeds = [(c, pairing_mask(c) << k) for c in cs]
    seeds += [(c << k, pairing_mask(c)) for c in cs]
    todo = [0]
    reached = {0}
    for t in todo:
        center, _ = seeds[t]
        for s, (_, f) in enumerate(seeds):
            if s not in reached and parity(center & f):
                reached.add(s)
                todo.append(s)
    assert len(reached) == len(seeds)
    p = coord(bit(0, 0) ^ bit(0, 1) ^ bit(1, 0) ^ bit(1, 1))
    q = coord(bit(0, 1) ^ bit(0, 2) ^ bit(1, 1) ^ bit(1, 2))
    assert p and q and p != q and (p ^ q) in cs
    if exhaustive:
        count = orbit(seeds, 2*k)
        assert count == (1 << (2*k)) - 1
        print(f'FULL_TRANSVECTION_CLOSURE a={a} b={b} centers={count} '
              f'expected={(1 << (2*k))-1} Sp_order={symplectic_order(k)}')
    print(f'GRID a={a} b={b} n={n} k={k} distance=4 code_rank={len(code)} '
          f'logical_gram_rank={k} distinct_rectangles={len(cs)} '
          f'pairing_graph_connected=1 rank_one_identities={identities}')


def check_enlargement_schema():
    # U has dimension 4 and is saturated, v,z is an orthogonal hyperbolic plane.
    # For every nonzero u in U and every w in U, U* ∪ {v+u,z+w}
    # generates all 63 nonzero centers. This is a finite helper check, not a
    # replacement for the all-dimensional saturation proof.
    def mask(c):
        return ((c & 0b010101) << 1) | ((c & 0b101010) >> 1)
    cases = 0
    for u in range(1, 16):
        for w in range(16):
            centers = list(range(1, 16)) + [16 ^ u, 32 ^ w]
            assert orbit([(c, mask(c)) for c in centers], 6) == 63
            cases += 1
    print(f'SATURATION_ENLARGEMENT dim_U=4 cases={cases} closure=63')


def main():
    admissible = [(a, b) for a in range(4, 121, 2) for b in range(a, 121, 2)
                  if a*b <= 120 and (a+b) % 4 == 2]
    assert len(admissible) == 14
    print('SOURCE_FIDELITY arxiv=2608.05688v2 admissible_unordered_grids_ab_le_120=14')
    for a, b in admissible:
        check_grid(a, b, exhaustive=(a, b) == (4, 6))
    check_grid(6, 4, exhaustive=True)
    check_grid(4, 4, exhaustive=True)  # Outside the stated parity range; control.
    check_enlargement_schema()
    print('SCOPE finite checks only; all-parameter fullness remains a proof obligation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
