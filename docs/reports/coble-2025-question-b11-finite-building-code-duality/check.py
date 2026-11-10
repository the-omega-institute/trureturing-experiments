#!/usr/bin/env python3
"""Exact GF(2) checks for Coble's dissertation, Appendix B, Question B.11.

No external dependencies. Vectors are integer bit masks; row operations use XOR.
All projective coordinates and incidence tests are integers modulo a prime.
The finite test list is evidence about these buildings, not a universal proof.
"""
from collections import deque
from itertools import combinations, permutations, product


def mask(xs):
    return sum(1 << x for x in set(xs))


def members(v):
    while v:
        b = v & -v
        yield b.bit_length() - 1
        v ^= b


def basis(vectors):
    pivots = {}
    for v in vectors:
        while v:
            p = v.bit_length() - 1
            if p in pivots:
                v ^= pivots[p]
            else:
                pivots[p] = v
                break
    return pivots


def contains(pivots, v):
    while v:
        p = v.bit_length() - 1
        if p not in pivots:
            return False
        v ^= pivots[p]
    return True


def components(n, panels, types, support=None):
    support = (1 << n) - 1 if support is None else support
    adj = {c: set() for c in members(support)}
    for s in members(types):
        for panel in panels[s]:
            part = set(members(panel & support))
            for c in part:
                adj[c].update(part - {c})
    result, unseen = [], set(adj)
    while unseen:
        a = unseen.pop()
        todo, comp = [a], {a}
        while todo:
            for b in adj[todo.pop()] & unseen:
                unseen.remove(b)
                comp.add(b)
                todo.append(b)
        result.append(mask(comp))
    return sorted(result)


def analyze(name, m, n, panels, apartments):
    full = (1 << n) - 1
    assert n > 0 and len(panels) == m
    for family in panels:
        assert sum(p.bit_count() for p in family) == n
        assert all(p.bit_count() >= 2 for p in family)
        assert sum(family) == full  # disjoint nonempty panel supports
    residues = {j: components(n, panels, j) for j in range(1 << m)}
    subapts = {j: set() for j in range(1 << m)}
    for a in apartments:
        for family in panels:
            assert all((a & p).bit_count() in (0, 2) for p in family)
        for j in subapts:
            subapts[j].update(components(n, panels, j, a))
    assert subapts[0] == {1 << c for c in range(n)}
    assert residues[(1 << m) - 1] == [full]

    # Graph distances from chamber 0; the unique nearest panel chamber is its gate.
    dist, todo = {0: 0}, deque([0])
    while todo:
        c = todo.popleft()
        for family in panels:
            panel = next(p for p in family if p >> c & 1)
            for d in members(panel):
                if d not in dist:
                    dist[d] = dist[c] + 1
                    todo.append(d)
    assert len(dist) == n
    descent = [0] * n
    for s, family in enumerate(panels):
        for p in family:
            ds = [dist[c] for c in members(p)]
            lo = min(ds)
            assert ds.count(lo) == 1 and set(ds) == {lo, lo + 1}
            for c in members(p):
                if dist[c] > lo:
                    descent[c] |= 1 << s

    # Construct the proposed gate/apartment triangular basis from raw subapartments.
    chosen = []
    for c in range(n):
        j = descent[c]
        r = next(r for r in residues[j] if r >> c & 1)
        lo = min(dist[d] for d in members(r))
        gates = [d for d in members(r) if dist[d] == lo]
        assert len(gates) == 1
        gate = gates[0]
        candidates = [a for a in subapts[j] if a >> c & 1 and a >> gate & 1]
        assert candidates
        a = min(candidates)
        assert all(d == c or dist[d] < dist[c] for d in members(a))
        chosen.append(a)
    assert len(basis(chosen)) == n

    # Simultaneous panel kernels, then every residue kernel and every threshold.
    for s, family in enumerate(panels):
        vs = [chosen[c] for c in range(n) if descent[c] >> s & 1]
        assert all((v & p).bit_count() % 2 == 0 for v in vs for p in family)
        assert len(basis(vs)) == n - len(family)
    for j, rs in residues.items():
        vs = [chosen[c] for c in range(n) if descent[c] & j]
        assert all((v & r).bit_count() % 2 == 0 for v in vs for r in rs)
        assert len(basis(vs)) == n - len(rs)
    for k in range(m + 1):
        generators = [a for j, aa in subapts.items() if j.bit_count() == k for a in aa]
        threshold = [chosen[c] for c in range(n) if descent[c].bit_count() >= k]
        b1, b2 = basis(generators), basis(threshold)
        assert len(b1) == len(b2)
        assert all(contains(b2, a) for a in generators)
        assert all(contains(b1, a) for a in threshold)

    # Direct source code definitions; no basis formula is used to generate the codes.
    acodes = {-1: []}
    bcodes = {-1: []}
    for r in range(m + 1):
        acodes[r] = [a for j, aa in subapts.items() if j.bit_count() == m - r for a in aa]
        bcodes[r] = [a for j, aa in residues.items() if j.bit_count() == m - r for a in aa]
    adims, bdims = [], []
    for r in range(-1, m + 1):
        a, b = acodes[r], bcodes[m - r - 1]
        # Proposition B.12 fidelity, followed by equality via complementary dimensions.
        assert all((x & y).bit_count() % 2 == 0 for x in a for y in b)
        assert len(basis(a)) + len(basis(b)) == n
        adims.append(len(basis(a)))
        bdims.append(len(basis(bcodes[r])))
    assert adims[0] == bdims[0] == 0 and adims[-1] == bdims[-1] == n
    print(f'{name}: rank={m} chambers={n} panels={[len(x) for x in panels]} '
          f'apartments={len(apartments)} A_dims={adims} B_dims={bdims} '
          f'B12=OK B11=OK triangular_panel_residue_threshold=OK')
    return adims, bdims


def product_building(sizes):
    chambers = list(product(*(range(q) for q in sizes)))
    index = {c: i for i, c in enumerate(chambers)}
    panels = []
    for s in range(len(sizes)):
        groups = {}
        for c in chambers:
            groups.setdefault(c[:s] + c[s + 1:], []).append(index[c])
        panels.append([mask(v) for v in groups.values()])
    apartments = {mask(index[c] for c in product(*pairs))
                  for pairs in product(*(list(combinations(range(q), 2)) for q in sizes))}
    return len(sizes), len(chambers), panels, sorted(apartments)


def projective_a2(q):
    assert q in (2, 3)
    def norm(v):
        a = next(a for a in v if a)
        return tuple((x * pow(a, -1, q)) % q for x in v)
    pts = sorted({norm(v) for v in product(range(q), repeat=3) if any(v)})
    lines = pts  # line normal vectors in dual coordinates
    incident = lambda p, l: sum(a * b for a, b in zip(p, l)) % q == 0
    flags = [(p, l) for p in pts for l in lines if incident(p, l)]
    index = {f: i for i, f in enumerate(flags)}
    panels = [[mask(i for i, f in enumerate(flags) if f[s] == v) for v in pts]
              for s in range(2)]
    def line(p, r):
        a, b, c = p
        d, e, f = r
        return norm(((b*f-c*e) % q, (c*d-a*f) % q, (a*e-b*d) % q))
    apartments = set()
    for triple in combinations(pts, 3):
        p, r, t = triple
        if incident(t, line(p, r)):
            continue
        fs = [(x, line(x, y)) for x, y in permutations(triple, 2)]
        apartments.add(mask(index[f] for f in fs))
    assert len(pts) == q*q + q + 1
    assert len(flags) == (q*q + q + 1)*(q + 1)
    return 2, len(flags), panels, sorted(apartments)


def thin_a(n):
    chambers = list(permutations(range(n + 1)))
    index = {c: i for i, c in enumerate(chambers)}
    panels = []
    for s in range(n):
        family = set()
        for c in chambers:
            d = list(c)
            d[s], d[s + 1] = d[s + 1], d[s]
            family.add(mask([index[c], index[tuple(d)]]))
        panels.append(sorted(family))
    return n, len(chambers), panels, [(1 << len(chambers)) - 1]


def main():
    count = 0
    for sizes in [(), (2,), (3,), (4,), (7,), (2, 2), (3, 3), (2, 3, 4), (2, 3, 5)]:
        a, b = analyze(f'product_A1{sizes}', *product_building(sizes))
        if sizes == (3, 3):
            assert a == [0, 4, 8, 9] and b == [0, 1, 5, 9]
        count += 1
    for n in (2, 3):
        a, b = analyze(f'thin_A{n}', *thin_a(n))
        assert a == b  # Proposition B.10 and Coble--Barg Theorem 1.7 sanity check
        count += 1
    for q in (2, 3):
        a, b = analyze(f'GL3_F{q}_flags', *projective_a2(q))
        assert a == [0, q**3, (q*q+q+1)*(q+1)-1, (q*q+q+1)*(q+1)]
        count += 1
    print(f'Checked {count} explicitly constructed buildings, every order -1 <= r <= rank.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
