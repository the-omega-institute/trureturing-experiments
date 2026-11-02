#!/usr/bin/env python3
"""Exact prime-field checks of arXiv:2606.28611v1, Question 4.

No floating-point arithmetic or third-party packages are used. Codes are computed
on unordered edges from local parity checks, independently of the arc reduction.
"""

from itertools import product
from collections import Counter


def rref(rows, width, p):
    a = [[x % p for x in row] for row in rows if any(x % p for x in row)]
    pivots = []
    for j in range(width):
        k = next((k for k in range(len(pivots), len(a)) if a[k][j]), None)
        if k is None:
            continue
        i = len(pivots)
        a[i], a[k] = a[k], a[i]
        inv = pow(a[i][j], -1, p)
        a[i] = [(x * inv) % p for x in a[i]]
        for k in range(len(a)):
            if k != i and a[k][j]:
                c = a[k][j]
                a[k] = [(x - c * y) % p for x, y in zip(a[k], a[i])]
        pivots.append(j)
        if len(pivots) == len(a):
            break
    return a[:len(pivots)], pivots


def rank(rows, width, p):
    return len(rref(rows, width, p)[1])


def kernel(rows, width, p):
    a, pivots = rref(rows, width, p)
    result = []
    for j in range(width):
        if j in pivots:
            continue
        x = [0] * width
        x[j] = 1
        for row, k in zip(a, pivots):
            x[k] = -row[j] % p
        result.append(x)
    return result


def kron(a, b, p):
    return [(x * y) % p for x in a for y in b]


def tensor_basis(a, b, p):
    return [kron(x, y, p) for x in a for y in b]


def graph(moduli, generators):
    vertices = list(product(*(range(n) for n in moduli)))
    assert len(set(generators)) == len(generators)
    assert all(any(s) for s in generators)
    assert all(tuple(-x % n for x, n in zip(s, moduli)) in generators
               for s in generators)
    edge_index, arc_edges, arcs = {}, {}, []
    for g in vertices:
        for s in generators:
            h = tuple((x + y) % n for x, y, n in zip(g, s, moduli))
            e = tuple(sorted((g, h)))
            edge_index.setdefault(e, len(edge_index))
            arc_edges[g, s] = edge_index[e]
            arcs.append((g, s))
    edges = list(edge_index)
    return vertices, generators, edges, edge_index, arc_edges, arcs


def cycle(n):
    return graph((n,), [(n - 1,), (1,)])


def circulant(n):
    return graph((n,), [(1,), (n - 1,), (2,), (n - 2,)])


def complete(n):
    return graph((n,), [(i,) for i in range(1, n)])


def checks_on_edges(g, local, p):
    vertices, generators, edges, _, incidence, _ = g
    parity = kernel(local, len(generators), p)
    rows = []
    for v in vertices:
        for h in parity:
            row = [0] * len(edges)
            for s, c in zip(generators, h):
                row[incidence[v, s]] = c
            rows.append(row)
    return rows


def factor(g, local, p):
    vertices, generators, edges, _, incidence, arcs = g
    code = kernel(checks_on_edges(g, local, p), len(edges), p)
    parity = kernel(local, len(generators), p)
    arc_index = {a: i for i, a in enumerate(arcs)}
    reversal = []
    for v, s in arcs:
        # Determine the other endpoint from the unordered edge, without a formula
        # for the group law, then find the inverse generator by its index.
        e = edges[incidence[v, s]]
        w = e[1] if e[0] == v else e[0]
        reverse_s = next(t for t in generators if incidence[w, t] == incidence[v, s])
        reversal.append(arc_index[w, reverse_s])
    assert all(reversal[reversal[i]] == i and reversal[i] != i for i in range(len(arcs)))
    uchecks = []
    for v in vertices:
        for h in parity:
            row = [0] * len(arcs)
            for s, c in zip(generators, h):
                row[arc_index[v, s]] = c
            uchecks.append(row)
    reverse_checks = [[row[reversal[i]] for i in range(len(arcs))] for row in uchecks]
    d = kernel(uchecks + reverse_checks, len(arcs), p)
    nimages = [[(x[reversal[i]] - x[i]) % p for i in range(len(arcs))] for x in d]
    r = rank(nimages, len(arcs), p)
    assert len(d) == len(code) + r
    return code, r


def compare(g, a, h, b, p):
    cg, rg = factor(g, a, p)
    ch, rh = factor(h, b, p)
    # This constructor uses endpoint pairs of actual factor edges. It is
    # independent of the arc-space calculation used to compute r_G and r_H.
    vertices = [(x, y) for x in g[0] for y in h[0]]
    generators = [(s, t) for s in g[1] for t in h[1]]
    edge_index, incidence, psi = {}, {}, []
    for x, y in vertices:
        for s, t in generators:
            eg, eh = g[2][g[4][x, s]], h[2][h[4][y, t]]
            xx = eg[1] if eg[0] == x else eg[0]
            yy = eh[1] if eh[0] == y else eh[0]
            e = tuple(sorted(((x, y), (xx, yy))))
            if e not in edge_index:
                edge_index[e] = len(edge_index)
                psi.append((g[4][x, s], h[4][y, t]))
            incidence[(x, y), (s, t)] = edge_index[e]
    delta = (vertices, generators, list(edge_index), edge_index, incidence, [])
    assert set(Counter(psi).values()) == {2}  # Source Lemma 7.4(c).
    assert len(psi) == 2 * len(g[2]) * len(h[2])
    ab = tensor_basis(a, b, p)
    rows = checks_on_edges(delta, ab, p)
    dim_delta = len(psi) - rank(rows, len(psi), p)
    pullbacks = [[x[i] * y[j] % p for i, j in psi] for x in cg for y in ch]
    # Source Proposition 7.5(b): check every tensor basis generator against
    # every defining constraint, not just an inequality of dimensions.
    assert all(sum(c * v for c, v in zip(row, word)) % p == 0
               for word in pullbacks for row in rows)
    dim_p = rank(pullbacks, len(psi), p)
    assert dim_p == len(cg) * len(ch)
    assert dim_delta == dim_p + rg * rh, (p, len(g[0]), len(h[0]), dim_delta, dim_p, rg, rh)
    assert (dim_delta == dim_p) == (rg == 0 or rh == 0)
    return dim_delta, dim_p, len(cg), len(ch), rg, rh


def subspaces2(p):
    return [('zero', []), ('top', [[1, 0], [0, 1]])] + [
        (str(a), [[1, a]]) for a in range(p)] + [('infinity', [[0, 1]])]


def augmentation(d, p):
    return kernel([[1] * d], d, p)


def main():
    print('Exact arithmetic over prime fields F_2, F_3, F_5; no floating point.')
    fidelity = 0
    for p in (2, 3, 5):
        for n in range(3, 9):
            for label, local in subspaces2(p):
                code, _ = factor(cycle(n), local, p)
                if label == 'zero' or label == 'infinity':
                    expected = 0
                elif label == 'top':
                    expected = n
                else:
                    expected = int(pow(int(label), n, p) == 1)
                assert len(code) == expected, (p, n, label)
                fidelity += 1
    print(f'FIDELITY: {fidelity} checks reproduce Proposition 6.2 and the zero/full cases.')
    count = 0
    for p in (2, 3, 5):
        for n, m in product((3, 4), repeat=2):
            for _, a in subspaces2(p):
                for _, b in subspaces2(p):
                    compare(cycle(n), a, cycle(m), b, p)
                    count += 1
    print(f'EXHAUSTIVE: {count} comparisons, every local subspace on C3/C4 ordered pairs.')
    for p in (2, 3, 5):
        spaces = subspaces2(p)
        for k in range(30):
            n, m = 3 + k % 6, 3 + (k * 5 + 1) % 6
            a = spaces[(k * 7 + 2) % len(spaces)][1]
            b = spaces[(k * 11 + 3) % len(spaces)][1]
            compare(cycle(n), a, cycle(m), b, p)
            count += 1
    print('DETERMINISTIC: 90 further cycle comparisons, orders 3 through 8.')
    for p in (3, 5):
        for n, m in ((4, 4), (4, 5), (5, 5)):
            result = compare(complete(n), augmentation(n - 1, p),
                             complete(m), augmentation(m - 1, p), p)
            print(f'COMPLETE: p={p}, n={n}, m={m}, delta={result[0]}, P={result[1]}')
            count += 1
        for n, m in ((5, 6), (6, 7)):
            result = compare(circulant(n), augmentation(4, p),
                             circulant(m), augmentation(4, p), p)
            assert result == (n*m+(n+1)*(m+1), n*m, n, m, n+1, m+1)
            print(f'STRICT_FAMILY: p={p}, n={n}, m={m}, delta={result[0]}, P={result[1]}')
            count += 1
    # Nondegenerate equality-family representatives: an odd repetition cycle
    # paired with a four-regular augmentation circulant over an odd field.
    for p in (3, 5):
        for n in (3, 5, 7):
            result = compare(cycle(n), [[1, 1]], circulant(5), augmentation(4, p), p)
            assert result == (5, 5, 1, 5, 0, 6)
    print('EQUALITY_FAMILY: 6 nondegenerate odd-cycle/augmentation comparisons.')
    print(f'PRODUCT_COMPARISONS: {count} systematic + 6 equality = {count+6}; failures=0.')
    print('Each comparison checks Lemma 7.4, Proposition 7.5(b), the defect identity, and its zero criterion.')
    print('Finite checks support the general argument; they are not a universal or Lean proof.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
