#!/usr/bin/env python3
"""Exact certificates for Baek--Vargas arXiv:2609.01010v2, Remark 3.9.

Original check program. All coefficients are integers; the final certificate
identity is multiplied by four to represent the half-integral square completion.
This is experimental evidence, not a Lean proof or the proposed proof method.
"""
from collections import Counter, deque
from functools import cache
from itertools import combinations, product


def gamma(a, b, c):
    return (a+b+c-a*b-a*c-b*c
            -(c if a == b == 0 else 0)
            -(b if a == c == 0 else 0)
            -(a if b == c == 0 else 0))


def add_square(poly, variable, linear):
    """Add x_variable * (sum linear[j] x_j)^2 exactly."""
    nz = [(j, c) for j, c in enumerate(linear) if c]
    for j, a in nz:
        poly[tuple(sorted((variable, j, j)))] += a*a
    for (j, a), (k, b) in combinations(nz, 2):
        poly[tuple(sorted((variable, j, k)))] += 2*a*b


def lhs(A, scale=1):
    out = Counter()
    n = len(A)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                out[tuple(sorted((i, j, k)))] += scale*A[j][k]
    return +out, out  # the signed object is the second component


def clean(poly):
    return {m: c for m, c in poly.items() if c}


@cache
def small_decomposition(A):
    """Return integer component sign vectors and entrywise nonnegative N."""
    n = len(A)
    assert n <= 4
    signs = [0]*n
    vectors = []
    for root in range(n):
        if signs[root]:
            continue
        signs[root] = 1
        q = deque([root])
        component = []
        while q:
            i = q.popleft()
            component.append(i)
            for j in range(n):
                if i != j and A[i][j] == -1:
                    if signs[j]:
                        assert signs[j] == -signs[i]
                    else:
                        signs[j] = -signs[i]
                        q.append(j)
        vectors.append(tuple(signs[j] if j in component else 0 for j in range(n)))
    N = tuple(tuple(A[i][j]-sum(s[i]*s[j] for s in vectors)
                    for j in range(n)) for i in range(n))
    assert all(c >= 0 for row in N for c in row)
    return tuple(vectors), N


def base_parts(A):
    """Row squares plus auxiliary weighted squares plus nonnegative cubic."""
    n = len(A)
    aux = []
    nonneg = Counter()
    for i in range(n):
        Z = [j for j in range(n) if j != i and A[i][j] == 0]
        assert len(Z) <= 4
        B = tuple(tuple(A[j][k] for k in Z) for j in Z)
        vectors, N = small_decomposition(B)
        for s in vectors:
            v = [0]*n
            for j, c in zip(Z, s):
                v[j] = c
            aux.append((i, v))
        for j, vj in enumerate(Z):
            for k, vk in enumerate(Z):
                nonneg[tuple(sorted((i, vj, vk)))] += N[j][k]
    for i, j, k in combinations(range(n), 3):
        g = gamma(A[i][j], A[i][k], A[j][k])
        assert g >= 0
        nonneg[(i, j, k)] += 2*g
    assert all(c >= 0 for c in nonneg.values())
    return aux, nonneg


def check_certificate(A):
    n = len(A)
    isolated = next((i for i in range(n)
                     if all(A[i][j] == 0 for j in range(n) if j != i)), None)
    extend = n == 6 and isolated is not None
    keep = [i for i in range(n) if not extend or i != isolated]
    B = tuple(tuple(A[i][j] for j in keep) for i in keep)
    aux, nonneg = base_parts(B)
    out = Counter()
    for i, row in enumerate(B):
        v = [0]*n
        for j, c in enumerate(row):
            v[keep[j]] = 2*c
        if extend:
            v[isolated] = 1
        add_square(out, keep[i], v)
    for i, v in aux:
        w = [0]*n
        for j, c in enumerate(v):
            w[keep[j]] = 2*c
        add_square(out, keep[i], w)
    for m, c in nonneg.items():
        out[tuple(sorted(keep[j] for j in m))] += 4*c
    if extend:
        for i in keep:
            out[tuple(sorted((i, isolated, isolated)))] += 3
        out[(isolated, isolated, isolated)] += 4
    assert clean(out) == clean(lhs(A, 4)[1]), A
    return extend


def hp_matrices(n):
    """All labeled HP matrices, by negative graph and free 0/+1 pairs.

    Every negative graph is enumerated once. Triangles are excluded. Pairs
    joined by two negative edges are forced +1; remaining nonedges are free.
    Thus this is exhaustive without relying on graph-isomorphism software.
    """
    pairs = list(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        adj = [0]*n
        for e, (i, j) in enumerate(pairs):
            if mask >> e & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
        if any((mask >> e & 1) and (adj[i] & adj[j])
               for e, (i, j) in enumerate(pairs)):
            continue
        A = [[int(i == j) for j in range(n)] for i in range(n)]
        free = []
        for e, (i, j) in enumerate(pairs):
            if mask >> e & 1:
                A[i][j] = A[j][i] = -1
            elif adj[i] & adj[j]:
                A[i][j] = A[j][i] = 1
            else:
                free.append((i, j))
        for values in product((0, 1), repeat=len(free)):
            for (i, j), v in zip(free, values):
                A[i][j] = A[j][i] = v
            yield tuple(map(tuple, A))


def fidelity():
    expected = {(1,1,1):0, (0,1,1):1, (0,0,1):0, (0,0,0):0,
                (-1,1,1):2, (-1,0,1):1, (-1,0,0):0, (-1,-1,1):0}
    for triple, value in expected.items():
        assert gamma(*triple) == value
    # Theorem 3.4's identity holds even without HP: all ternary triples.
    for a, b, c in product((-1,0,1), repeat=3):
        A = ((1,a,b),(a,1,c),(b,c,1))
        out = Counter()
        for i in range(3):
            add_square(out, i, A[i])
            Z = [j for j in range(3) if j != i and A[i][j] == 0]
            for j in Z:
                for k in Z:
                    out[tuple(sorted((i,j,k)))] += A[j][k]
        out[(0,1,2)] += 2*gamma(a,b,c)
        assert clean(out) == clean(lhs(A)[1])
    print('FIDELITY: Lemma 3.3 gamma table: 8/8; Theorem 3.4 ternary triples: 27/27')
    # Source Example 3.11: the seven-cycle HP matrix.
    A = tuple(tuple(1 if i == j or (i-j)%7 in (2,5)
                    else -1 if (i-j)%7 in (1,6) else 0
                    for j in range(7)) for i in range(7))
    check_certificate(A)
    print('FIDELITY: Example 3.11 standard order-seven matrix: exact Q1 certificate OK')


def main():
    fidelity()
    for n in range(1,7):
        count = extensions = 0
        for A in hp_matrices(n):
            extensions += check_certificate(A)
            count += 1
        print(f'EXHAUSTIVE: order={n} labeled_HP_matrices={count} square_completion_branch={extensions}',
              flush=True)
    print(f'SMALL_BLOCK_DECOMPOSITIONS: {small_decomposition.cache_info().currsize}')
    print('ARITHMETIC: integers only; identities scaled by 4; no floating point')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
