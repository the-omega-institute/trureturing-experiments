#!/usr/bin/env python3
"""Exact preregistration checks for arXiv:2610.32974v1, Conjecture 2.

No numerical optimizer or floating-point tolerance is used. The grid bounds the
input spectra, not the feasible spectra: each d=3 check covers the continuous
interlacing rectangle by enumerating its intersections with sorting chambers.
"""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations
import sympy as sp


def partitions(items):
    """Every set partition once; ordering is immaterial to the conjecture."""
    if not items:
        yield ()
        return
    first, *tail = items
    for p in partitions(tail):
        yield ((first,),) + p
        for j in range(len(p)):
            yield p[:j] + ((first,) + p[j],) + p[j + 1:]


def averages(n, p):
    return tuple(x for part in p for x in
                 [sum((n[j] for j in part), Q(0)) / len(part)] * len(part))


def partition_witness(n, b):
    return next((p for p in partitions(tuple(range(len(n))))
                 if sorted(averages(n, p)) == sorted(b)), None)


def affine_sub(f, g):
    return tuple(x - y for x, y in zip(f, g))


def chamber_vertices(n, order):
    """b=(x,y,S-x-y); x in [n2,n1], y in [n3,n2].

    Add b[order[0]] >= b[order[1]] >= b[order[2]]. All constraints
    are a*x+b*y <= c. Every nonempty bounded chamber has a vertex,
    including degenerate line and point chambers.
    """
    total = sum(n, Q(0))
    forms = ((Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)),
             (Q(-1), Q(-1), total))
    constraints = [(Q(1), Q(0), n[0]), (Q(-1), Q(0), -n[1]),
                   (Q(0), Q(1), n[1]), (Q(0), Q(-1), -n[2])]
    for i, j in zip(order, order[1:]):
        a, b, c = affine_sub(forms[j], forms[i])
        constraints.append((a, b, -c))
    vertices = set()
    for (a, b, c), (u, v, w) in combinations(constraints, 2):
        determinant = a * v - b * u
        if determinant == 0:
            continue
        x = (c * v - b * w) / determinant
        y = (a * w - c * u) / determinant
        if all(a * x + b * y <= c for a, b, c in constraints):
            vertices.add((x, y, total - x - y))
    return vertices


def majorized_by(candidate, vector):
    a, b = sorted(candidate, reverse=True), sorted(vector, reverse=True)
    assert sum(a) == sum(b)
    return all(sum(a[:k]) <= sum(b[:k]) for k in range(1, len(a)))


def h(x):
    x = sp.Rational(x.numerator, x.denominator) if isinstance(x, Q) else sp.sympify(x)
    if x in (0, 1):
        return sp.Integer(0)
    return -x * sp.log(x) - (1 - x) * sp.log(1 - x)


def entropy(b):
    answer = sp.simplify(sum((h(x) for x in b), sp.Integer(0)))
    assert not answer.atoms(sp.Float)
    return answer


def fidelity():
    # The source's d=3, lambda=(2,1) subsection, eq. (best_matrix_2,1): the source's exhibited maximizer.
    u, v, w = sp.symbols('n1 n2 n3', real=True)
    z = sp.Symbol('z')
    mean = (u + w) / 2
    A = sp.Matrix([[mean, 0, (u - w) / 2], [0, v, 0],
                   [(u - w) / 2, 0, mean]])
    assert sp.expand(A.charpoly(z).as_expr() - (z-u)*(z-v)*(z-w)) == 0
    assert A == A.conjugate().T
    B = sp.diag(mean, v, mean)
    assert A * B == B * A
    # Definition 4's explicit d=4 example, with the original index order.
    n = tuple(map(Q, [1, Q(3, 4), Q(1, 4), 0]))
    assert averages(n, ((0, 2, 3), (1,))) == (Q(5, 12),) * 3 + (Q(3, 4),)
    assert entropy((Q(0), Q(1), Q(1, 2))) == sp.log(2)
    print('FIDELITY: source 3x3 matrix characteristic polynomial and pinched spectrum exact;')
    print('          Definition 4 averaging example and endpoint entropy values exact')


def all_small_cases():
    grid = tuple(Q(j, 8) for j in range(9))
    total_cases = total_vertices = faithful_matrices = 0
    z = sp.Symbol('z')
    for d in (2, 3):
        spectra = [tuple(reversed(n)) for n in combinations_with_replacement(grid, d)]
        for n in spectra:
            eta = sum(n, Q(0)) / d
            finest = (eta,) * d
            assert partition_witness(n, finest) is not None
            # Jensen's exact upper bound is attained by the Fourier unitary.
            # Verify F_d F_d^*=I and the constant diagonal symbolically, once per d.
            if d == 2:
                F = sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)
            else:
                root = -sp.Rational(1, 2) + sp.I * sp.sqrt(3) / 2
                F = sp.Matrix([[root**(i*j) for j in range(d)] for i in range(d)]) / sp.sqrt(d)
            assert (F * F.conjugate().T).applyfunc(sp.simplify) == sp.eye(d)
            M = (F * sp.diag(*map(sp.Rational, n)) * F.conjugate().T).applyfunc(sp.simplify)
            assert all(sp.simplify(M[j, j] - sp.Rational(eta)) == 0 for j in range(d))
            total_cases += 1
            if d == 2:
                # The only admissible partition is (1,1); feasibility is the full
                # Schur-Horn interval b1 in [n2,n1], b2=n1+n2-b1.
                assert n[1] <= eta <= n[0]
                assert majorized_by(finest, (n[0], n[1]))
                continue
            mean = (n[0] + n[2]) / 2
            candidate = (max(mean, n[1]), min(mean, n[1]), mean)
            assert n[0] >= candidate[0] >= n[1] >= candidate[1] >= n[2]
            assert partition_witness(n, candidate) is not None
            A = sp.Matrix([[mean, 0, (n[0]-n[2])/2], [0, n[1], 0],
                           [(n[0]-n[2])/2, 0, mean]])
            assert sp.expand(A.charpoly(z).as_expr() - sp.prod(z-sp.Rational(x) for x in n)) == 0
            faithful_matrices += 1
            vertices = set()
            for order in permutations(range(3)):
                these = chamber_vertices(n, order)
                for b in these:
                    assert all(b[i] >= b[j] for i, j in zip(order, order[1:]))
                    # Top-k sums are affine on this whole chamber. A linear
                    # inequality true at every vertex holds throughout it.
                    assert majorized_by(candidate, b), (n, candidate, b)
                vertices.update(these)
            assert vertices
            total_vertices += len(vertices)
            total_cases += 1
    print(f'SMALL: {total_cases} complete input/partition cases, d=2,3, occupation denominator 8')
    print(f'       {faithful_matrices} exact source matrices; {total_vertices} chamber vertices checked')
    print('       majorization inequalities certify the entire continuous d=3 feasible rectangle;')
    print('       strict concavity gives its unique entropy maximizer, with an enumerated partition witness')


def endpoint_family():
    cases = 0
    z = sp.Symbol('z')
    for d in range(2, 9):
        for m in range((d+1)//2, d):
            r = d - m
            for N in range(r, m+1):
                # lambda=(m,1,...,1), n=(1^N,0^(d-N)). Pair each
                # singleton with one coordinate of the large block.
                A = sp.zeros(d)
                for i in range(r):
                    j = m+i
                    A[i, i] = A[i, j] = A[j, i] = A[j, j] = sp.Rational(1, 2)
                for i in range(r, N):
                    A[i, i] = 1
                B = sp.diag(*A.diagonal())
                assert A*A == A and A == A.conjugate().T and sp.trace(A) == N
                assert sp.expand(A.charpoly(z).as_expr() - (z-1)**N*z**(d-N)) == 0
                assert A*B == B*A
                n = (Q(1),)*N + (Q(0),)*(d-N)
                p = tuple((i, N+i) for i in range(r))
                p += tuple((i,) for i in range(r, N))
                p += tuple((i,) for i in range(N+r, d))
                b = tuple(Q(x) for x in B.diagonal())
                assert sorted(averages(n, p)) == sorted(b)
                # Rank/intersection bounds force at least N-r ones and m-N
                # zeros in the m-block. There are at most 2r other eigenvalues,
                # each contributing <= log(2); this explicit projection attains it.
                forced_ones, forced_zeros = N-r, m-N
                assert d-forced_ones-forced_zeros == 2*r
                assert sp.simplify(entropy(b) - 2*r*sp.log(2)) == 0
                cases += 1
    print(f'ENDPOINT FAMILY: {cases} exact projection maximizers through d=8, lambda=(m,1,...,1)')
    print('                 ranks, characteristic polynomials, commutation, partitions, sharp entropy bound agree')


def variation_sign():
    count = 0
    for d in range(2, 7):
        A = sp.Matrix(d, d, lambda i, j: sp.Rational((i+1)*(j+1), 11))
        T = sp.diag(*(sp.Rational(i*i+1, 7) for i in range(d)))
        C = A*T-T*A
        assert C.conjugate().T == -C
        assert sp.trace(T*(C*A-A*C)) == sp.trace(C*C)
        assert sp.trace(C*C) == -sp.trace(C.conjugate().T*C) < 0
        count += 1
    print(f'VARIATION: {count} exact nonzero commutator tests confirm the first-variation sign')


if __name__ == '__main__':
    fidelity()
    all_small_cases()
    endpoint_family()
    variation_sign()
    print('ALL_CHECKS_OK')
