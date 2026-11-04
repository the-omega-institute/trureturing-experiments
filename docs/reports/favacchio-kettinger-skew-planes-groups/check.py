#!/usr/bin/env python3
"""Exact checks for Favacchio--Kettinger arXiv:2607.11259v1, Conjecture 3.19.

All coefficients are integers, rationals, or Gaussian rationals. The finite
checks support the general repeated-block argument described in README.md;
they do not enumerate all dimensions or constitute a Lean settlement.
"""
from collections import deque
from itertools import combinations, permutations
import sympy as s

I = s.I
Q = s.Rational


def exact(z):
    return s.expand(z)


def normalized(v):
    v = tuple(map(exact, v))
    pivot = next(z for z in v if z != 0)
    return tuple(s.cancel(z / pivot).expand() for z in v)


def cls(m):
    return normalized(m)


def equal(m, n):
    return all(exact(z) == 0 for z in m - n)


def group(generators):
    d = generators[0].rows
    identity = s.ImmutableMatrix(s.eye(d))
    seen = {cls(identity): identity}
    queue = deque([identity])
    while queue:
        h = queue.popleft()
        for g in generators:
            product = s.ImmutableMatrix((g * h).applyfunc(exact))
            key = cls(product)
            if key not in seen:
                # Store normalized representatives, so repeated multiplication
                # cannot accumulate arbitrary scalar factors.
                seen[key] = s.ImmutableMatrix(d, d, key)
                queue.append(seen[key])
                assert len(seen) <= 120, 'unexpected projective closure growth'
    return seen


def delta(m, k):
    return s.ImmutableMatrix(s.diag(*([m] * k)))


def monomials(d, degree):
    if d == 1:
        yield (degree,)
    else:
        for first in range(degree + 1):
            for rest in monomials(d - 1, degree - first):
                yield (first,) + rest


def source_fidelity():
    # Source Section 2, example labelled ex:12-points-P5: M=diag(2,3,5).
    m, eye, zero = s.diag(2, 3, 5), s.eye(3), s.zeros(3)
    assert m.det() == 30 and (m - eye).det() == 8
    assert m.charpoly().as_expr().factor() == (
        (m.charpoly().gen - 2) * (m.charpoly().gen - 3) * (m.charpoly().gen - 5)
    ).expand().factor()
    frames = [zero.col_join(eye), eye.col_join(zero), eye.col_join(eye), eye.col_join(m)]
    points = [normalized(frame[:, j]) for frame in frames for j in range(3)]
    assert len(set(points)) == 12
    hilbert = []
    for degree in range(5):
        powers = list(monomials(6, degree))
        evaluation = s.Matrix([
            [s.prod(z ** exponent for z, exponent in zip(p, a)) for a in powers]
            for p in points
        ])
        hilbert.append(evaluation.rank())
    h_vector = [hilbert[0]] + [b - a for a, b in zip(hilbert, hilbert[1:])]
    assert hilbert == [1, 6, 9, 12, 12]
    assert h_vector == [1, 5, 3, 3, 0]
    print('FIDELITY source ex:12-points-P5: 12 points, Hilbert=[1,6,9,12,12], h=(1,5,3,3)')


A = s.ImmutableMatrix([[Q(1, 2) + I / 2, 0], [0, Q(1, 2) - I / 2]])
B = s.ImmutableMatrix([[Q(1, 2), I / 2], [I / 2, Q(1, 2)]])
S = [s.ImmutableMatrix(s.zeros(2)), s.ImmutableMatrix(s.eye(2)), A, B]
D = A - B


def primitive(matrices, i, j, k):
    """Source f_ijk, in graph coordinates; index -1 is infinity."""
    d = matrices[0].rows
    if k == -1:
        return s.ImmutableMatrix(s.eye(d))
    if j == -1:
        return matrices[i] - matrices[k]
    if i == -1:
        return (matrices[j] - matrices[k]).inv()
    return (matrices[j] - matrices[k]).inv() * (matrices[i] - matrices[k])


def geometry_fidelity(k):
    """Independently solve the ambient incidence equations for every projection.

    Frame_j * q - Frame_k * t = Frame_i * u gives the intersection of
    span([Frame_i*u], L_k) with L_j, without assuming the source formula.
    """
    matrices = [delta(m, k) for m in S]
    d = 2 * k
    eye, zero = s.eye(d), s.zeros(d)
    frames = {-1: zero.col_join(eye)}
    frames.update({i: eye.col_join(m) for i, m in enumerate(matrices)})
    count = 0
    for i, j, center in permutations(frames, 3):
        solution = frames[j].row_join(-frames[center]).inv() * frames[i]
        t = solution[:d, :]
        expected = primitive(matrices, i, j, center)
        assert equal(t, expected)
        assert t.det() != 0
        assert equal(primitive(matrices, j, i, center) * t, eye)
        count += 1
    assert count == 60
    print(f'GEOMETRY k={k}: all 60 source primitive-projection formulas and inverses checked by ambient incidence')


def permutation_closure(generators):
    identity = tuple(range(6))
    seen, queue = {identity}, deque([identity])
    while queue:
        h = queue.popleft()
        for g in generators:
            product = tuple(g[h[j]] for j in range(6))
            if product not in seen:
                seen.add(product)
                queue.append(product)
    return seen


def main():
    source_fidelity()
    dets = [(S[j] - S[i]).det().expand() for i, j in combinations(range(4), 2)]
    assert dets == [1, Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2)]
    assert equal(A * (A - s.eye(2)), -s.eye(2) / 2)
    assert equal(B * (B - s.eye(2)), -s.eye(2) / 2)
    print('SKEWNESS six determinants:', dets)
    print('IDENTITIES A(A-I)=B(B-I)=-I/2')
    H = group([A, B, D])
    assert len(H) == 24
    assert cls(A * B) != cls(B * A)
    X = [(1, 0), (0, 1), (1, 1), (-1, 1), (I, 1), (-I, 1)]
    keys = [normalized(x) for x in X]

    def action(m):
        result = tuple(keys.index(normalized(m * s.Matrix(x))) for x in X)
        assert sorted(result) == list(range(6))
        return result

    action_H = {action(m) for m in H.values()}
    independently_generated = permutation_closure(list(map(action, [A, B, D])))
    assert action_H == independently_generated and len(action_H) == 24
    assert len({p[0] for p in action_H}) == 6
    assert sum(p[0] == 0 for p in action_H) == 4
    assert normalized(A * B * s.Matrix([0, 1])) == normalized((-1, 1))
    assert normalized(B * A * s.Matrix([0, 1])) == normalized((I, 1))
    print('GROUP matrix closure=24; independent permutation closure=24; faithful action; orbit=6; stabilizer=4; noncommuting A,B')
    geometry_fidelity(1)
    geometry_fidelity(2)
    geometry_fidelity(3)
    indices = [-1, 0, 1, 2, 3]
    edges = {(i, j, k): s.ImmutableMatrix(primitive(S, i, j, k).applyfunc(exact))
             for i, j, k in permutations(indices, 3)}
    assert all(cls(m) in H for m in edges.values())
    # Enumerate the complete connected groupoid in graph-coordinate classes,
    # retaining source/target indices. This is stronger than closing A and B.
    identity_key = cls(s.eye(2))
    states = {(-1, identity_key)}
    queue = deque([(-1, s.eye(2))])
    while queue:
        i, h = queue.popleft()
        for (source, target, center), t in edges.items():
            if source != i:
                continue
            key = cls(t * h)
            state = (target, key)
            if state not in states:
                states.add(state)
                queue.append((target, s.ImmutableMatrix(2, 2, key)))
                assert len(states) <= 120
    isotropy = {key for target, key in states if target == -1}
    assert isotropy == set(H) and len(states) == 120
    print('FULL GROUPOID from infinity: 120 morphisms, 24 in each target fiber; isotropy equals difference group')
    for k in range(2, 7):
        matrices = [delta(m, k) for m in S]
        assert all(exact((matrices[j] - matrices[i]).det()) == dets[t] ** k
                   for t, (i, j) in enumerate(combinations(range(4), 2)))
        assert len({cls(delta(h, k)) for h in H.values()}) == 24
        for g in [A, B, D]:
            for h in H.values():
                assert equal(delta(g * h, k), delta(g, k) * delta(h, k))
                assert equal(delta(g.inv(), k), delta(g, k).inv())
        assert all(equal(primitive(matrices, i, j, center), delta(t, k))
                   for (i, j, center), t in edges.items())
        print(f'AMPLIFICATION k={k}, n={2*k-1}: skewness; 24 distinct projective blocks; primitive formulas preserved')
    v = s.Matrix([1, 0, 0, 1])
    orbit = [delta(h, 2) * v for h in H.values()]
    assert len({normalized(x) for x in orbit}) == 24
    determinant = s.Matrix.hstack(v, delta(A, 2) * v, delta(B, 2) * v, delta(A * B, 2) * v).det().expand()
    assert determinant == -Q(1, 4)
    assert s.Matrix.hstack(*orbit).rank() == 4
    frames = [s.zeros(4).col_join(s.eye(4))]
    frames += [s.eye(4).col_join(delta(m, 2)) for m in S]
    points = [frame * u for frame in frames for u in orbit]
    assert len({normalized(x) for x in points}) == 120
    assert s.Matrix.hstack(*points).rank() == 8
    print('NONDEGENERACY k=2: free orbit=24; determinant=-1/4; 120 distinct points; ambient vector rank=8 (P^7)')
    print('SCOPE finite checks k=2..6; all-k family requires the stated symbolic block-homomorphism argument')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
