#!/usr/bin/env python3
"""Exact checks for Yoeli--Gour's fixed-superchannel question (arXiv:2602.18323v1).

Requires Python 3 and SymPy. No floating-point arithmetic or numerical tolerance
is used. These finite checks support the stated dimension-uniform argument;
they do not certify an infinite family or supply a Lean proof.
"""
from itertools import permutations, combinations
import sympy as s
from sympy.polys.matrices import DomainMatrix


def equal(a, b):
    assert a.shape == b.shape
    assert all(s.expand(z) == 0 for z in a - b)


def unit(d, i, j):
    m = s.zeros(d)
    m[i, j] = 1
    return m


def dephase(d):
    return s.diag(*[int(i == j) for i in range(d) for j in range(d)])


def superop(u):
    return s.kronecker_product(u, u.conjugate())


def choi(n, d):
    return s.Matrix(d*d, d*d,
                    lambda p, q: n[(p % d)*d + q % d,
                                   (p // d)*d + q // d])


def vec(u):
    d = u.rows
    return s.Matrix([u[b, a] for a in range(d) for b in range(d)])


def monomial(p, phase_at=None):
    d = len(p)
    return s.Matrix(d, d, lambda b, a:
                    (s.I if a == phase_at else 1) if b == p[a] else 0)


def source_fidelity(d):
    delta = dephase(d)
    equal(delta*delta, delta)
    # Source dephaser: Kraus E_jj, with sum K* K = I.
    kraus = [unit(d, j, j) for j in range(d)]
    equal(sum((k.H*k for k in kraus), s.zeros(d)), s.eye(d))
    equal(sum((superop(k) for k in kraus), s.zeros(d*d)), delta)
    rho = s.eye(d)/s.Integer(d)
    equal(s.Matrix(d*d, 1, list(rho)), delta*s.Matrix(d*d, 1, list(rho)))
    assert s.trace(rho) == 1 and all(rho[j, j] > 0 for j in range(d))
    # Eq. (29): N -> Delta N Delta. Its Choi action is complete dephasing.
    for a in range(d*d):
        for b in range(d*d):
            e = unit(d*d, a, b)
            projected = delta*e*delta
            equal(delta*projected*delta, projected)
            j = choi(e, d)
            equal(choi(projected, d), s.diag(*j.diagonal()))
    # All classical deterministic channels are fixed, with diagonal Choi matrices.
    for f in permutations(range(d)):
        n = s.zeros(d*d)
        for j in range(d):
            n[f[j]*d+f[j], j*d+j] = 1
        equal(delta*n*delta, n)
    # Identity is DIO but is not classical (d >= 2).
    ident = s.eye(d*d)
    equal(ident*delta, delta*ident)
    if d >= 2:
        assert delta*ident*delta != ident
    depol = s.Matrix(d*d, d*d, lambda p, q:
                     s.Rational(1, d) if p//d == p % d and q//d == q % d else 0)
    equal(delta*depol, depol*delta)
    equal(choi(depol, d), s.eye(d*d)/s.Integer(d))
    print(f'FIDELITY d={d}: dephaser CPTP/idempotent/full-rank fixed state; '
          f'Eq29 Choi dephasing; identity DIO; depolarizing Choi=I/{d}')


def family(d):
    ps = list(permutations(range(d)))
    vectors, us, labels = [], [], []
    for p in ps:
        for j in [None] + list(range(d)):
            u = monomial(p, j)
            equal(u.H*u, s.eye(d))
            n = superop(u)
            equal(n*dephase(d), dephase(d)*n)
            v = vec(u)
            equal(choi(n, d), v*v.H)
            equal((v*v.H)*(v*v.H), d*(v*v.H))
            labels.append((p, j)); vectors.append(v); us.append(u)
    # Explicitly produce each vectorized matrix unit as a difference.
    # This establishes the checked family's span without a numerical rank test.
    for a in range(d):
        for b in range(d):
            p = list(range(d)); p[a], p[b] = p[b], p[a]; p = tuple(p)
            difference = vec(monomial(p, a)) - vec(monomial(p))
            target = s.zeros(d*d, 1); target[a*d+b] = s.I-1
            equal(difference, target)
    # Generate the actual nonzero-overlap graph, exactly.
    edges = [set() for _ in vectors]
    for i in range(len(vectors)):
        for j in range(i):
            overlap = s.expand((vectors[i].H*vectors[j])[0])
            if overlap != 0:
                edges[i].add(j); edges[j].add(i)
    remaining = set(range(len(vectors))); sizes = []
    while remaining:
        stack = [remaining.pop()]; count = 0
        while stack:
            i = stack.pop(); count += 1
            for j in edges[i]:
                if j in remaining:
                    remaining.remove(j); stack.append(j)
        sizes.append(count)
    assert len(sizes) == (2 if d == 2 else 1)
    # Check the two overlap formulas used for every dimension in the outline.
    for p in ps:
        for j in range(d):
            assert s.expand((vec(monomial(p)).H*vec(monomial(p, j)))[0]) == d-1+s.I
        for a, b in combinations(range(d), 2):
            q = list(p); q[a], q[b] = q[b], q[a]
            assert (vec(monomial(p)).H*vec(monomial(tuple(q))))[0] == d-2
    # Independent exact joint-commutant rank calculation at d=2,3.
    if d <= 3:
        n = d*d; rows = []
        for v in vectors:
            q = v*v.H
            for a in range(n):
                for b in range(n):
                    row = [s.S.Zero]*(n*n)
                    for k in range(n):
                        row[a*n+k] += q[k, b]
                        row[k*n+b] -= q[a, k]
                    rows.append(row)
        rank = DomainMatrix.from_Matrix(s.Matrix(rows)).to_field().rank()
        dimension = n*n-rank
        assert dimension == (2 if d == 2 else 1)
        print(f'COMMUTANT d={d}: exact rank={rank}, dimension={dimension}')
    # A rational unitary outside DIO; its E00 output has off-diagonal 12/25.
    u = s.eye(d)
    u[:2, :2] = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5)],
                          [s.Rational(4, 5), s.Rational(3, 5)]])
    equal(u.H*u, s.eye(d))
    n = superop(u)
    assert n*dephase(d) != dephase(d)*n
    assert (u*unit(d, 0, 0)*u.H)[0, 1] == s.Rational(12, 25)
    print(f'FAMILY d={d}: {len(vectors)} DIO monomial vectors; span={d*d}; '
          f'overlap component sizes={sizes}; non-DIO rational rotation verified')


def qubit_exception():
    z = s.diag(1, -1); c = superop(z); delta = dephase(2)
    equal(c, 2*delta-s.eye(4))
    # Theta is a physical equal mixture of two unitary pre/post processes.
    # Test the complete linear superoperator on every matrix-unit channel map.
    for a in range(4):
        for b in range(4):
            e = unit(4, a, b)
            t = (e+c*e*c)/2
            equal((t+c*t*c)/2, t)
            assert (t == e) == (e*delta == delta*e)
    w = s.kronecker_product(z, z)
    for a in range(4):
        for b in range(4):
            e = unit(4, a, b)
            equal(choi((e+c*e*c)/2, 2), (choi(e, 2)+w*choi(e, 2)*w.H)/2)
    # Direct fixed-space equality: masks are identical, hence this holds for
    # arbitrary complex linear maps, rather than just sampled CPTP channels.
    theta_mask = s.Matrix(4, 4, lambda a, b: (1+c[a, a]*c[b, b])/2)
    dio_mask = s.Matrix(4, 4, lambda a, b: int(delta[a, a] == delta[b, b]))
    equal(theta_mask, dio_mask)
    print('QUBIT: physical involution twirl idempotent; fixed space exactly DIO')
    # The same channel involution fails CP in dimension 3.
    h = choi(2*dephase(3)-s.eye(9), 3)
    v = vec(s.eye(3))
    assert (v.H*h*v)[0] == -3
    print('BOUNDARY d=3: Choi quadratic form of 2*Delta-id is -3 (not CP)')


def main():
    for d in (1, 2, 3): source_fidelity(d)
    for d in (2, 3, 4): family(d)
    qubit_exception()
    print('Finite exact checks support the general argument; no infinite-family certification claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
