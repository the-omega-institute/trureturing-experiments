#!/usr/bin/env python3
"""Exact checks for Merhav arXiv:2606.26748v1, Section 6.

No floating-point arithmetic is used. r**U denotes exp(beta*U), beta=log(r).
The matrix checks cover every real local kernel at the listed finite sizes;
they do not prove the unrestricted measurable-space conjecture.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product


def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def transpose(a):
    return [list(x) for x in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum((x*y for x, y in zip(row, col)), F(0))
             for col in bt] for row in a]


def quad(v, a):
    return sum((v[i]*a[i][j]*v[j] for i in range(len(v))
                for j in range(len(v))), F(0))


def finite_model(energies, r, k, d):
    assert k > 0 and d > 0 and k % d == 0
    q = k // d
    alphabet = range(len(energies))
    z = sum((F(1, r**u) for u in energies), F(0))
    single = [F(1, r**u)/z for u in energies]
    chunks = list(product(alphabet, repeat=d))
    hs = [sum(energies[x] for x in a) for a in chunks]
    nu = [product_value(single[x] for x in a) for a in chunks]
    states = list(product(range(len(chunks)), repeat=q))
    lookup = {a: i for i, a in enumerate(states)}
    n = len(states)
    probs = [product_value(nu[j] for j in a) for a in states]
    totals = [sum(hs[j] for j in a) for a in states]
    fibres = defaultdict(list)
    for i, a in enumerate(states):
        fibres[tuple(hs[j] for j in a)].append(i)
    projection = zero(n)
    for ids in fibres.values():
        # The Boltzmann density is constant on each chunk-energy fibre.
        assert len({probs[i] for i in ids}) == 1
        for i in ids:
            for j in ids:
                projection[i][j] = F(1, len(ids))
    assert sum(probs) == 1
    assert matmul(projection, projection) == projection
    for j in range(n):
        assert sum(probs[i]*projection[i][j] for i in range(n)) == probs[j]
        # lambda is counting measure: this is M-preservation for every kernel.
        assert sum(F(1, r**totals[i])*projection[i][j]
                   for i in range(n)) == F(1, r**totals[j])

    # Exact local covariance expression C, so w^T C w = R_0 + 2 sum R_l.
    c = [[-F(2*q-1)*probs[i]*probs[j] for j in range(n)]
         for i in range(n)]
    for i in range(n):
        c[i][i] += probs[i]
    for lag in range(1, q):
        for a in product(range(len(chunks)), repeat=q+lag):
            p = product_value(nu[j] for j in a)
            i, j = lookup[a[:q]], lookup[a[lag:lag+q]]
            c[i][j] += p
            c[j][i] += p

    # Independent global ring calculation: t >= 2q, n_scalar=d*t.
    t = 2*q
    ring_cov = zero(n)
    residual_gram = zero(n)
    conditional_residual = defaultdict(lambda: [F(0)]*n)
    for a in product(range(len(chunks)), repeat=t):
        p = product_value(nu[j] for j in a)
        counts = [0]*n
        for start in range(t):
            counts[lookup[tuple(a[(start+j) % t] for j in range(q))]] += 1
        v = [F(x, t) for x in counts]
        pv = [sum(v[i]*projection[i][j] for i in range(n)) for j in range(n)]
        residual = [v[i]-pv[i] for i in range(n)]
        centre = [v[i]-probs[i] for i in range(n)]
        stat = tuple(hs[j] for j in a)
        for i in range(n):
            conditional_residual[stat][i] += p*residual[i]
            for j in range(n):
                ring_cov[i][j] += t*p*centre[i]*centre[j]
                residual_gram[i][j] += t*p*residual[i]*residual[j]
    assert all(all(x == 0 for x in v) for v in conditional_residual.values())
    assert ring_cov == c  # Checks source equation (35), including cyclic counting.
    projected_c = matmul(transpose(projection), matmul(c, projection))
    difference = [[c[i][j]-projected_c[i][j] for j in range(n)] for i in range(n)]
    assert difference == residual_gram
    # The last equality is a positive weighted Gram factorization, valid for ALL w.
    print(f'ALL_KERNEL_MATRIX energies={energies} exp_beta={r} k={k} d={d} '
          f'kernel_dimension={n} ring_configurations={len(chunks)**t} '
          'normalization=OK global_conditional_projection=OK '
          'covariance_pair_count=OK variance_gap_gram=OK')
    return states, totals, probs, projection, c, chunks, z


def product_value(values):
    answer = F(1)
    for value in values:
        answer *= value
    return answer


def fidelity():
    # Proposition 2 and equations (8)-(10): total-energy averaging at NOL.
    energies, r, k = (0, 1, 1), 2, 3
    xs = list(product(range(len(energies)), repeat=k))
    z = sum((F(1, r**u) for u in energies), F(0))
    u = [sum(energies[x] for x in a) for a in xs]
    m = [F(1+(i+1)**2) for i in range(len(xs))]
    groups = defaultdict(list)
    for i, v in enumerate(u):
        groups[v].append(i)
    bar = [F(0)]*len(xs)
    for ids in groups.values():
        for i in ids:
            bar[i] = sum((m[j] for j in ids), F(0))/len(ids)
    mass = sum(m)
    assert sum(bar) == mass
    p = [F(1, r**v)/z**k for v in u]
    w = [m[i]*r**u[i] for i in range(len(xs))]
    wb = [bar[i]*r**u[i] for i in range(len(xs))]
    mean = sum(p[i]*w[i] for i in range(len(xs)))
    second = sum(p[i]*w[i]**2 for i in range(len(xs)))
    var = second-mean**2
    ratio = z**k*sum(m[i]**2*r**u[i] for i in range(len(xs)))/mass**2
    assert mean == mass/z**k
    assert ratio == 1+var/mean**2
    assert var/(k*mean**2) == (ratio-1)/k
    second_bar = sum(p[i]*wb[i]**2 for i in range(len(xs)))
    assert second_bar < second
    print(f'SOURCE_FIDELITY Proposition_2 M={mass} mean={mean} '
          f'Q_k={ratio} V_nol={(ratio-1)/k} strict_total_energy_improvement=OK')


def worked_example(model):
    states, totals, probs, projection, c, chunks, z = model
    m = [F(1+chunks[a[0]][0]*chunks[a[1]][0]) for a in states]
    w = [m[i]*2**totals[i] for i in range(len(states))]
    wb = [sum(projection[i][j]*w[j] for j in range(len(states)))
          for i in range(len(states))]
    bar = [wb[i]/2**totals[i] for i in range(len(states))]
    for i, a in enumerate(states):
        assert bar[i] == 1+F(sum(chunks[a[0]])*sum(chunks[a[1]]), 4)
    mass = sum(m)
    mu = sum(probs[i]*w[i] for i in range(len(states)))
    r0 = sum(probs[i]*w[i]**2 for i in range(len(states)))-mu**2
    b0 = sum(probs[i]*wb[i]**2 for i in range(len(states)))-mu**2
    total, btotal = quad(w, c), quad(wb, c)
    r1, b1 = (total-r0)/2, (btotal-b0)/2
    v, bv = 2*total/(16*mu**2), 2*btotal/(16*mu**2)
    assert (mass,z,mu) == (20,F(3,2),F(320,81))
    assert (r0,r1,b0,b1) == (F(142544,6561),F(49664,6561),
                              F(128288,6561),F(47360,6561))
    assert (v,bv,v-bv) == (F(15117,51200),F(6969,25600),F(1179,51200))
    assert (total/4,btotal/4) == (F(20156,2187),F(18584,2187))
    print(f'NONDEGENERATE_EXAMPLE M={mass} Z={z} mu={mu} '
          f'R=({r0},{r1}) projected_R=({b0},{b1}) '
          f'V={v} projected_V={bv} improvement={v-bv} '
          f'n8_ring_variances=({total/4},{btotal/4})')


if __name__ == '__main__':
    fidelity()
    example = None
    for energies, r, maximum in [((0, 1), 2, 4), ((0, 0, 1), 3, 3)]:
        for k in range(1, maximum+1):
            for d in range(1, k+1):
                if k % d == 0:
                    model = finite_model(energies, r, k, d)
                    if energies == (0, 1) and k == 4 and d == 2:
                        example = model
    worked_example(example)
    print('ALL_CHECKS_OK')
