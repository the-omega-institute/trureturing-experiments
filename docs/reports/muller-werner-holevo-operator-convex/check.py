#!/usr/bin/env python3
"""Exact finite checks for the Muller Werner--Holevo preregistration.

No floating point is used. Matrix-unit polarization checks the purity
identity on the complete matrix algebra for every listed dimension tuple.
This finite check is evidence about the implementation, not a proof of
the arbitrary-dimension or arbitrary-function conjecture.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod
import sympy as S
from sympy.polys.matrices import DomainMatrix


def indices(ds):
    return list(product(*(range(d) for d in ds)))


def unit_outputs(ds):
    """Tensor channel on E_ab, from Phi_d(E_ab)=(delta_ab I-E_ba)/(d-1)."""
    ix = indices(ds)
    outputs = {}
    for a in ix:
        for b in ix:
            local = []
            for d, ai, bi in zip(ds, a, b):
                entries = {}
                if ai == bi:
                    for j in range(d):
                        entries[j, j] = 1
                entries[bi, ai] = entries.get((bi, ai), 0) - 1
                local.append([(r, c, v) for (r, c), v in entries.items() if v])
            out = {}
            for terms in product(*local):
                row = tuple(t[0] for t in terms)
                col = tuple(t[1] for t in terms)
                out[row, col] = Q(prod(t[2] for t in terms), prod(d - 1 for d in ds))
            outputs[a, b] = out
    return ix, outputs


def trace_pair(a, b):
    return sum((v * b.get((c, r), 0) for (r, c), v in a.items()), Q(0))


def reduced_unit(a, b, subset):
    if any(a[i] != b[i] for i in range(len(a)) if i not in subset):
        return None
    return tuple(a[i] for i in subset), tuple(b[i] for i in subset)


def check_all_bilinear_coefficients(ds):
    ix, out = unit_outputs(ds)
    r = prod(d - 1 for d in ds)
    subsets = [tuple(i for i in range(len(ds)) if mask & (1 << i))
               for mask in range(1 << len(ds))]
    reduced = {(a, b, sub): reduced_unit(a, b, sub)
               for a in ix for b in ix for sub in subsets}
    weights = {sub: prod(ds[i] - 2 for i in range(len(ds)) if i not in sub)
               for sub in subsets}
    count = 0
    for a, b, c, d in product(ix, repeat=4):
        rhs = 0
        for sub in subsets:
            e, f = reduced[a, b, sub], reduced[c, d, sub]
            if e is not None and f is not None and e[0] == f[1] and e[1] == f[0]:
                rhs += weights[sub]
        assert trace_pair(out[a, b], out[c, d]) == Q(rhs, r * r), (ds, a, b, c, d)
        count += 1
    print(f'PURITY_IDENTITY dims={ds} all_matrix_unit_pairs={count}', flush=True)
    return count


def tensor_channel(ds, rho):
    ix, out = unit_outputs(ds)
    lookup = {x: i for i, x in enumerate(ix)}
    result = S.zeros(len(ix))
    for a in ix:
        for b in ix:
            coefficient = rho[lookup[a], lookup[b]]
            if coefficient:
                for (x, y), v in out[a, b].items():
                    result[lookup[x], lookup[y]] += coefficient * S.Rational(v.numerator, v.denominator)
    return result.applyfunc(S.expand)


def density(v):
    v = S.Matrix(v)
    return (v * v.conjugate().T / (v.conjugate().T * v)[0]).applyfunc(S.expand)


def check_source_spectra():
    for d in range(2, 7):
        p = S.zeros(d * d)
        p[0, 0] = 1
        out = tensor_channel((d, d), p)
        r = (d - 1) ** 2
        assert out.eigenvals() == {S.Integer(0): 2 * d - 1, S.Rational(1, r): r}
        v = [int(i // d == i % d) for i in range(d * d)]
        ent = tensor_channel((d, d), density(v))
        high = S.Rational(2 * d - 2, d * r)
        low = S.Rational(d - 2, d * r)
        assert ent.eigenvals() == {high: 1, low: d * d - 1}
        print(f'SOURCE_SPECTRA d={d} product_rank={r} entangled_eigenvalues={high},{low}', flush=True)


def trace_tests(ds, rho):
    d, r = prod(ds), prod(x - 1 for x in ds)
    sigma = tensor_channel(ds, rho)
    assert sigma == sigma.conjugate().T
    assert S.trace(sigma) == 1
    purity = S.expand(S.trace(sigma * sigma))
    assert 0 < purity <= S.Rational(1, r)
    assert purity <= S.Rational(1, r)  # f(x)=x^2
    for t in [S.Rational(1, 7), S.Rational(1, 2), S.Integer(1), S.Integer(3)]:
        # f_t(x)=x^2/(x+t), h_t(x)=x/(x+t).
        val = S.expand(S.trace(sigma * sigma * (sigma + t * S.eye(d)).inv(method="DM")))
        jensen = purity / (purity + t)
        bound = S.Rational(1, r) / (S.Rational(1, r) + t)
        assert val <= jensen <= bound
    for lam in [S.Rational(-4, 5), S.Rational(-1, 2), S.Integer(0)]:
        # Theorem 10's f_lambda(x)=(2x-1)^2/(1-lambda(2x-1)).
        a = 2 * sigma - S.eye(d)
        value = S.expand(S.trace(a * a * (S.eye(d) - lam * a).inv(method="DM")))
        f = lambda x: (2 * x - 1) ** 2 / (1 - lam * (2 * x - 1))
        assert value <= r * f(S.Rational(1, r)) + (d - r) * f(S.Integer(0))
    # Endpoint discontinuity: f(x)=x^2+1_{x=0}; spectral value uses exact nullity.
    nullity = d - DomainMatrix.from_Matrix(sigma).rank()
    assert purity + nullity <= S.Rational(1, r) + d - r
    return sigma, purity


def check_states():
    total = 0
    for ds in [(2, 2), (2, 3), (3, 3), (2, 2, 2), (2, 3, 2), (3, 3, 3), (4, 4)]:
        d, r = prod(ds), prod(x - 1 for x in ds)
        samples = []
        v = [S.Integer(0)] * d
        v[0] = 1
        samples.append(density(v))
        v[-1] = S.I
        samples.append(density(v))
        samples.append(density([1 if i % 3 == 0 else -1 if i % 3 == 1 else S.I for i in range(d)]))
        samples.append(S.eye(d) / d)
        samples.append((samples[0] + 2 * samples[1] + 3 * samples[2]) / 6)
        for sample_index, rho in enumerate(samples):
            print(f'TRACE_CASE dims={ds} sample={sample_index}', flush=True)
            _, purity = trace_tests(ds, rho)
            total += 1
        # Non-coordinate pure product: every local vector has entries 1,i,1,i,... .
        local = [S.Matrix([S.I if j % 2 else 1 for j in range(k)]) for k in ds]
        vector = local[0]
        for a in local[1:]:
            vector = S.kronecker_product(vector, a)
        sigma = tensor_channel(ds, density(vector))
        assert (sigma * sigma).applyfunc(S.expand) == sigma / r
        assert DomainMatrix.from_Matrix(sigma).rank() == r
        _, purity = trace_tests(ds, density(vector))
        assert purity == S.Rational(1, r)
        total += 1
        print(f'TRACE_FUNCTIONS dims={ds} states=6 pure_product_equality=OK', flush=True)
    # All real ternary vectors in dimension four, modulo sign and normalization.
    seen = set()
    for v in product([-1, 0, 1], repeat=4):
        if not any(v):
            continue
        rho = density(v)
        key = tuple(rho)
        if key in seen:
            continue
        seen.add(key)
        trace_tests((2, 2), rho)
        total += 1
    print(f'TRACE_FUNCTIONS exhaustive_real_ternary_two_qubit_projectors={len(seen)} total_states={total}', flush=True)
    return total


def main():
    print(f'Exact arithmetic: Python fractions and SymPy {S.__version__}; no floating point', flush=True)
    check_source_spectra()
    shapes = [(d,) for d in range(2, 5)]
    shapes += list(product(range(2, 4), repeat=2))
    shapes += list(product(range(2, 4), repeat=3))
    pairs = sum(check_all_bilinear_coefficients(ds) for ds in shapes)
    states = check_states()
    print(f'SUMMARY dimension_tuples={len(shapes)} bilinear_coefficients={pairs} density_inputs={states}', flush=True)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
