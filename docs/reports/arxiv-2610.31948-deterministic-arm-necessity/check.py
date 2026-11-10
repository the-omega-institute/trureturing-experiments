#!/usr/bin/env python3
"""Exact finite checks for arXiv:2610.31948v1, Theorem 3 and its open question.

Python >= 3.12; SymPy >= 1.14. No floating-point calculations are used.
These checks support a preregistration; they do not prove the arbitrary-law claim.
"""

from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from math import prod
import sympy as sp


def mean(law):
    return sum((z * p for z, p in law.items()), F(0))


def second(law):
    assert sum(law.values()) == 1 and mean(law) == 0
    assert all(p > 0 for p in law.values())
    return sum((z * z * p for z, p in law.items()), F(0))


def expected_max(laws):
    return sum((max(z for z, _ in row) * prod(p for _, p in row)
                for row in product(*(list(law.items()) for law in laws))), F(0))


def binary(x, u):
    return {u: x, -x * u / (1 - x): 1 - x}


def mixture(law):
    """Return (weight, centered kernel) using the positive/negative product formula."""
    b = sum((z * p for z, p in law.items() if z > 0), F(0))
    if not b:
        assert law == {F(0): F(1)}
        return [(F(1), law)]
    out = []
    if law.get(F(0), 0):
        out.append((law[F(0)], {F(0): F(1)}))
    for u, pu in law.items():
        if u <= 0:
            continue
        for z, pv in law.items():
            if z >= 0:
                continue
            v = -z
            out.append(((u + v) * pu * pv / b,
                        {u: v / (u + v), -v: u / (u + v)}))
    reconstructed = {}
    for weight, kernel in out:
        for z, p in kernel.items():
            reconstructed[z] = reconstructed.get(z, F(0)) + weight * p
    assert reconstructed == law
    assert sum(w for w, _ in out) == 1
    assert sum(w * second(k) for w, k in out) == second(law)
    return out


def symbolic_checks():
    t, x, p = sp.symbols('t x P', positive=True)
    q = 1 - x
    g = t / (2 * (t + 1))
    psi = (1 - p**2 - (1 - p)**2 / t) / 2
    lhs = 2 * t * (t + 1) * (g - x * q * (1 - p)**2 - q**2 * psi)
    rhs = (t - (t + 1) * q * (1 - p))**2 + 2*t*(t+1)*x*p*q*(1-p)
    assert sp.factor(lhs - rhs) == 0
    for n in range(1, 6):
        qs = sp.symbols(f'q0:{n}')
        r = [sp.Integer(1)]
        for z in qs:
            r.append(r[-1] * z)
        f = sum((1-qs[j])*qs[j]*r[j]**2 for j in range(n))
        telescope = (1-r[-1]**2-sum((r[j]-r[j+1])**2 for j in range(n)))/2
        assert sp.expand(f - telescope) == 0
    print('SYMBOLIC: source Appendix E SOS identity; telescoping t=1..5: exact')


def cascade_checks():
    for m in range(2, 11):
        laws = [{F(0): F(1)}] + [
            binary(F(1, k+1), F(k, m-1)) for k in range(1, m)]
        e = expected_max(laws)
        v = sum(map(second, laws), F(0))
        assert e == F(1, 2) and v == F(m, 2*(m-1))
        assert e**2 == F(m-1, 2*m) * v
        print(f'FIDELITY M={m}: Emax={e}; V={v}; Emax^2/V={e**2/v}')
        if m <= 5:
            for scale in [F(1, 3), F(1), F(7, 2)]:
                scaled = [{z*scale: p for z, p in law.items()} for law in laws]
                for shift in range(m):
                    moved = scaled[shift:] + scaled[:shift]
                    assert expected_max(moved)**2 == F(m-1, 2*m)*sum(map(second, moved))
    print('EQUALITY: scaled and cyclically permuted source cascades M=2..5: exact')


def grid_checks():
    grid = [binary(x, u) for x in [F(1, 3), F(1, 2), F(2, 3)]
            for u in [F(1, 2), F(1), F(2)]]
    grid.append({F(0): F(1)})
    for m in range(2, 6):
        count = strict = equal = 0
        min_gap = None
        c2 = F(m-1, 2*m)
        for indices in combinations_with_replacement(range(len(grid)), m):
            laws = [grid[i] for i in indices]
            v = sum(map(second, laws), F(0))
            e = expected_max(laws)
            gap = c2*v-e**2
            assert e >= 0 and gap >= 0
            count += 1
            all_nonzero = all(second(law) > 0 for law in laws)
            if all_nonzero:
                assert gap > 0
                strict += 1
                min_gap = gap if min_gap is None else min(min_gap, gap)
                pairs = sorted(((max(law), law[max(law)]) for law in laws), reverse=True)
                xs = [x for _, x in pairs]
                qs = [1-x for x in xs]
                ws = [u*x for u, x in pairs]
                ss = [prod(qs[:i]) for i in range(m)]
                all_down = prod(qs)
                ls = []
                for k in range(m):
                    coeff = [ss[i]-(all_down/qs[k] if i == k else 0) for i in range(m)]
                    lk = sum(ws[i]*coeff[i] for i in range(m))
                    hk = sum(xs[i]*qs[i]*coeff[i]**2 for i in range(m))
                    assert lk**2 <= v*hk
                    if k < m-1:
                        assert hk < c2
                    else:
                        assert coeff[-1] == 0 and hk <= c2 and lk**2 < v*hk
                    ls.append(lk)
                assert e == max(ls)
            if v and gap == 0:
                equal += 1
                assert any(law == {F(0): F(1)} for law in laws)
        print(f'GRID M={m}: tuples={count}; all-nondegenerate strict={strict}; '
              f'positive-V equalities={equal}; min squared strict gap={min_gap}')


def mixture_checks():
    laws = [binary(F(1, 3), F(1)),
            {F(-2): F(1, 4), F(0): F(1, 4), F(1): F(1, 2)},
            {F(-3): F(1, 12), F(-1): F(1, 4), F(0): F(1, 3),
             F(1): F(1, 6), F(2): F(1, 6)},
            {F(-2): F(1, 6), F(-1): F(1, 3), F(1): F(1, 3), F(2): F(1, 6)}]
    decomposed = [mixture(law) for law in laws]
    for m in range(2, 6):
        original = [laws[i % len(laws)] for i in range(m)]
        mixes = [decomposed[i % len(laws)] for i in range(m)]
        mixed_e = mixed_v = nondegenerate_mass = F(0)
        components = 0
        c2 = F(m-1, 2*m)
        for row in product(*mixes):
            weight = prod(w for w, _ in row)
            kernels = [k for _, k in row]
            e = expected_max(kernels)
            v = sum(map(second, kernels), F(0))
            assert e >= 0 and e**2 <= c2*v
            if all(second(k) > 0 for k in kernels):
                assert e**2 < c2*v
                nondegenerate_mass += weight
            mixed_e += weight*e
            mixed_v += weight*v
            components += 1
        assert mixed_e == expected_max(original)
        assert mixed_v == sum(map(second, original), F(0))
        assert nondegenerate_mass == prod(1-law.get(F(0), 0) for law in original)
        assert nondegenerate_mass > 0 and mixed_e**2 < c2*mixed_v
        print(f'MIXTURE M={m}: components={components}; Emax={mixed_e}; '
              f'V={mixed_v}; all-nondegenerate mass={nondegenerate_mass}')


if __name__ == '__main__':
    symbolic_checks()
    cascade_checks()
    grid_checks()
    mixture_checks()
    print('ALL_CHECKS_OK')
