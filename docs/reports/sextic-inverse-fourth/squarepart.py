#!/usr/bin/env python3
"""Original small-product-squarepart sextic correlations on finite annuli.

Flat weights give exact rational moments; the original smooth bump uses
binary64. No finite value establishes an eventual sign or an RH estimate.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
sys.dont_write_bytecode = True
import numpy as np
from eisenstein import Engine, mul, norm
from moment import bump, square_plans

# zeta_6**j=A[j]+B[j]*zeta_6 and zeta_6**2=zeta_6-1.
A = np.array([1, 0, -1, -1, 0, 1], dtype=np.int64)
B = np.array([0, 1, 1, 0, -1, -1], dtype=np.int64)
ROOT = np.exp(math.pi*1j/3)


def norm_square(a, b):
    return a*a+a*b+b*b


def product_squarepart(fa, fb):
    return math.prod(norm(p)**((fa.get(p, 0)+fb.get(p, 0))//2)
                     for p in fa.keys() | fb.keys())


def inverse_plan(engine, columns, weights, D):
    """Literal integer flat counts; reuse the original smooth coefficients."""
    flat = {}
    for i, d in enumerate(columns):
        for j in range(i, len(columns)):
            a = mul(d, columns[j])
            if a not in flat:
                flat[a] = [i, j, 0]
            flat[a][2] += (1 if i == j else 2)*engine.mu(d)*engine.mu(columns[j])
    smooth = {mul(columns[i], columns[j]): c
              for i, j, c in square_plans(engine, columns, weights, D)[2]}
    if set(smooth) != set(flat):
        raise ArithmeticError('different inverse-square supports')
    result = []
    for a, (i, j, count) in flat.items():
        fa = engine._factor(a)
        if any(v > 2 for v in fa.values()):
            raise ArithmeticError('non-cube-free inverse coefficient')
        core = tuple(p for p, v in sorted(fa.items()) if v % 2)
        result.append((a, i, j, count, smooth[a], core, engine.unit_weight(a)))
    return result


def pair_plan(engine, plan, U):
    pairs = []
    for i, a in enumerate(plan):
        for j in range(i+1, len(plan)):
            b = plan[j]
            if a[5] == b[5] or a[6] != b[6]:
                continue
            fa, fb = engine._factor(a[0]), engine._factor(b[0])
            P = product_squarepart(fa, fb)
            parity_norm = math.prod(norm(p) for p in fa.keys() | fb.keys()
                                    if (fa.get(p, 0)+fb.get(p, 0)) % 2)
            if P*P*parity_norm != norm(a[0])*norm(b[0]):
                raise ArithmeticError('squarepart identity failed')
            pairs.append((i, j, P, P**4 < U))
    return pairs


def run_scale(engine, U):
    started = time.perf_counter()
    D = math.sqrt(U)
    columns = [d for d in engine.ideals
               if U < 4*norm(d)**2 < 9*U and engine.mu(d)]
    rows = [u for u in engine.ideals
            if U < norm(u) < 2*U and engine.is_sixth_free(u)]
    if not columns or not rows:
        raise ValueError('empty arithmetic annulus')
    if 6*len(rows)*len(columns)**4 >= 2**62:
        raise ValueError('requested integer array accumulation exceeds safe bound')
    weights = bump([norm(d)/D for d in columns], 1, .5)
    plan = inverse_plan(engine, columns, weights, D)
    pairs = pair_plan(engine, plan, U)
    labels = np.array([[engine.symbol_exponent(d, u) for d in columns]
                       for u in rows], dtype=np.int16)
    i1, i2 = np.array([a[1] for a in plan]), np.array([a[2] for a in plan])
    zero = (labels[:, i1] < 0) | (labels[:, i2] < 0)
    expanded = (labels[:, i1]+labels[:, i2]) % 6
    ex_a, ex_b = A[expanded].copy(), B[expanded].copy()
    ex_a[zero] = 0
    ex_b[zero] = 0
    counts = np.array([a[3] for a in plan], dtype=np.int64)
    coeff = np.array([a[4] for a in plan])
    omega = bump([norm(u)/U for u in rows], 1.5, .5)
    groups = defaultdict(list)
    for j, item in enumerate(plan):
        groups[(item[5], item[6])].append(j)
    flat_h2 = 0
    smooth_h2_rows = np.zeros(len(rows))
    chars = ex_a+ex_b*ROOT
    for indices in groups.values():
        ia, ib = ex_a[:, indices]@counts[indices], ex_b[:, indices]@counts[indices]
        flat_h2 += 6*int(norm_square(ia, ib).sum())
        smooth_h2_rows += 6*np.abs(chars[:, indices]@coeff[indices])**2
    flat_fourth = 0
    smooth_fourth_rows = np.zeros(len(rows))
    mu = np.array([engine.mu(d) for d in columns], dtype=np.int64)
    unit = np.array([engine.unit_weight(d) for d in columns])
    for k in range(6):
        phases = (labels+k*unit) % 6
        aa, bb = A[phases].copy(), B[phases].copy()
        aa[labels < 0] = 0
        bb[labels < 0] = 0
        flat_fourth += int((norm_square(aa@mu, bb@mu)**2).sum())
        m = (aa+bb*ROOT)@(mu*weights/math.sqrt(D))
        smooth_fourth_rows += np.abs(m)**4
    flat_small = flat_large = 0
    smooth_small = smooth_large = 0.0
    for i, j, P, small in pairs:
        trace = (2*ex_a[:, i]*ex_a[:, j]+ex_a[:, i]*ex_b[:, j]
                 +ex_b[:, i]*ex_a[:, j]+2*ex_b[:, i]*ex_b[:, j])
        value = 6*int(counts[i])*int(counts[j])*int(trace.sum())
        numerical = 6*coeff[i]*coeff[j]*float(np.dot(omega, trace))
        if small:
            flat_small += value
            smooth_small += numerical
        else:
            flat_large += value
            smooth_large += numerical
    flat_e6 = flat_fourth-flat_h2
    if flat_e6 != flat_small+flat_large:
        raise ArithmeticError('exact signed decomposition failed')
    smooth_fourth = float(np.dot(omega, smooth_fourth_rows))
    smooth_h2 = float(np.dot(omega, smooth_h2_rows))
    smooth_e6 = smooth_fourth-smooth_h2
    error = abs(smooth_e6-smooth_small-smooth_large)
    if error > 1e-10*max(1.0, abs(smooth_fourth), abs(smooth_h2)):
        raise ArithmeticError('smooth signed decomposition failed')
    def exact(n):
        value = Fraction(n, U)
        return {'numerator': value.numerator, 'denominator': value.denominator,
                'decimal': float(value)}
    return {'U': U, 'D': D, 'physical_rows': 6*len(rows),
            'squarefree_inverse_columns': len(columns), 'inverse_square_columns': len(plan),
            'small_P_unordered_unit_compatible_pairs': sum(p[3] for p in pairs),
            'all_E6_unordered_unit_compatible_pairs': len(pairs),
            'flat': {name: exact(n) for name, n in
                     [('fourth', flat_fourth), ('h2', flat_h2), ('e6', flat_e6),
                      ('small_P', flat_small), ('large_P', flat_large)]},
            'bump': {'fourth': smooth_fourth, 'h2': smooth_h2, 'e6': smooth_e6,
                     'small_P': smooth_small, 'large_P': smooth_large,
                     'decomposition_absolute_error': error},
            'seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scales', type=int, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if min(args.scales) < 64:
        parser.error('scales must be >=64')
    engine = Engine(2*max(args.scales))
    results = []
    for U in args.scales:
        row = run_scale(engine, U)
        results.append(row)
        print(json.dumps({'U': U, 'flat_small_P': row['flat']['small_P'],
                          'bump_small_P': row['bump']['small_P'],
                          'seconds': row['seconds']}), flush=True)
    source = Path(__file__).resolve()
    report = {'scope': 'complete original good sixth-power-free rows; six units; trivial fixed ray twist',
              'mask': 'unequal parity cores; exact P^4<U; P=root norm of squarepart(ab)',
              'arithmetic': 'existing exact Eisenstein factorization and sextic labels; integer flat traces, normalization denominator U',
              'flat_profile': 'unit weights on D/2<Nd<3D/2 and U<Nu<2U; realizable by smooth annular functions on each finite grid, not one asserted fixed smooth profile across scales',
              'bump_profile': 'same fixed smooth bump as moment.py; binary64, no interval certificate',
              'limitations': 'finite sign checks only; no eventual sign, asymptotic exponent, inverse-fourth gain, Robin or RH result',
              'python': platform.python_version(), 'numpy': np.__version__,
              'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [source, source.with_name('moment.py'), source.with_name('eisenstein.py')]},
              'results': results}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
