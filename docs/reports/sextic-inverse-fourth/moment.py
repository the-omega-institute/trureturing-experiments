#!/usr/bin/env python3
"""Finite full-row Eisenstein sextic moments, with exact character labels.

Character arithmetic is exact; smooth weights and moments use binary64.
No finite output is an asymptotic bound or a Robin/RH certificate.
"""
import argparse
from collections import defaultdict
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

ROOTS = np.exp((math.pi * 1j / 3) * np.arange(6))
UNIT_DFT = ROOTS[(np.arange(6)[:, None] * np.arange(6)) % 6]


def bump(x, center, radius):
    t = (np.asarray(x, dtype=float) - center) / radius
    out = np.zeros_like(t)
    inside = np.abs(t) < 1
    out[inside] = np.exp(1 - 1 / (1 - t[inside] ** 2))
    return out


def square_plans(engine, columns, weights, D):
    """Literal unordered d1,d2 expansion, grouped by parity core and unit phase."""
    grouped = defaultdict(list)
    diagonal = {}
    head = {}
    active = [i for i, d in enumerate(columns) if engine.mu(d)]
    for k, i in enumerate(active):
        d = columns[i]
        f1 = engine.factors[d]
        for j in active[k:]:
            e = columns[j]
            f2 = engine.factors[e]
            coefficient = (1 if i == j else 2) * engine.mu(d) * engine.mu(e)
            coefficient *= weights[i] * weights[j] / D
            a = mul(d, e)
            core = tuple(sorted(set(f1) ^ set(f2)))
            h = (engine.unit_weight(d) + engine.unit_weight(e)) % 6
            grouped[(core, h)].append((i, j, coefficient))
            if a not in diagonal:
                diagonal[a] = [i, j, 0.0]
            diagonal[a][2] += coefficient
            if not set(f1).intersection(f2):
                if a not in head:
                    head[a] = [i, j, 0.0, h]
                head[a][2] += coefficient
    single = []
    multiple = []
    for (core, h), terms in grouped.items():
        if len(terms) == 1:
            single.append((*terms[0], h))
        else:
            multiple.append((core, h, terms))
    return single, multiple, list(diagonal.values()), list(head.values())


def pair_values(chars, terms):
    if not terms:
        return np.zeros((len(chars), 0), dtype=complex)
    i, j, c = zip(*(term[:3] for term in terms))
    return chars[:, i] * chars[:, j] * np.asarray(c)


def character_matrix(engine, columns, rows):
    primes = sorted({p for d in columns for p in engine.factors[d]})
    prime_data = {
        p: np.fromiter((engine.prime_symbol_exponent(p, u) for u in rows),
                       dtype=np.int16, count=len(rows))
        for p in primes
    }
    chars = np.ones((len(rows), len(columns)), dtype=complex)
    collisions = 0
    for j, d in enumerate(columns):
        exponent = np.zeros(len(rows), dtype=np.int16)
        zero = np.zeros(len(rows), dtype=bool)
        for p, v in engine.factors[d].items():
            data = prime_data[p]
            zero |= data < 0
            exponent += v * np.maximum(data, 0)
        chars[:, j] = ROOTS[exponent % 6]
        chars[zero, j] = 0
        collisions += int(np.count_nonzero(zero))
    return chars, collisions


def run_scale(engine, U, chunk):
    started = time.perf_counter()
    D = math.sqrt(U)
    columns = [d for d in engine.ideals if D / 2 < norm(d) < 1.5 * D]
    rows = [u for u in engine.ideals if U < norm(u) < 2 * U
            and engine.is_sixth_free(u)]
    if not columns or not rows:
        raise ValueError('empty arithmetic annulus')
    v = bump([norm(d) / D for d in columns], 1, .5)
    mu = np.array([engine.mu(d) for d in columns])
    lam = np.array([engine.liouville(d) for d in columns])
    unit = np.array([engine.unit_weight(d) for d in columns])
    unit_matrix = UNIT_DFT[unit]
    coefficients = {
        'mobius': mu * v / math.sqrt(D),
        'squarefree_positive': np.abs(mu) * v / math.sqrt(D),
        'liouville': lam * v / math.sqrt(D),
        'ordinary': v / math.sqrt(D),
    }
    plans = square_plans(engine, columns, v, D)
    single, multiple, diagonal, head = plans
    stats = defaultdict(float)
    error = 0.0
    collisions = 0
    for start in range(0, len(rows), chunk):
        part = rows[start:start + chunk]
        w = bump([norm(u) / U for u in part], 1.5, .5)
        chars, nzero = character_matrix(engine, columns, part)
        collisions += nzero
        values = {}
        for name, coeff in coefficients.items():
            values[name] = (chars * coeff) @ unit_matrix
            z = np.abs(values[name]) ** 2
            stats[name + '_second'] += float(np.dot(w, z.sum(axis=1)))
            stats[name + '_fourth'] += float(np.dot(w, (z * z).sum(axis=1)))
        # Six units are integrated by their exact finite Fourier orthogonality.
        harmonics = np.zeros((len(part), 6), dtype=complex)
        h2 = np.zeros(len(part))
        terms = pair_values(chars, single)
        if single:
            h2 += 6 * np.sum(np.abs(terms) ** 2, axis=1)
            hs = np.array([x[3] for x in single])
            for h in range(6):
                harmonics[:, h] += terms[:, hs == h].sum(axis=1)
        for core, h, pairs in multiple:
            amplitude = pair_values(chars, pairs).sum(axis=1)
            h2 += 6 * np.abs(amplitude) ** 2
            harmonics[:, h] += amplitude
        squares = harmonics @ UNIT_DFT
        expected = values['mobius'] ** 2
        error = max(error, float(np.max(np.abs(squares - expected)) /
                                 max(1.0, float(np.max(np.abs(expected))))))
        if error > 2e-11:
            raise ArithmeticError('literal inverse-square expansion failed')
        diag = 6 * np.sum(np.abs(pair_values(chars, diagonal)) ** 2, axis=1)
        head_terms = pair_values(chars, head)
        head_harmonics = np.zeros_like(harmonics)
        if head:
            head_hs = np.array([x[3] for x in head])
            for h in range(6):
                head_harmonics[:, h] = head_terms[:, head_hs == h].sum(axis=1)
        head_values = head_harmonics @ UNIT_DFT
        head_energy = (np.abs(head_values) ** 2).sum(axis=1)
        head_diag = 6 * np.sum(np.abs(head_terms) ** 2, axis=1)
        tail_energy = (np.abs(squares - head_values) ** 2).sum(axis=1)
        stats['h2'] += float(np.dot(w, h2))
        stats['diagonal'] += float(np.dot(w, diag))
        stats['t1_energy'] += float(np.dot(w, head_energy))
        stats['t1_diagonal'] += float(np.dot(w, head_diag))
        stats['t_gt1_energy'] += float(np.dot(w, tail_energy))
        stats['row_mass'] += 6 * float(w.sum())
    stats['e6'] = stats['mobius_fourth'] - stats['h2']
    stats['h2_off_diagonal'] = stats['h2'] - stats['diagonal']
    stats['t1_e6'] = stats['t1_energy'] - stats['t1_diagonal']
    stats['t1_tail_interference'] = (stats['mobius_fourth'] -
                                   stats['t1_energy'] - stats['t_gt1_energy'])
    return {
        'U': U, 'D': D, 'primary_rows': len(rows),
        'physical_rows': 6 * len(rows), 'columns': len(columns),
        'squarefree_columns': int(np.count_nonzero(mu)),
        'inverse_square_columns': len(diagonal), 't1_columns': len(head),
        'collision_zeros_primary_matrix': collisions,
        'square_expansion_relative_error': error,
        'coefficient_norm_square': float(sum(x[2] ** 2 for x in diagonal)),
        'seconds': time.perf_counter() - started,
        **dict(stats),
        'per_U': {k: value / U for k, value in stats.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scales', type=int, nargs='+', required=True)
    parser.add_argument('--chunk', type=int, default=256)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if min(args.scales) < 64 or args.chunk < 1:
        parser.error('scales must be >=64 and chunk positive')
    started = time.perf_counter()
    engine = Engine(2 * max(args.scales))
    results = []
    for scale in args.scales:
        result = run_scale(engine, scale, args.chunk)
        results.append(result)
        print(json.dumps({k: result[k] for k in
                          ['U', 'physical_rows', 'columns', 'mobius_fourth',
                           'h2', 'e6', 't1_e6', 'seconds']}), flush=True)
    source = Path(__file__).resolve()
    report = {
        'scope': 'finite complete good sixth-power-free element rows; fixed trivial ray twist',
        'bad_rational_primes': [2, 3], 'units_per_ideal': 6,
        'character': 'chi_d(u), exact sextic residue labels in 0..5; collision zero',
        'profiles': {'V': 'bump centered 1 radius 1/2',
                     'Omega': 'bump centered 3/2 radius 1/2',
                     'bump': 'exp(1-1/(1-t^2)) for |t|<1, zero otherwise'},
        'normalization': 'D^-1/2 sum_d coefficient(d) V(Nd/D) chi_d(u), D=sqrt(U)',
        'h2': 'sum over parity cores of the squared inverse-square block; all6 units',
        'e6': 'mobius_fourth-h2, signed; not assumed nonnegative',
        't1': 'actual shared-gcd-one coefficient, with interference retained in full moment',
        'arithmetic': 'integer factorization and exact finite-field character exponents',
        'numerics': 'binary64 weights/moments, no certified inequality or asymptotic exponent',
        'python': platform.python_version(), 'numpy': np.__version__,
        'sources': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [source, source.with_name('eisenstein.py')]},
        'total_seconds': time.perf_counter() - started,
        'results': results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
