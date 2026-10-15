"""Parameter application of the existing weighted Schur and Gamma suppliers.

This checks the c=27 scalar inputs, not a low-block sign or RH. The fixed
project arithmetic source is supplied read-only; it is not redistributed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import isqrt
from pathlib import Path
import platform

from mpmath import iv


POWERS = {2: (2, 1), 3: (3, 1), 4: (2, 2), 5: (5, 1),
          7: (7, 1), 8: (2, 3), 9: (3, 2), 11: (11, 1),
          13: (13, 1), 16: (2, 4), 17: (17, 1), 19: (19, 1),
          23: (23, 1), 25: (5, 2)}
BASE_SHA = '8bb067fc5499b0f2e1e48836e7a82237a15504109f82a856c72478d1096d69d0'
B, ROW = F(29, 50), F(59, 10)


def rat(q):
    q = F(q)
    return iv.mpf(q.numerator) / q.denominator


def endpoint(q):
    sign, mantissa, exponent, _ = q
    return (-1 if sign else 1) * F(mantissa) * F(2) ** exponent


def rational_minimum(aa, bb, cc, left, right):
    points = [left, right]
    if aa > 0 and left < -bb / (2 * aa) < right:
        points.append(-bb / (2 * aa))
    return min(aa * x * x + bb * x + cc for x in points)


def run(arithmetic_source, digits):
    raw = arithmetic_source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != BASE_SHA:
        raise ValueError('Wrong existing arithmetic source.')
    spec = importlib.util.spec_from_file_location('prime3_arithmetic', arithmetic_source)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)  # Defines routines; does not execute its old run().
    iv.dps = digits
    grid = 10 ** 12
    logs = {}
    for p in sorted({p for p, _ in POWERS.values()}):
        value = iv.ln(p)
        lo, hi = map(endpoint, value._mpi_)
        lo = F((lo * grid).__floor__(), grid)
        hi = F((hi * grid).__ceil__(), grid)
        if not (bool(rat(lo) < value) and bool(value < rat(hi))):
            raise ArithmeticError('Logarithm enclosure failed.')
        logs[p] = lo, hi
    weights, shifts = {}, {}
    for n, (p, j) in POWERS.items():
        if n != p ** j or not (2 <= n < 27):
            raise ValueError('Prime-power list mismatch.')
        root = isqrt(n * grid * grid)
        if not root * root <= n * grid * grid < (root + 1) ** 2:
            raise ArithmeticError('Root enclosure failed.')
        weights[n] = logs[p][1] / F(root, grid)
        shifts[n] = ((F(2 * j, 3),) * 2 if p == 3 else
                     (2 * j * logs[p][0] / (3 * logs[3][1]),
                      2 * j * logs[p][1] / (3 * logs[3][0])))
    # Every actual support switch lies in its rational outward enclosure.
    ends = {F(0), F(1)}
    for lo, hi in shifts.values():
        ends.update(x for x in [lo - 1, hi - 1, 1 - hi, 1 - lo] if 0 < x < 1)
    ends = sorted(ends)
    rows = []
    for left, right in zip(ends, ends[1:]):
        negative = [n for n, (lo, hi) in shifts.items() if right > lo - 1]
        positive = [n for n, (lo, hi) in shifts.items() if left < 1 - lo]
        ww = sum((weights[n] for n in negative + positive), F(0))
        bb = (sum((weights[n] * shifts[n][0] for n in negative), F(0))
              - sum((weights[n] * shifts[n][1] for n in positive), F(0)))
        dd = sum((weights[n] * shifts[n][1] ** 2 for n in negative + positive), F(0))
        margin = rational_minimum(B * (ROW - ww), 2 * B * bb,
                                  ROW - ww - B * dd, left, right)
        if margin <= 0:
            raise ArithmeticError('Outward weighted row bound failed.')
        rows.append({'left': str(left), 'right': str(right), 'negative': negative,
                     'positive': positive, 'polynomial_minimum': str(margin)})
    margin = min(F(row['polynomial_minimum']) for row in rows)
    if margin <= F(1, 400):
        raise ArithmeticError('Uniform rational row margin failed.')
    length = iv.ln(27)

    def floor(n, debt):
        z = iv.mpc(rat(F(1, 4)), iv.pi * (n + 1) / length)
        return base.psi_iv(z).real - iv.ln(iv.pi) - debt

    # Classical path adjacency norm for each paired clipped translation.
    separate = iv.mpf(0)
    chains = {}
    for n, (p, j) in POWERS.items():
        k = 1
        while n ** k < 27:
            k += 1
        if not n ** (k - 1) < 27 <= n ** k:
            raise ArithmeticError('Chain length failed.')
        chains[n] = k
        separate += iv.ln(p) / iv.sqrt(n) * 2 * iv.cos(iv.pi / (k + 1))
    checks = {
        'separate_16972_negative': bool(floor(16972, separate) < 0),
        'separate_16973_positive': bool(floor(16973, separate) > 0),
        'joint_1202_negative': bool(floor(1202, rat(ROW)) < 0),
        'joint_1203_positive': bool(floor(1203, rat(ROW)) > 0),
        'joint_1280_above_one_sixteenth': bool(floor(1280, rat(ROW)) > rat(F(1, 16))),
    }
    if not all(checks.values()):
        raise ArithmeticError('Gamma floor comparison failed.')
    return {
        'cutoff': 27, 'half_width': '3log(3)/2', 'prime_powers': POWERS,
        'weight': '1+(29/50)(u/a)^2', 'row_upper': str(ROW),
        'rational_grid_denominator': grid, 'interval_digits': digits,
        'logarithm_intervals': {p: list(map(str, pair)) for p, pair in logs.items()},
        'support_cells': rows, 'minimum_polynomial_margin': str(margin),
        'all_polynomial_minima_above': '1/400',
        'separate_chain_lengths': chains, 'separate_debt_interval': str(separate),
        'separate_floor_16972': str(floor(16972, separate)),
        'separate_floor_16973': str(floor(16973, separate)),
        'joint_floor_1202': str(floor(1202, rat(ROW))),
        'joint_floor_1203': str(floor(1203, rat(ROW))),
        'joint_floor_1280': str(floor(1280, rat(ROW))), 'checks': checks,
        'arithmetic_sha256': BASE_SHA, 'python': platform.python_version(),
        'mpmath': base.mpmath.__version__, 'numpy': base.np.__version__,
        'sympy': base.sympy.__version__,
        'scope': 'Scalar input certificate conditional on existing paper-level '
                 'weighted-Schur and Neumann-Gamma/form-domain identifications. '
                 'No low-block sign, full coupling certificate, Lean proof or RH claim.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arithmetic-source', type=Path, required=True)
    parser.add_argument('--digits', type=int, default=70)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.digits < 55:
        parser.error('At least 55 interval digits are required.')
    result = run(args.arithmetic_source, args.digits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ['support_cells', 'logarithm_intervals']}, indent=2))
