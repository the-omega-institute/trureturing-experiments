"""Consume archived w0 data for the necessary surplus-shift quadratic.

One small native solve supplies the missing exact-Xi contractions. Exact
rational coefficient bounds then enclose the minimum on [0, 1], without
parameter search or any replay of source actions or scientific suppliers.
"""
import argparse
from fractions import Fraction
import hashlib
from pathlib import Path
import platform
import sys
import time

import flint
from flint import arb_mat, ctx

from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


HERE = Path(__file__).resolve().parent
PINS = {
    'sources': '34fca6954cdd706a11479cb68fa30995dd7586e0207c03ee6c3b68742d0896ed',
    'projection': 'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'gram': '0d28fdab5adb68146ae44735809cbfeb2cd0de55bac91143bd4548d5d69501c7',
}


def minimum(coefficients):
    """Exact endpoint/convex-vertex minimum of a rational quadratic."""
    a, b, c = coefficients
    candidates = [(Fraction(0), a), (Fraction(1), a+b+c)]
    if c > 0:
        vertex = -b/(2*c)
        if 0 < vertex < 1:
            candidates.append((vertex, a+b*vertex+c*vertex*vertex))
    return min(candidates, key=lambda pair: pair[1])


def rational(value):
    return [str(value.numerator), str(value.denominator)]


def consume(args):
    started = time.monotonic()
    guards = {}
    packets = {}
    for role, digest in PINS.items():
        path = getattr(args, role)
        packet, actual = read_json(path)
        if actual != digest:
            raise ValueError('Require the accepted '+role+' packet.')
        packets[role], guards[path] = packet, actual
    # Invocation guards only; no cache or reuse decision depends on code bytes.
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            path = Path(filename).resolve()
            guards[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    source = Path(__file__).resolve()
    guards[source] = hashlib.sha256(source.read_bytes()).hexdigest()
    sources, projection, gram = (packets[k] for k in ('sources', 'projection', 'gram'))
    if (sources['schema'], projection['schema'], gram['schema']) != (
            'ACTUAL_C27_TWO_SOURCE_JOINT_LOWER_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1', 'ACTUAL_C27_LIFT_FINITE_GRAM_V1'):
        raise ValueError('Different accepted objects.')
    modes = [0, 1, 2, 3, 4, 5, 6, 9, 10, 12, 14, 15, 18, 19, 24, 29, 33, 39]
    if any(p['grid_bits'] != 640 or p['remaining_modes'] != modes
           for p in packets.values()):
        raise ValueError('Different exact lift or entry grid.')
    if (projection['finite_gram_sha256'] != PINS['gram'] or
            projection['packet_sha256']['sources'] != PINS['sources'] or
            projection['exact_source_bits'][0] != 384 or
            projection['exact_source_integer_sha256'][0] !=
            '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'):
        raise ValueError('Different original w0 or packet correspondence.')
    if any((p['N'], p['M'], p['finite_first'], p['finite_last']) !=
           (1280, 32768, 2049, 32768) for p in (projection, gram)):
        raise ValueError('Different finite support.')
    if any(p['input_sha256'] != gram['input_sha256'] for p in packets.values()):
        raise ValueError('Different actual scientific inputs.')
    for role, digest in gram['packet_sha256'].items():
        if projection['packet_sha256'].get(role) != digest:
            raise ValueError('Different shared packet role: '+role)
    for role in ('base', 'boundary', 'actions', 'core', 'projection', 'gap'):
        if sources['packet_sha256'][role] != gram['packet_sha256'][role]:
            raise ValueError('Different actual base/lift/far packet: '+role)

    kJ = interval_ball(sources['finite_charge_matrix'][0][0], 640)
    h = interval_ball(projection['matrices']['source_gram'][0][0], 640)
    b = matrix_in([sources['finite_row_on_lift'][0]], 640)
    r = matrix_in([projection['matrices']['source_on_lift'][0]], 640)
    Xi = matrix_in(gram['matrices']['far_shifted_core'], 640)
    if (b.nrows(), b.ncols(), r.nrows(), r.ncols(), Xi.nrows(), Xi.ncols()) != (
            1, 18, 1, 18, 18, 18) or not h > 0:
        raise ValueError('Require the same eighteen-dimensional rows and positive h.')
    if any(Xi[i, j] != Xi[j, i] for i in range(18) for j in range(18)):
        raise ValueError('Require the symmetric archived Xi enclosure.')
    rhs = arb_mat([[b[0, i], r[0, i]] for i in range(18)])
    solution = Xi.solve(rhs, algorithm='precond')
    contraction = rhs.transpose()*solution
    if any(not x.is_finite() for x in solution.entries()+contraction.entries()):
        raise ArithmeticError('Nonfinite solve; stop without a precision retry.')
    if not contraction[0, 1].overlaps(contraction[1, 0]):
        raise ArithmeticError('Inconsistent common Hermitian cross contraction.')
    coefficients = [kJ+contraction[0, 0],
                    h+contraction[0, 1]+contraction[1, 0], contraction[1, 1]]
    pairs = [envelope(value, 640) for value in coefficients]
    if any(pair is None for pair in pairs):
        raise ArithmeticError('Nonfinite coefficient enclosure.')
    # t and t^2 are nonnegative: coefficient endpoints bound scalar
    # polynomials on this interval, without any matrix endpoint ordering.
    low = [Fraction(int(pair[0]), 1 << 640) for pair in pairs]
    high = [Fraction(int(pair[1]), 1 << 640) for pair in pairs]
    low_t, low_min = minimum(low)
    high_t, high_min = minimum(high)
    if low_min > high_min:
        raise ArithmeticError('Reversed exact minimum enclosure.')
    verdict = ('fixed-source-shift-strict-sufficiency-excluded' if low_min >= 0 else
               'negative-necessary-region-only' if high_min < 0 else
               'inconclusive-necessary-quadratic-minimum')
    result = dict(
        schema='ACTUAL_C27_SURPLUS_SHIFT_NECESSARY_GATE_V1',
        packet_sha256=PINS, source_sha256=guards[source],
        consumed_sha256={path.name: digest for path, digest in guards.items()},
        remaining_modes=modes, N=1280, M=32768, finite_first=2049, finite_last=32768,
        source_row=0, source_bits=384, grid_bits=640, working_bits=args.precision,
        coefficient_order=['constant', 'linear', 'quadratic'],
        coefficient_intervals=pairs, inverse_solution=matrix_out(solution, 640),
        common_contractions=matrix_out(contraction, 640),
        minimum_lower_rational=rational(low_min), minimum_upper_rational=rational(high_min),
        lower_polynomial_minimizer=rational(low_t),
        upper_polynomial_minimizer=rational(high_t),
        minimum_domain=['0', '1'], verdict=verdict,
        native_solve_shape=[18, 18, 2], native_solve_calls=1, parameter_scan=False,
        old_scientific_productions_replayed=0, new_source_actions=0,
        runtime=dict(python=platform.python_version(), flint=flint.__version__,
                     machine=platform.machine()),
        elapsed_seconds=time.monotonic()-started,
        conditions=[
            'Same exact archived w0, complete-T effective action, retained residual and lift.',
            'Common actual forms/domains/inverses, same base and equal-inertia complement order.',
            'Same complete paid far Gram, unscaled actual u far energy, and sound enclosures.',
            'An actual shifted construction additionally needs its SAME-u directional surplus budget.',
        ],
        scope='ell(t)=kJ+t*h+(bJ+t*rV)Xi^-1(bJ+t*rV)* is a necessary LOWER gate, '
              'not a sufficient upper for the actual charge/good-block/Schur score. Nonnegative '
              'minimum excludes strict success only for this fixed-w0 family. Strict positivity '
              'also excludes its positive-denominator PSD upper-compression update. A negative '
              'minimum supplies no actual budget, update or positivity certificate. The construction '
              'uses 0<t<1; closed-interval endpoints only bound its polynomial. No generic novelty, '
              'Lean, actual retained/cofinal positivity, Robin-tail or RH result.')
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('A consumed input/helper changed: '+path.name)
    write_packet(args.out, result)
    print('coefficients', [v.str(18) for v in coefficients], flush=True)
    print('minimum lower', float(low_min), 'upper', float(high_min), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for role, filename in dict(sources='full-gamma-two-source.json.zlib',
                               projection='full-gamma-projection-upper.json.zlib',
                               gram='full-gamma-lift-finite-gram.json.zlib').items():
        parser.add_argument('--'+role, type=Path, default=HERE/filename)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.precision < 832:
        parser.error('Use a fresh .zlib output and at least192 guardbits above grid640.')
    if platform.python_version() != '3.13.12' or flint.__version__ != '0.9.0':
        parser.error('Require Python3.13.12 and python-flint0.9.0.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
