"""Test the cached complete-comparison projection ceiling after one shift.

Only U2 + tau*rV*rV/h and one native18determinant are new. The
sources, projection, actions, energies, inverses and far Gram are reused.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, ctx, fmpq

from full_gamma_correlated_sign import symmetric
from full_gamma_small_shift_exclusion import MODES
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
TAU = fmpq(1, 48)
PINS = {
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
    'full-gamma-prime-potential-resolvent.json':
        '4903fe5c85831b7290c0b713c7a689341d38425f6abd394eef69ae9aac230ca8',
}


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    guards = {source: hashlib.sha256(source.read_bytes()).hexdigest()}
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            path = Path(filename).resolve()
            guards[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted packet: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    projection, budget, eligibility = packets
    if ([p['schema'] for p in packets] != [
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_REFINED_RESOLVENT_BUDGET_V1']
            or projection['grid_bits'] != BITS or projection['remaining_modes'] != MODES
            or projection['determinant_sign'] != -1
            or eligibility['grid_bits'] != BITS or eligibility['fixed_parameter'] != ['1', '48']
            or eligibility['verdict'] != 'conditional-fixed-shift-budget-admitted'
            or budget['packet_sha256']['full-gamma-projection-upper.json.zlib'] !=
                PINS['full-gamma-projection-upper.json.zlib']
            or eligibility['packet_sha256']['full-gamma-surplus-budget.json'] !=
                PINS['full-gamma-surplus-budget.json']
            or projection['exact_source_integer_sha256'][0] != budget['source_integer_sha256']
            or eligibility['source_integer_sha256'] != budget['source_integer_sha256']
            or budget['source_integer_sha256'] !=
                '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
            or projection['exact_source_bits'][0] != 384
            or budget['source_bits'] != 384 or budget['source_row'] != 0):
        raise ValueError('Different common projection/source or admitted fixed shift.')
    U_pairs = projection['matrices']['threshold_upper']
    r_pairs = projection['matrices']['source_on_lift'][0]
    if len(U_pairs) != 18 or any(len(row) != 18 for row in U_pairs) or len(r_pairs) != 18:
        raise ValueError('Incomplete saved18compression or source row.')
    for row in U_pairs+[r_pairs]:
        for pair in row:
            if len(pair) != 2 or int(pair[0]) > int(pair[1]):
                raise ValueError('Invalid saved endpoints.')
    hlo, hhi = (fmpq(*map(int, budget['exact_rationals'][key]))
                 for key in ('h_lower', 'h_upper'))
    if not 0 < hlo <= hhi:
        raise ValueError('Require the SAME positive true-h enclosure.')
    h = arb(arb((hlo+hhi)/2), arb((hhi-hlo)/2))
    U = symmetric(matrix_in(U_pairs, BITS))
    r = matrix_in([r_pairs], BITS)
    shifted = symmetric(U+(arb(TAU)/h)*(r.transpose()*r))
    if any(not value.is_finite() for value in shifted.entries()):
        raise ArithmeticError('Nonfinite new matrix; stop without retry.')
    determinant = shifted.det()
    det_bits = 18*BITS+64
    pair = envelope(determinant, det_bits)
    if pair is None:
        raise ArithmeticError('Nonfinite new determinant; stop without retry.')
    sign = -1 if int(pair[1]) < 0 else (1 if int(pair[0]) > 0 else 0)
    verdict = {-1: 'conditional-shifted-comparison-family-excluded',
               1: 'shifted-projection-upper-positive-comparison-undecided',
               0: 'inconclusive-shifted-projection-determinant'}[sign]
    result = dict(schema='ACTUAL_C27_SHIFTED_PROJECTION_UPPER_V1',
                  packet_sha256=PINS, input_sha256=projection['input_sha256'],
                  remaining_modes=MODES, grid_bits=BITS, working_bits=1024,
                  determinant_grid_bits=det_bits, determinant_interval=pair,
                  determinant_sign=sign, shifted_threshold_upper=matrix_out(shifted, BITS),
                  fixed_parameter=['1', '48'], h_lower=str(hlo), h_upper=str(hhi),
                  source_integer_sha256=budget['source_integer_sha256'],
                  source_bits=384, source_row=0, source_sha256=guards[source],
                  invocation_source_sha256={p.name: digest for p, digest in guards.items()},
                  verdict=verdict,
                  conditions=list(dict.fromkeys(eligibility['conditions']))+[
                      'The inherited common positive inverse forms, legal E elimination and unitary source/lift identifications apply to the SAME full effective M and actual exterior short.',
                      'The accepted two-source projection gives V*(Q-R*M^-1*R)V<=U2 for the full comparison, not only a finite exterior truncation.',
                      'The SAME full u=Mw and true h=w*u give M_tau^-1=M^-1-tau*w*w*/h, rV=w*R*V; conditional eligibility at tau1/48 is inherited.',
                      'The archived U2 and rV enclosures contain one common exact symmetric matrix and row; the accepted U2 has at least17positive directions.',
                      'Every admitted inverse-majorant/residual-trial lower construction for M_tau is dominated by Q-R*M_tau^-1*R on these SAME forms.'],
                  scope='One new cached18matrix and determinant at fixed tau1/48. A negative '
                      'upper determinant excludes whole-retained positivity certification by '
                      'this fixed M_tau comparison, including the fully transported majorant. '
                      'Positive/zero-crossing upper gives no comparison feasibility or actual '
                      'positivity. No actual Weil negativity, cofinal positivity, Robin, RH, '
                      'Lean or generic novelty claim.',
                  native_determinant_calls=1, new_inverse_or_solve_calls=0,
                  old_producers_replayed=0, new_source_action_energy_Gram_far_calls=0,
                  parameter_source_precision_scan=False, precision_retry=False,
                  elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed during consumption: '+name)
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed during consumption: '+path.name)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':')).encode(), 9))
    print('new shifted determinant', determinant.str(18), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use fresh .zlib output; accepted archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
