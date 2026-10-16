"""Transport the cached exact trial to the supplied band-budget endpoint.

Reuse the old first coordinate, mixed source scalar, curvature, signed
trial lower and full far energy. Only one missing18term mixed Gram
contraction is new. No matrix solve, inverse or old producer is run.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, ctx, fmpq

from diagonal_floor_witness import interval_ball
from full_gamma_small_shift_exclusion import MODES, PINS as OLD_PINS
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = dict(OLD_PINS, **{
    'full-gamma-small-shift-exclusion.json':
        '73d72f7dd45b44bcc3a919b9112e37450070eb82e937fdababcbfe90e4f2bedc',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
})


def ball(value):
    return arb(fmpq(value.numerator, value.denominator))


def parameters(a, curvature, r0, h, tau):
    """No complement or trial is rebuilt by this scalar rank-one update."""
    mu = ball(tau)/h
    eta = 1+mu*curvature
    pivot = a+mu*r0**2
    values = dict(shift_over_true_h=mu, pivot_ratio=eta,
                  transported_negative_pivot=pivot)
    if not (a < 0 and h > 0 and 0 <= tau < 1 and eta > 0 and pivot < 0):
        return values, 'inconclusive-transported-negative-pivot'
    return values, None


def complete(values, mixed, g00, old_energy, cross, qlower, charge):
    """Complete the SAME far square after the exact trial transport."""
    mu, eta = values['shift_over_true_h'], values['pivot_ratio']
    pivot = values['transported_negative_pivot']
    kappa = arb(fmpq(33, 32))
    coefficient = kappa-kappa**2*g00/(-pivot)
    values['transported_far_quadratic_coefficient'] = coefficient
    if not coefficient > 0:
        return values, 'inconclusive-transported-far-coefficient'
    theta = values['transport_coefficient']
    energy = old_energy-2*theta*cross+theta**2*g00
    if energy.upper() < 0:
        raise ArithmeticError('A PSD Gram energy has a negative upper; invalid inputs.')
    # Intersect a zero-crossing enclosure with the known nonnegative energy.
    energy_upper = max(arb(0), energy.upper())
    correction = mu*mixed**2/eta
    penalty = kappa**2*energy_upper/coefficient
    lower = charge+qlower-correction.upper()-penalty.upper()
    values.update(transported_far_gradient_energy=energy,
                  finite_signed_trial_decrement=correction,
                  transported_signed_trial_lower=qlower-correction.upper(),
                  transported_completed_square_penalty_upper=penalty.upper(),
                  transported_family_exclusion_lower=lower.lower())
    verdict = ('conditional-transported-supplied-interval-excluded' if lower > 0
               else 'inconclusive-transported-fixed-trial-lower')
    return values, verdict


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and Path(path).resolve().parent == HERE:
            invoked.add(Path(path).resolve())
    guards = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    gap, joint, lift, projection, actual, core, budget, small, band = packets.values()
    schemas = [
        'ACTUAL_C27_FULL_GAMMA_BASELINE_GAP_V1',
        'ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1',
        'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
        'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
        'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
        'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
        'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1',
        'ACTUAL_C27_SMALL_SHIFT_EXCLUSION_V1',
        'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1']
    if [p['schema'] for p in packets.values()] != schemas:
        raise ValueError('Different accepted packet roles.')
    if (small['packet_sha256'] != OLD_PINS
            or small['exact_trial_source'] != 'joint.new_reference_profile'
            or gap['exact_trial_source'] != small['exact_trial_source']
            or any(p['grid_bits'] != BITS or p['remaining_modes'] != MODES
                   or p['input_sha256'] != gap['input_sha256']
                   for p in (gap, joint, lift, projection, actual, core, small))):
        raise ValueError('Different common source, lift, grid or exact trial.')
    for key, filename in dict(joint='full-gamma-actual-joint-lower.json.zlib',
                              projection='full-gamma-lift-projection.json.zlib',
                              actual='full-gamma-actual-channel-projection.json.zlib',
                              core='joint-gamma-6765-rho32-v2-schur.json').items():
        if gap['packet_sha256'][key] != PINS[filename]:
            raise ValueError('Different baseline dependency: '+key)
    if (band['packet_sha256']['full-gamma-surplus-budget.json'] != PINS['full-gamma-surplus-budget.json']
            or band['packet_sha256']['full-gamma-projection-upper.json.zlib'] != PINS['full-gamma-projection-upper.json.zlib']
            or projection['packet_sha256']['gap'] != PINS['full-gamma-baseline-gap.json.zlib']
            or projection['packet_sha256']['projection'] != PINS['full-gamma-lift-projection.json.zlib']
            or projection['packet_sha256']['core'] != PINS['joint-gamma-6765-rho32-v2-schur.json']
            or budget['packet_sha256']['full-gamma-projection-upper.json.zlib'] != PINS['full-gamma-projection-upper.json.zlib']
            or band['verdict'] != 'conditional-band-budget-below-projection-capacity-target'
            or band['grid_bits'] != BITS):
        raise ValueError('Different accepted band-budget or same-source correspondence.')
    for p in (band, budget):
        if (p['source_bits'] != 384 or p['source_row'] != 0
                or p['source_integer_sha256'] != projection['exact_source_integer_sha256'][0]
                or p['source_bits'] != projection['exact_source_bits'][0]):
            raise ValueError('Different original unscaled384-bit row0 source.')
    if small['source_sha256'] != guards[HERE/'full_gamma_small_shift_exclusion.py']:
        raise ValueError('Different saved scalar producer.')
    G = lift['far_gram_upper']
    r = projection['matrices']['source_on_lift'][0]
    x = joint['new_reference_profile']
    if (len(G) != 18 or any(len(row) != 18 for row in G) or len(r) != 18
            or len(x) != 17 or any(len(row) != 1 for row in x)):
        raise ValueError('Incomplete raw paid Gram, source row or exact17trial.')
    rq = budget['exact_rationals']
    hlo, hhi = (F(*map(int, rq[key])) for key in ('h_lower', 'h_upper'))
    tau = F(*map(int, band['exact_rationals']['supplied_parameter_cap']))
    if not (0 < hlo <= hhi and 0 < tau < 1):
        raise ValueError('Invalid true h interval or supplied endpoint.')
    h = ball((hlo+hhi)/2)+arb(0, ball((hhi-hlo)/2))
    saved = small['scalar_bounds']
    old = gap['scalar_bounds']
    read = lambda pair: interval_ball(pair, BITS)
    r0 = read(r[0])
    values, verdict = parameters(read(old['shifted_negative_pivot']),
        read(saved['new_signed_curvature']), r0, h, tau)
    contractions = 0
    if verdict is None:
        values['transport_coefficient'] = (
            values['shift_over_true_h']*r0*read(saved['new_mixed_scalar'])
            /values['transported_negative_pivot'])
        g0 = [read(saved['new_gradient_first_coordinate'])]+[arb(fmpq(row[0])) for row in x]
        cross = sum((read(gij)*gj for gij, gj in zip(G[0], g0)), arb(0))
        contractions = 1
        values['new_mixed_gram_scalar'] = cross
        charge = ball(F(int(joint['scalar_bounds']['actual_common_charge_floor'][0]), 1 << BITS))
        qlower = ball(F(int(old['shifted_finite_signed_trial_lower'][0]), 1 << BITS))
        values, verdict = complete(values, read(saved['new_mixed_scalar']), read(G[0][0]),
            read(old['shifted_far_trial_gradient_energy']), cross, qlower, charge)
    exported = {key: envelope(value, BITS) for key, value in values.items()}
    if any(pair is None or len(pair) != 2 or int(pair[0]) > int(pair[1])
           for pair in exported.values()):
        raise ArithmeticError('Nonfinite or unordered scalar export; stop without retry.')
    result = dict(schema='ACTUAL_C27_SCALAR_TRANSPORT_V1', packet_sha256=PINS,
        source_sha256=guards[source],
        invocation_source_sha256={p.name: digest for p, digest in guards.items()},
        input_sha256=gap['input_sha256'], remaining_modes=MODES,
        source_integer_sha256=budget['source_integer_sha256'], source_bits=384, source_row=0,
        working_bits=1024, grid_bits=BITS, scalar_bounds=exported,
        exact_rationals=dict(h_lower=[str(hlo.numerator), str(hlo.denominator)],
            h_upper=[str(hhi.numerator), str(hhi.denominator)],
            supplied_parameter_cap=[str(tau.numerator), str(tau.denominator)]),
        exact_trial_source='joint.new_reference_profile', verdict=verdict,
        new_mixed_gram_contractions=contractions, new_mixed_gram_terms=18*contractions,
        old_producers_replayed=0, new_inverse_or_solve_calls=0,
        new_matrix_or_determinant_calls=0, parameter_source_precision_scan=False,
        conditions=gap['conditions']+band['conditions']+[
            'The SAME actual transported majorant has unchanged full u, residual row B and charge K; no extra tau*h or shifted finite row is used.',
            'The saved small-shift first coordinate, mixed scalar and curvature describe the SAME Sigma0, bJ and exact17trial; its held-majorant verdict and penalty are not transferred.',
            'The raw G=lift.far_gram_upper is the SAME complete paid PSD Gram, not the energy-scaled joint Gram.',
            'The signed trial is reused only as a lower; its rank-one decrement uses the TRUE h interval.',
            'Strict positive endpoint lower excludes every smaller supplied shift by Loewner monotonicity for each D>0,D>=K; a nonpositive lower is inconclusive.'],
        scope='Conditional fixed-trial exclusion for the SAME fully transported inverse-majorant '
            'lower family on the supplied sufficient interval only. No actual Weil negativity '
            'or positivity, maximal capacity, good-block score, retained/cofinal positivity, '
            'Robin-tail, RH, Lean or generic novelty conclusion.',
        elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed: '+name)
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed: '+path.name)
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    for key, value in values.items():
        print(key, value.str(18), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.json':
        parser.error('Use a fresh .json output; accepted packets are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
