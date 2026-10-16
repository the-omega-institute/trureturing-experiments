"""Bound the actual same-source far energy with one signed Hankel Gram.

Use the cached q and exact w without recomputing an action or inverse.
The off-diagonal harmonic telescope is a linear Toeplitz convolution;
only the diagonal is padded by its integral bound. Cached smooth and
mixed-H remainders, charge identities and full-band budget are reused.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_poly, ctx, fmpq

from diagonal_floor_witness import interval_ball
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
FIRST, LAST, SPLIT = 1281, 32768, 2048
CHECK_INDICES = (0, 767, 31487)
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'full-gamma-dual-32768-finite.json.zlib':
        '4085fa7048bafc987c75b2340b56d98d556534ae74d304e61469d4e6e52af1b5',
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-common-column.json.zlib':
        '314e341189f72516aa9f251c6b873587f1bfae54d902c35f6c87c9b03130dfe7',
    'full-gamma-actual-channel-projection.json.zlib':
        '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'inverse-residual-pilot.json.zlib':
        '6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
}
PRODUCERS = ('full_gamma_dual_middle.py', 'full_gamma_dual_defect.py',
    'full_gamma_primal_tail.py', 'full_gamma_common_column.py',
    'full_gamma_actual_channel_projection.py', 'arb_inverse_residual_pilot.py',
    'full_gamma_band_surplus_budget.py', 'full_gamma_prime_potential_refined.py',
    'full_gamma_surplus_budget.py')


def principal_tail(z, first, cutoff, check_indices):
    """Real signed source; first omitted output mode is cutoff+1."""
    size = len(z)
    if not size or first <= 0 or first+size-1 > cutoff:
        raise ValueError('Require a nonempty source before the far cutoff.')
    if any(not value.is_finite() for value in z):
        raise ArithmeticError('Nonfinite source; stop without retry.')
    if any(i < 0 or i >= size for i in check_indices):
        raise ValueError('Invalid preselected convolution coordinate.')
    kernel = [arb(fmpq(1, j)) if j else arb(0)
              for j in range(1-size, size)]
    product = arb_poly(z)*arb_poly(kernel)
    kz = [product[size-1+i] for i in range(size)]
    diagonal, cross, h = arb(0), arb(0), arb(0)
    for i, value in enumerate(z):
        n = first+i
        diagonal += value*value/(cutoff+n)
        cross += 2*h*value*kz[i]
        h += arb(fmpq(1, cutoff+n+1))
    padded = diagonal+cross
    if not padded.is_finite() or padded.upper() < 0:
        raise ArithmeticError('Nonfinite or negative energy upper; stop without retry.')
    return dict(hankel_diagonal_padding=diagonal, hankel_signed_cross=cross,
                hankel_padded_energy=padded), {str(i): kz[i] for i in check_indices}


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    finite, defect, primal, common, actual, pilot, band, refined, budget = packets
    if [p['schema'] for p in packets] != [
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1', 'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_INVERSE_RESIDUAL_PILOT_V1', 'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1',
            'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1']:
        raise ValueError('Different cached packet roles.')
    if (any(p['metadata'] != finite['metadata'] for p in (defect, primal, common))
            or any(p['grid_bits'] != BITS for p in
                   (finite, defect, primal, common, actual, band, refined))
            or (finite['metadata']['N'], finite['metadata']['M'],
                finite['metadata']['last_mode'], finite['dual_first'], finite['dual_last'])
                != (FIRST-1, LAST, SPLIT, SPLIT+1, LAST)
            or finite['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different source support, basis, metadata or grid.')
    fpin = PINS['full-gamma-dual-32768-finite.json.zlib']
    dpin = PINS['full-gamma-dual-32768-defect.json.zlib']
    ppin = PINS['inverse-residual-pilot.json.zlib']
    if (any(p['finite_action_sha256'] != fpin or p['finite_defect_sha256'] != dpin
            or p['pilot_sha256'] != ppin for p in (primal, common))
            or common['primal_tail_sha256'] != PINS['full-gamma-primal-tail.json']
            or defect['finite_action_sha256'] != fpin
            or finite['metadata']['input_sha256']['inverse-residual-pilot.json.zlib'] != ppin
            or actual['packet_sha256']['common'] != PINS['full-gamma-common-column.json.zlib']
            or actual['packet_sha256']['defect'] != dpin):
        raise ValueError('Different actual source, action, far weight or channel.')
    for name in ('full-gamma-dual-32768-finite.json.zlib',
                 'full-gamma-dual-32768-defect.json.zlib', 'inverse-residual-pilot.json.zlib',
                 'full-gamma-common-column.json.zlib', 'full-gamma-primal-tail.json'):
        if budget['packet_sha256'][name] != PINS[name]:
            raise ValueError('Different inherited budget source: '+name)
    for name in ('full-gamma-prime-potential-refined.json', 'full-gamma-surplus-budget.json'):
        if band['packet_sha256'][name] != PINS[name]:
            raise ValueError('Different band weight or budget: '+name)
    integers = finite['dual_integers']
    identifier = hashlib.sha256(('\n'.join(str(int(x)) for x in integers)+'\n').encode()).hexdigest()
    if (finite['dual_bits'] != 384 or len(integers) != LAST-SPLIT or identifier != SOURCE_ID
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256'])
                   != (384, 0, SOURCE_ID) for p in (budget, band))
            or refined['source_integer_sha256'] != SOURCE_ID
            or any(len(defect[key]) != SPLIT-FIRST+1
                   or any(len(row) != 1 for row in defect[key])
                   for key in ('D_inverse_difference', 'F_inverse_z'))):
        raise ValueError('Different original exact source or saved inverse-profile shape.')
    producers = {}
    for packet, name in zip(packets, PRODUCERS):
        path = directory/name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != packet['source_sha256']:
            raise ValueError('Different cached producer source: '+name)
        producers[path] = digest
    z = [-(-1)**(FIRST+i)*interval_ball(row[0], BITS)
         for i, row in enumerate(defect['D_inverse_difference'])]
    z += [(-1)**(SPLIT+1+i)*arb(fmpq(int(x), 1 << 384)) for i, x in enumerate(integers)]
    scalars, selected = principal_tail(z, FIRST, LAST, CHECK_INDICES)
    far = fmpq(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator']))
    if not far > 0:
        raise ArithmeticError('Require the original strictly positive far weight.')
    upper = lambda packet, field, group: interval_ball(packet[group][field], BITS).upper()
    rw = upper(primal, 'gp_w_remainder_term', 'bounds')
    rq = upper(primal, 'gp_q_remainder_term', 'bounds')
    hmixed = upper(primal, 'actual_mixed_tail_norm_upper', 'bounds')
    if any(not v.is_finite() or v < 0 for v in (rw, rq, hmixed)):
        raise ArithmeticError('Invalid cached remainder upper.')
    principal_norm = arb(far).sqrt()*scalars['hankel_padded_energy'].upper().sqrt()/2
    remainder_norm = arb(far).sqrt()*(rw+rq)
    full_norm = (principal_norm+remainder_norm+hmixed).upper()
    energy = (full_norm*full_norm).upper()
    old_energy = upper(common, 'far_action_energy_upper', 'bounds')
    hlo = interval_ball(defect['denominator'], BITS).lower()
    channel_lo = interval_ball(actual['scalar_bounds']['actual_channel_energy'], BITS).lower()
    k0 = (upper(common, 'middle_action_energy', 'bounds')-hlo-channel_lo).upper()
    fee = (k0+arb(fmpq(33, 32))*energy).upper()
    old_fee = upper(actual, 'refined_common_charge_upper', 'scalar_bounds')
    if not fee.is_finite() or not fee > 0:
        raise ArithmeticError('No finite positive source-valid fee; stop without retry.')
    unweighted = (energy/arb(far)).upper()
    a1 = upper(refined, 'finite_weighted_first_moment', 'scalar_bounds')
    b2 = upper(band, 'high_band_weighted_source_norm_squared_upper', 'scalar_bounds')
    maximum = fmpq(*map(int, band['exact_rationals']['maximum_weight']))
    h_exact = fmpq(*map(int, budget['exact_rationals']['h_lower']))
    old_cap = fmpq(*map(int, band['exact_rationals']['supplied_parameter_cap']))
    if not b2 >= 0 or not a1 > 0 or not maximum > 0 or not h_exact > 0:
        raise ArithmeticError('Invalid saved full-band budget components.')
    weighted_budget = (a1+2*(unweighted*b2).sqrt()+arb(maximum)*unweighted).upper()
    budget_pair = envelope(weighted_budget, BITS)
    if budget_pair is None or not weighted_budget > 0:
        raise ArithmeticError('No finite positive full budget; stop without retry.')
    budget_upper = fmpq(int(budget_pair[1]), 1 << BITS)
    cap = h_exact/(h_exact+budget_upper)
    if not 0 < cap < 1 or not 0 < old_cap < 1:
        raise ArithmeticError('Invalid supplied eligibility cap; stop without retry.')
    scalars.update(principal_weighted_norm_upper=principal_norm,
        smooth_remainder_weighted_norm_upper=remainder_norm, mixed_H_weighted_norm_upper=hmixed,
        full_far_weighted_norm_upper=full_norm, far_energy_upper=energy,
        old_far_energy_upper=old_energy, far_energy_improvement_lower=(old_energy-energy).lower(),
        charge_constant_upper=k0, source_valid_fee_upper=fee, old_source_valid_fee_upper=old_fee,
        fee_improvement_lower=(old_fee-fee).lower(), far_unweighted_energy_upper=unweighted,
        complete_weighted_budget_upper=weighted_budget)
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    checks_export = {key: envelope(value, BITS) for key, value in selected.items()}
    if any(pair is None or int(pair[0]) > int(pair[1])
           for pair in list(exported.values())+list(checks_export.values())):
        raise ArithmeticError('Invalid outward export; stop without retry.')
    improved = energy < old_energy
    fee_improved, cap_improved = fee < old_fee, cap > old_cap
    result = dict(schema='ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={p.name: h for p, h in guards.items()},
        cached_producer_source_sha256={p.name: h for p, h in producers.items()},
        metadata=finite['metadata'], source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        first_source_mode=FIRST, split_source_mode=SPLIT, last_source_mode=LAST,
        first_far_mode=LAST+1, grid_bits=BITS, working_bits=1024,
        exact_rationals={key: [str(v.p), str(v.q)] for key, v in dict(
            original_far_weight=far, kappa=fmpq(33, 32), maximum_reflected_weight=maximum,
            h_lower=h_exact, supplied_parameter_cap=cap, old_supplied_parameter_cap=old_cap).items()},
        scalar_bounds=exported, preselected_convolution_indices=list(CHECK_INDICES),
        selected_convolution_intervals=checks_export,
        strict_far_energy_improvement=bool(improved), strict_fee_improvement=bool(fee_improved),
        strict_supplied_cap_improvement=bool(cap_improved), positive_source_valid_fee=True,
        verdict=('conditional-same-source-far-energy-improved' if improved
                 else 'inconclusive-fixed-signed-hankel-tail-recipe'),
        new_polynomial_products=1, new_joint_source_coefficients=len(z),
        old_producers_replayed=0, new_inverse_or_solve_calls=0, new_far_modes_generated=0,
        new_weight_projection_or_determinant_calls=0, parameter_source_precision_scan=False,
        elapsed_seconds=time.monotonic()-started,
        conditions=primal['conditions']+actual['conditions']+band['conditions']+[
            'The real SAME-source action is uF=G_FJ*w-G_FE*q+H_FE*(q-p), with the original constant positive Wfar.',
            'The signed principal Gram retains all E/J cross terms; only its diagonal is padded, and its energy normalization is1/4.',
            'Saved gp_w/q_remainder_term bounds are unweighted and already include gamma; the mixed-H bound is already weighted.',
            'The source charge is K=middle-h+(33/32)*s-Ephi. Its independent constant is rebuilt using the channel and h LOWER endpoints.',
            'The saved full-band Bessel projection applies unchanged to the SAME uJ and reflected potential weight; an improved cap supplies eligibility only.'],
        scope='Numerical application of classical partial fractions, diagonal integral bounds, '
            'FLINT linear convolution, Minkowski and existing charge/Bessel budget identities. '
            'No new generic theorem or Lean verification. The complete signed far row, joint good-block '
            'metric, actual retained/cofinal positivity, Robin tail and RH remain unresolved.')
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed: '+name)
    for path, digest in {**guards, **producers}.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked helper or cached producer changed: '+path.name)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':')).encode(), 9))
    for key, value in scalars.items():
        print(key, value.str(20), flush=True)
    print('supplied cap', arb(cap).str(20), 'old cap', arb(old_cap).str(20), flush=True)
    print(result['verdict'], 'fee improved', fee_improved, 'cap improved', cap_improved,
          'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh .zlib output; accepted archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
