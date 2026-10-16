"""Project one saved coherent direction onto the complete short residual jet.

Price a forward jet error and the same-source principal-proxy error before
supplying any long-source signed moment. Reuse all old certificates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, ctx, fmpq

from diagonal_floor_witness import interval_ball
from full_gamma_coherent_split_capacity import freeze_upper
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-coherent-mixed-capacity.json.zlib':
        'b082eaad102c5f7ca38d1c49c19970952288ce630198a4ea76b48b96d8d23b65',
    'full-gamma-coherent-split-capacity.json.zlib':
        '6e4232d6079a4279511a67c13c6501d51a3b4fdaf410ccfee957adce29e1d7bc',
    'full-gamma-lift-projection.json.zlib':
        'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
    'far-coupling.json.zlib':
        'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7',
}


def project_and_price(features, direction, lift_diagonal, allowance,
                      ncap, principal, proxy_error, kappa, clin):
    """Enclose the ONE common direction; no Gram-error-to-vector shortcut."""
    count = direction.nrows()
    if (direction.ncols() != 1 or features.ncols() != count or count < 1
            or features.nrows() < 2 or features.nrows() % 2
            or len(lift_diagonal) != count):
        raise ValueError('Incompatible short-jet and lifted direction dimensions.')
    values = list(features.entries())+list(direction.entries())+list(lift_diagonal)
    values += [arb(v) for v in (allowance, ncap, principal, proxy_error, kappa, clin)]
    if any(not v.is_finite() for v in values):
        raise ArithmeticError('Nonfinite interface input; stop without retry.')
    if (not allowance >= 0 or not ncap >= 0 or not principal >= 0
            or not proxy_error >= 0 or not kappa > 0
            or any(v.upper() < 0 for v in lift_diagonal)):
        raise ArithmeticError('Missing positive scalar or lift-norm bound.')
    h = features*direction
    if any(not v.is_finite() for v in h.entries()):
        raise ArithmeticError('Nonfinite short-jet coefficients; stop without retry.')
    lift_l1 = sum((abs(direction[j, 0]).upper()*lift_diagonal[j].upper().sqrt().upper()
                   for j in range(count)), arb(0)).upper()
    lift_upper = freeze_upper(lift_l1)
    jet_error = freeze_upper((2*arb(allowance)).sqrt().upper()*arb(lift_upper))
    error = freeze_upper(arb(kappa)*(arb(proxy_error)*arb(ncap)
                                    +arb(principal)*arb(jet_error)))
    threshold = -clin/2-error
    scalars = dict(lift_norm_l1_upper=arb(lift_upper),
        cached_complete_far_allowance=arb(allowance),
        jet_forward_error_upper=arb(jet_error),
        cached_exact_direction_norm_upper=arb(ncap),
        cached_principal_proxy_norm_upper=arb(principal),
        principal_proxy_error_upper=arb(proxy_error),
        signed_proxy_functional_error_upper=arb(error),
        frozen_linear_envelope_upper=arb(clin),
        strict_theta_upper_threshold=arb(-clin/2),
        strict_mu_upper_threshold=arb(threshold))
    if any(not v.is_finite() for v in scalars.values()):
        raise ArithmeticError('Nonfinite priced error; stop without retry.')
    return h, scalars, lift_upper, jet_error, error, threshold


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted interface packet: '+name)
        guards[directory/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    mixed = packets['full-gamma-coherent-mixed-capacity.json.zlib']
    split = packets['full-gamma-coherent-split-capacity.json.zlib']
    lift = packets['full-gamma-lift-projection.json.zlib']
    proxy = packets['full-gamma-hankel-tail-energy.json.zlib']
    original = packets['structured-joint-2048-target.json.zlib']
    far = packets['far-coupling.json.zlib']
    if [p['schema'] for p in (mixed, split, lift, proxy, original, far)] != [
            'ACTUAL_C27_COHERENT_MIXED_CAPACITY_V1',
            'ACTUAL_C27_COHERENT_SPLIT_CAPACITY_V1',
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'ACTUAL_C27_STRUCTURED_TARGET_V1', 'ACTUAL_C27_FAR_COUPLING_GRAM_V1']:
        raise ValueError('Different accepted packet roles.')
    for key in ('source_integer_sha256', 'source_bits', 'source_row'):
        if not mixed[key] == split[key] == proxy[key]:
            raise ValueError('Different same-source identity: '+key)
    if (any(p['grid_bits'] != BITS for p in (mixed, split, lift, proxy, original, far))
            or not mixed['remaining_modes'] == split['remaining_modes'] == lift['remaining_modes']
            or len(mixed['remaining_modes']) != 18
            or not mixed['input_sha256'] == split['input_sha256'] == lift['input_sha256']
                == original['metadata']['input_sha256']
            or lift['original_target_sha256'] != PINS['structured-joint-2048-target.json.zlib']
            or (lift['N'], lift['M'], lift['last_mode'], far['jet_order'], far['moment_diagonals'])
                != (1280, 32768, 2048, 56, 56)
            or lift['basis'] != original['metadata']['basis']
            or lift['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or proxy['first_far_mode'] != 32769
            or proxy['split_source_mode'] != 2048 or proxy['last_source_mode'] != 32768
            or (original['metadata']['N'], original['metadata']['M'],
                original['metadata']['last_mode']) != (1280, 32768, 2048)
            or (far['c'], far['N'], far['M']) != (27, 1280, 32768)
            or (mixed['source_bits'], mixed['source_row']) != (384, 0)):
        raise ValueError('Different common lift, physical support, basis or grid.')
    inherited = dict(split['inherited_dependency_sha256'])
    for packet in (split, mixed, proxy):
        for field in ('invocation_source_sha256', 'cached_producer_source_sha256',
                      'packet_sha256'):
            for name, digest in packet.get(field, {}).items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting inherited pin: '+name)
                inherited[name] = digest
    producer_pins = {
        'full_gamma_lift_projection.py': lift['source_sha256'],
        'structured_inverse_residual_interval.py': original['producer_source_sha256'],
        'arb_far_coupling.py': far['source_sha256'],
        'arb_tail_moments.py': far['moment_supplier_sha256'],
    }
    for name, digest in producer_pins.items():
        if name in inherited and inherited[name] != digest:
            raise ValueError('Different producer pin: '+name)
        inherited[name] = digest
    for name, digest in inherited.items():
        path = directory/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Different accepted read dependency: '+name)
        if path in guards and guards[path] != digest:
            raise ValueError('Different invocation dependency: '+name)
        guards[path] = digest
    if any(mixed['exact_rationals'][k] != v for k, v in split['exact_rationals'].items()
           if k in mixed['exact_rationals']):
        raise ValueError('Different frozen common parameters.')
    exact = {k: fmpq(*map(int, split['exact_rationals'][k])) for k in
             ('frozen_linear_envelope_upper', 'strict_theta_upper_threshold',
              'selected_direction_norm_upper', 'kappa')}
    clin = exact['frozen_linear_envelope_upper']
    if (exact['strict_theta_upper_threshold'] != -clin/2
            or exact['kappa'] != fmpq(33, 32)
            or exact['kappa'] != fmpq(*map(int, proxy['exact_rationals']['kappa']))
            or fmpq(*map(int, proxy['exact_rationals']['original_far_weight']))
                != fmpq(int(far['D_numerator']), 1 << far['D_grid_bits'])):
        raise ValueError('Different same-source normalization or signed threshold.')
    allowance = fmpq(lift['exact_far_allowance'])
    # Metadata allowances retain their ORIGINAL512 grid; matrix arrays use640.
    metadata_bits = lift['original_metadata_grid_bits']
    original_allowance = fmpq(int(original['metadata']['complete_far_operator_allowance'][1]),
                              1 << metadata_bits)
    if metadata_bits != 512 or not allowance >= original_allowance >= 0:
        raise ValueError('Missing original complete allowance or metadata grid.')
    features = matrix_in(lift['far_features_on_lift'], BITS)
    direction = matrix_in(mixed['reference_direction_on_lift'], BITS)
    if (features.nrows(), features.ncols(), direction.nrows(), direction.ncols()) != (112, 18, 18, 1):
        raise ValueError('Incomplete full56 feature direction.')
    gram_rows = lift['lift_gram']
    if len(gram_rows) != 18 or any(len(row) != 18 for row in gram_rows):
        raise ValueError('Incomplete saved lift diagonal.')
    diagonal = [interval_ball(gram_rows[j][j], BITS) for j in range(18)]
    upper = lambda key: interval_ball(proxy['scalar_bounds'][key], BITS).upper()
    principal = upper('principal_weighted_norm_upper')
    proxy_error = freeze_upper(upper('smooth_remainder_weighted_norm_upper')
                               +upper('mixed_H_weighted_norm_upper'))
    h, scalars, l1, eps_b, error, threshold = project_and_price(
        features, direction, diagonal, allowance, exact['selected_direction_norm_upper'],
        principal, proxy_error, exact['kappa'], clin)
    exports = {k: envelope(v, BITS) for k, v in scalars.items()}
    coefficients = matrix_out(h, BITS)
    if (any(p is None or int(p[0]) > int(p[1]) for p in exports.values())
            or any(p is None or int(p[0]) > int(p[1]) for row in coefficients for p in row)):
        raise ArithmeticError('Invalid forward-interface exports; stop without retry.')
    result = dict(schema='ACTUAL_C27_SIGNED_JET_INTERFACE_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: h for p, h in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: h for p, h in guards.items() if p not in invoked},
        input_sha256=mixed['input_sha256'], remaining_modes=mixed['remaining_modes'],
        source_integer_sha256=mixed['source_integer_sha256'], source_bits=384, source_row=0,
        physical_last_mode=2048, first_far_mode=32769, jet_order=56,
        feature_order='A0..A55,B0..B55', grid_bits=BITS, working_bits=1024,
        original_metadata_grid_bits=metadata_bits, short_jet_coefficients=coefficients,
        exact_rationals=split['exact_rationals'] | {
            'complete_far_allowance': [str(allowance.p), str(allowance.q)],
            'lift_norm_l1_upper': [str(l1.p), str(l1.q)],
            'jet_forward_error_upper': [str(eps_b.p), str(eps_b.q)],
            'principal_proxy_error_upper': [str(proxy_error.p), str(proxy_error.q)],
            'signed_proxy_functional_error_upper': [str(error.p), str(error.q)],
            'strict_mu_upper_threshold': [str(threshold.p), str(threshold.q)]},
        scalar_bounds=exports, verdict='conditional-signed-jet-interface-only',
        actual_mu_supplied=False, actual_signed_theta_supplied=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False,
        new_cached_feature_actions=1, new_lift_diagonal_scalar_terms=18,
        new_reference_vector_actions=0, new_reference_quadratic_forms=0,
        new_full_far_gram_forms=0, new_trace_contractions=0, new_inverse_or_solve_calls=0,
        new_source_or_far_suppliers=0, old_producers_replayed=0,
        parameter_source_precision_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=split['conditions']+[
            'The full56 jet uses exact2/pi and actual same s_m. Its features are the fixed P0 midpoints transported through the SAME exact physical2048 residual and lift.',
            'Forward error uses the original geometric and feature-rounding component bounds: each squared operator error is covered by its nonnegative allowance term. Gram closeness alone would not imply it.',
            'sqrt(2eps_complete)*||VvR|| pays both component errors; eps_complete retains residual amplification and its original512-grid allowance. All18 diagonal terms bound the ONE common direction.',
            'The full PSD moment Gram controls the feature error; the sparse masked moment core is not assumed PSD and is not used as a whole-jet norm.',
            'Same u0 uses minus q on E and original w on J. Its independent principal norm cap and smooth/mixedH residual payments remain valid; no q remainder is renamed a p remainder.',
            'theta-mu=kappa*(<uF-u0,BvR>+<u0,BvR-j>) gives error<=kappa*(eps_u*Ncap+P*eps_b). All saved coefficient intervals must propagate; no midpoint for free.',
            'The exact future condition is UP(Clin+2*(I_mu+error))<0; a priced threshold without ACTUAL same-source I_mu supplies no sign, favorable realization or expensive-supplier justification. The good-block metric remains unpaid.'
        ], scope='Source-specific cached forward-jet and principal-proxy error interface using existing bounds. No generic method or Lean claim. Actual mu/theta, center/full-border sign, good-block metric, retained/cofinal positivity, Robin tail and RH remain unresolved.')
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError('Read dependency changed: '+str(path))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':'), allow_nan=False).encode(), 9))
    for key, value in scalars.items():
        print(key, value.str(20), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


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
