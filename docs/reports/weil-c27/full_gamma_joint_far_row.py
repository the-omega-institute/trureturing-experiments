"""Supply the complete fixed SAME-source far row from accepted caches.

Reuse signed coordinate0; evaluate each remaining original lift column
once. Apply the existing finite form and square completion to the joint
row. No source recovery/interpolation, inverse/solve or direction scan.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_mat, arb_poly, ctx, fmpq
from diagonal_floor_witness import interval_ball
from full_gamma_cached_source_response import finite_pairing, cached_tail
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_signed_coordinate import coordinate_errors, PINS as COORDINATE_PINS
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from full_gamma_zero_channel_jet_refine import scaled_forward
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, R, SIZE = 640, 32768, 56, 18
PINS = dict(COORDINATE_PINS)
PINS.update({'full-gamma-signed-coordinate.json.zlib': '12cd1ff9de20ab9083d27ddae812d415ec37c95520de619a6f9813304c2081a5', 'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd', 'full-gamma-directional-metric.json.zlib': 'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2', 'full-gamma-projection-upper.json.zlib': 'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3', 'full-gamma-actual-channel-projection.json.zlib': '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745', 'full-gamma-joint-plane-paid-floor.json.zlib': '3bdf5f3f81daeeda99c5fa71f2a48ee64c84c09da0f20070b45b0ad781aeec22'})


def joint_row_forms(f, sigma, reference, rho, y0, gamma, b, g, h):
    """Evaluate new SAME-row forms; common-form identities are premises."""
    n = f.ncols()
    if (n < 3 or f.nrows() != 1 or sigma.nrows() != n or sigma.ncols() != n
            or (reference.nrows(), reference.ncols()) != (n-1, n-1)
            or (rho.nrows(), rho.ncols(), y0.nrows(), y0.ncols()) != (1, n, 1, n)):
        raise ValueError('Incompatible common row, core or reference dimensions.')
    require_finite(list(f.entries())+list(sigma.entries())+list(reference.entries())
        +list(rho.entries())+list(y0.entries())+[arb(x) for x in (gamma, b, g, h)])
    delta = -sigma[0, 0]
    a = h+g*g/b if b > 0 else arb(0)
    if not delta > 0 or not 0 < gamma <= 1 or not b > 0 or not a > 0:
        raise ValueError('Missing negative-pivot/common-plane domains.')
    # Same accepted elimination map; no inverse or old reference quadratic.
    T = arb_mat([[-sigma[0, j+1]/sigma[0, 0] if i == 0
                  else arb(int(i == j+1)) for j in range(n-1)] for i in range(n)])
    z = f*T
    action = reference*z.transpose()
    pff = (z*action)[0, 0]
    pfy = ((y0*T)*action)[0, 0]
    pfr = ((rho*T)*action)[0, 0]
    qff = arb(gamma)*pff-f[0, 0]**2/delta
    theta = arb(gamma)*pfy-f[0, 0]*y0[0, 0]/delta
    response = arb(gamma)*pfr-f[0, 0]*rho[0, 0]/delta
    s = response+(g/b)*theta
    residual = qff+theta**2/b-s**2/a
    result = dict(positive_row_energy=pff, row_Qeps_energy=qff,
        direct_Qeps_timelike_pairing=theta, direct_Qeps_rho_pairing=response,
        positive_plane_energy=a, direct_positive_plane_pairing=s,
        complete_residual_energy=residual)
    require_finite(list(result.values()))
    return result


def direct_paid_floor(forms, alpha0, b, g, h, k, fee):
    """Lower the all-real steering minimum of the same comparison form."""
    require_finite(list(forms.values())+[arb(x) for x in (alpha0, b, g, h, k, fee)])
    curvature = arb(k)+h
    if not k > 0 or not fee > 0 or not b > 0 or not curvature > 0:
        raise ValueError('Missing SAME-recipe positive curvature/fee.')
    slope = arb(k)*arb(alpha0)+g+forms['direct_Qeps_rho_pairing']
    constant = arb(k)*arb(alpha0)**2+arb(fee)-b \
        +2*forms['direct_Qeps_timelike_pairing']+forms['row_Qeps_energy']
    minimum = constant-slope**2/curvature
    require_finite([curvature, slope, constant, minimum])
    return dict(curvature=curvature, slope=slope, constant=constant,
                comparison_minimum=minimum), freeze_lower(minimum.lower())


def consume(out):
    started, source = time.perf_counter(), Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted joint-row input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    inherited = dict(PINS)
    for packet in packets.values():
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256',
                      'cached_producer_source_sha256'):
            for name, digest in packet.get(field, {}).items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting joint-row input binding: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        p = HERE/name
        if p in guards and guards[p] != digest:
            raise ValueError('Different invoked joint-row dependency: '+name)
        guards[p] = digest
    verify_guards(guards)
    get = lambda name: packets[name]
    cached = get('full-gamma-cached-source-response.json.zlib')
    cauchy = get('full-gamma-signed-mu-direct-sine.json.cauchy.zlib')
    mu = get('full-gamma-signed-mu-direct-sine.json.zlib')
    seed = get('tail-moments.json.zlib')
    lift = get('full-gamma-lift-projection.json.zlib')
    signed = get('full-gamma-signed-correction.json.zlib')
    old = get('full-gamma-negative-channel-capacity.json.zlib')
    hankel = get('full-gamma-hankel-tail-energy.json.zlib')
    first = get('full-gamma-signed-coordinate.json.zlib')
    core = get('joint-gamma-6765-rho32-v2-schur.json')
    cert = get('full-gamma-directional-metric.json.zlib')
    projection = get('full-gamma-projection-upper.json.zlib')
    actual = get('full-gamma-actual-channel-projection.json.zlib')
    plane = get('full-gamma-joint-plane-paid-floor.json.zlib')
    if (any(p['grid_bits'] != BITS for p in packets.values())
            or first['schema'] != 'ACTUAL_C27_SIGNED_RAW_COORDINATE_V1'
            or first['source_sha256'] != guards[HERE/'full_gamma_signed_coordinate.py']
            or cached['schema'] != 'ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1'
            or core['schema'] != 'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1'
            or cert['schema'] != 'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1'
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256']) !=
                   (384, 0, first['source_integer_sha256'])
                   for p in (cached, cauchy, mu, signed, old, hankel, plane))
            or projection['exact_source_integer_sha256'][0] != first['source_integer_sha256']
            or projection['exact_source_bits'][0] != 384
            or any(p['remaining_modes'] != lift['remaining_modes']
                   or p['input_sha256'] != lift['input_sha256']
                   for p in (core, cert, projection, actual, old))
            or len(lift['remaining_modes']) != SIZE
            or (lift['N'], lift['M'], lift['last_mode']) != (1280, M, 2048)
            or lift['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or cached['recovered_symbols_are_actual'] or not cached['inherited_symbol_error_paid']
            or seed['analytic_grid_bits'] != mu['symbol_metadata_grid_bits']
            or (seed['c'], seed['M']) != (27, M)
            or (cached['finite_first_mode'], cached['finite_last_mode'],
                cached['analytic_first_mode'], cached['normalization_scale']) != (M+1, 2*M, 2*M+1, M)
            or any(first['packet_sha256'][name] != COORDINATE_PINS[name]
                   for name in COORDINATE_PINS)
            or not first['positive_residual_lower_certified']):
        raise ValueError('Different SAME-source/lift/form or cached error contract.')
    shapes = ((lift['far_features_on_lift'], 112, SIZE), (lift['lift_gram'], SIZE, SIZE),
        (lift['far_gram_upper'], SIZE, SIZE), (core['reduced_schur_intervals'], SIZE, SIZE),
        (cert['reference_inverse'], SIZE-1, SIZE-1),
        (cached['physical_left_kernel_action'], 1, 112), (actual['actual_finite_row_on_lift'], 1, SIZE))
    if (any(len(rows) != n or any(len(row) != m for row in rows) for rows,n,m in shapes)
            or len(projection['matrices']['source_on_lift'][0]) != SIZE
            or len(cached['recovered_approximate_symbols']) != M
            or len(cauchy['finite_cauchy_source_outputs']) != M
            or not 1 <= len(cached['interpolated_source_coefficients']) <= 13
            or len(seed['normalized_polynomial']) != 64
            or any(str(q) not in mu['real_hurwitz_kernels'] for q in range(2,188))):
        raise ValueError('Incomplete fixed joint-row input dimensions.')
    rational = lambda p,key: fmpq(*map(int,p['exact_rationals'][key]))
    ball = lambda p,key: interval_ball(p['scalar_bounds'][key], BITS)
    upper = lambda p,key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    kappa, weight = rational(cached,'kappa'), rational(cached,'original_far_weight')
    alpha0, gamma = rational(plane,'alpha0'), rational(plane,'gamma')
    eta = upper(old,'neumann_eta_upper')
    if (kappa != fmpq(33,32) or kappa != rational(signed,'kappa')
            or kappa != rational(mu,'kappa') or not weight > 0
            or weight != rational(signed,'original_far_weight')
            or weight != rational(mu,'original_far_weight')
            or not 0 <= eta < 1 or gamma != 1/(1+eta)
            or alpha0 != rational(old,'fixed_alpha')):
        raise ValueError('Different raw weights, reference or comparison parameters.')
    # Saved scalars only: ceiling of coordinate-only LOWER certification.
    cross_cap = max(abs(rational(first,'cross_lower')),abs(rational(first,'cross_upper')))
    tlo = rational(first,'projected_norm_lower')
    if not tlo > 0:
        raise ValueError('Missing positive raw coordinate dual norm lower.')
    coordinate_ceiling = cross_cap**2/tlo
    exact = dict(coordinate_only_residual_certificate_ceiling=coordinate_ceiling)
    # First coordinate is REUSED, not evaluated again.
    flo, fhi = rational(first,'coordinate_lower'), rational(first,'coordinate_upper')
    if flo > fhi:
        raise ValueError('Unordered saved raw coordinate.')
    row = [arb(arb((flo+fhi)/2), arb((fhi-flo)/2))]
    samples = [interval_ball(v,BITS) for v in cauchy['finite_cauchy_source_outputs']]
    symbols = [interval_ball(v,BITS) for v in cached['recovered_approximate_symbols']]
    polynomial = arb_poly([interval_ball(v,BITS) for v in cached['interpolated_source_coefficients']])
    smooth = arb_poly([interval_ball(v,BITS) for v in seed['normalized_polynomial']])
    kernels = {q:interval_ball(v,BITS) for q,v in mu['real_hurwitz_kernels'].items()}
    feature_matrix = matrix_in([r[1:] for r in lift['far_features_on_lift']],BITS)
    mixed_row = arb(kappa)*(matrix_in(cached['physical_left_kernel_action'],BITS)*feature_matrix)
    old_lift = rational(signed,'lift_norm_upper')
    old_forward = rational(signed,'sharp_forward_error_upper')
    principal = upper(hankel,'principal_weighted_norm_upper')
    smooth_cap = upper(hankel,'smooth_remainder_weighted_norm_upper')
    sigma_error = interval_ball(seed['uniform_symbol_error'],seed['analytic_grid_bits']).upper()
    details = []
    for j in range(1,SIZE):
        features = [feature_matrix[i,j-1] for i in range(112)]
        lcap = freeze_upper(interval_ball(lift['lift_gram'][j][j],BITS).upper().sqrt())
        ncap2 = fmpq(int(lift['far_gram_upper'][j][j][1]),1 << BITS)
        if not ncap2 >= 0:
            raise ArithmeticError('Negative physical Gram diagonal cap; no retry.')
        ncap = freeze_upper(arb(ncap2).sqrt())
        forward = scaled_forward(old_forward,old_lift,lcap)
        mixed_error = freeze_upper(arb(rational(signed,'mixed_cross_error_upper'))*arb(lcap)/arb(old_lift))
        errors = coordinate_errors(features[:R],sigma_error,weight,
            rational(cached,'source_tail_norm_error_upper'),ncap,forward,
            principal,smooth_cap,mixed_error,kappa)
        finite = finite_pairing(samples,symbols,features,M,kappa,weight)
        tail, tail_cap, oscillatory = cached_tail(polynomial,features,smooth,kernels,
            upper(mu,'inverse_phase_gap_weight_cap'),M,kappa,weight)
        allowance = errors['nonoscillatory_error_upper']+oscillatory
        value = finite+tail+mixed_row[0,j-1]+arb(0,arb(allowance))
        require_finite([value])
        row.append(value)
        details.append(dict(coordinate=j, raw_coordinate=envelope(value,BITS),
            finite=envelope(finite,BITS), tail=envelope(tail,BITS),
            mixed=envelope(mixed_row[0,j-1],BITS),
            exact_rationals={key:[str(v.p),str(v.q)] for key,v in dict(errors,
                lift_norm_upper=lcap,physical_norm2_upper=ncap2,
                sharp_forward_error_upper=forward,oscillatory_error_upper=oscillatory,
                complete_error_upper=allowance).items()}))
    f = arb_mat([row])
    raw_sigma = matrix_in(core['reduced_schur_intervals'],BITS)
    sigma = (raw_sigma+raw_sigma.transpose())/2
    reference = arb_mat([[arb(fmpq(v)) for v in r] for r in cert['reference_inverse']])
    rho = matrix_in([projection['matrices']['source_on_lift'][0]],BITS)
    beta = matrix_in(actual['actual_finite_row_on_lift'],BITS)
    y0 = arb(alpha0)*rho+beta
    b,g,h = [ball(plane,key) for key in ('timelike_energy','rho_timelike_pairing','rho_reference_energy')]
    forms = joint_row_forms(f,sigma,reference,rho,y0,gamma,b,g,h)
    if (not forms['direct_Qeps_timelike_pairing'].overlaps(ball(cached,'fixed_transported_theta'))
            or not forms['positive_row_energy'].lower() <= rational(plane,'fixed_far_positive_energy_upper')
            or not forms['positive_plane_energy'].overlaps(ball(cached,'fixed_positive_plane_energy'))):
        raise ArithmeticError('New joint row conflicts with SAME-form saved observations; stop.')
    rlo = max(fmpq(0),rational(first,'residual_lower'),
              freeze_lower(forms['complete_residual_energy'].lower()))
    rhi = freeze_upper(forms['complete_residual_energy'].upper())
    if rlo > rhi:
        raise ArithmeticError('Joint residual contradicts accepted coordinate residual; stop.')
    exact.update(full_residual_lower=rlo,full_residual_upper=rhi)
    recipes = {}
    for recipe in ('fixed','energy_endpoint'):
        k,fee = [rational(plane,recipe+'_'+key) for key in ('reciprocal_coefficient','comparison_fee')]
        lower_F = rational(cached,recipe+'_all_alpha_floor_lower')
        previous = rational(first,recipe+'_retained_lower')
        minimum,direct_lower = direct_paid_floor(forms,alpha0,b,g,h,k,fee)
        projected_lower = lower_F+rlo
        retained = max(previous,projected_lower,direct_lower)
        ceiling = max(previous,lower_F+coordinate_ceiling)
        exact.update({recipe+'_coordinate_only_retained_certificate_ceiling':ceiling,
            recipe+'_saved_same_F_lower':lower_F,recipe+'_previous_retained_lower':previous,
            recipe+'_full_residual_augmented_lower':projected_lower,
            recipe+'_direct_all_alpha_lower':direct_lower,recipe+'_retained_lower':retained})
        forms.update({recipe+'_'+key:value for key,value in minimum.items()})
        recipes[recipe] = dict(coordinate_only_refinement_cannot_certify_positive=bool(ceiling < 0),
            strictly_improved=bool(retained > previous),positive_paid_lower_certified=bool(retained > 0))
    exported = {key:envelope(v,BITS) for key,v in forms.items()}
    if any(v is None or int(v[0]) > int(v[1]) for v in exported.values()):
        raise ArithmeticError('Invalid joint-form export; no retry.')
    result = dict(schema='ACTUAL_C27_JOINT_FAR_ROW_V1',packet_sha256=PINS,
        source_sha256=guards[source],inherited_dependency_sha256={p.name:h for p,h in guards.items() if p != source},
        source_integer_sha256=first['source_integer_sha256'],source_bits=384,source_row=0,
        grid_bits=BITS,working_bits=ctx.prec,remaining_modes=lift['remaining_modes'],
        actual_far_row=matrix_out(f,BITS),new_coordinate_details=details,
        scalar_bounds=exported,exact_rationals={k:[str(v.p),str(v.q)] for k,v in exact.items()},
        recipes=recipes,new_original_coordinate_columns=list(range(1,SIZE)),reused_coordinate0=True,
        new_finite_functionals=17,new_cached_tail_applications=17,new_mixed_batch_products=1,
        new_reference_actions=1,new_reference_scalar_contractions=3,new_direct_paid_minima=2,
        new_old_reference_quadratics=0,new_source_recovery_interpolation_kernel_row_inverse_solve=0,
        old_producers_or_fixtures_replayed=0,source_direction_parameter_precision_cutoff_scan=False,
        actual_center_sign_claim=False,actual_full_border_sign_claim=False,actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'Same source/lift/raw row f=kappa*u_F^*W_F*R_F*V, kappa33/32. Coordinate0 is reused; every remaining original basis coordinate is supplied once, not selected or scanned.',
            'Saved FULL physical Gram Loewner majorant bounds each ||B e_j||_W^2, with original W and no kappa. Uniform forward(actual-symbol jet) and raw mixed errors scale by the SAME lift norm cap; all symbol/source/Abel/smooth/forward/mixed allowances are paid.',
            'Coordinate intervals form a JOINT containing box for ONE source. Shared errors need not be independent; separate endpoints are never asserted jointly attainable.',
            'The original old core, reference inverse and finite rows give P=T R T*, Qeps=gamma P-e0e0*/delta <= Qtrue. Saved b,g,h and row data retain one common realization; overlap checks only detect inconsistency, not prove these identities.',
            'The complete residual is Qeps(f,f)+theta^2/b-s^2/a. Combine its lower with the previous coordinate lower by max, not addition. Add only to the saved SAME-F lower before retention.',
            'Direct square completion uses the SAME two paid recipes and all-real steering. A positive retained lower certifies only those paid comparisons for the fixed source under inherited premises, not smaller actual-center/full-border/cofinal positivity, Robin or RH.',
            'The coordinate-only ceiling bounds possible strength of a LOWER certificate with other inputs fixed; it is not an upper bound on actual residual or energy. A nonpositive numerical lower is inconclusive about actual sign. Existing methods, no generic novelty or Lean claim.',
        ])
    verify_guards(guards)
    digest = write_packet(out,result)
    print('saved complete joint far row',digest,'seconds',result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; never overwrite accepted data.')
    ctx.prec = 1024
    consume(args.out)
