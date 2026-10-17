"""Test ONE saved-row midpoint through a new common-source functional.

The fixed exact covector is Q0 times the saved complete-row midpoint.
Aggregate Hc before physical error bounds. Reuse the existing signed
functional, jet-norm and residual projection interfaces; no scan or solve.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_mat, arb_poly, ctx, fmpq, fmpq_mat
from diagonal_floor_witness import interval_ball
from full_gamma_cached_source_response import finite_pairing, cached_tail
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_directional_jet_norm import directional_norm
from full_gamma_joint_far_row import PINS as FULL_ROW_PINS
from full_gamma_mixed_norm_reprice import mixed_norm
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_signed_coordinate import coordinate_errors
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from full_gamma_zero_channel_jet_refine import scaled_forward
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, R, SIZE = 640, 32768, 56, 18
PINS = dict(FULL_ROW_PINS)
PINS.update({
    'full-gamma-joint-far-row.json.zlib':
        '00b9e4bdd8a99da4d6fdafc696d313bd43e4060e1ba7d8cec316f7cc5bcb7707',
    'full-gamma-residual-envelope.json.zlib':
        '7ca0ccd8a23dee5a0a60b793d21c07d287652d017afa387e28eb3a6243122d94',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
    'far-coupling.json.zlib':
        'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
})


def midpoint(pair, bits=BITS):
    """Exact midpoint of an ordered saved integer-grid interval."""
    if (type(bits) is not int or bits < 1 or len(pair) != 2
            or int(pair[0]) > int(pair[1])):
        raise ValueError('Require an ordered integer-grid interval.')
    return fmpq(int(pair[0])+int(pair[1]), 1 << (bits+1))


def fixed_trial(row, raw_core, reference, gamma, bits=BITS):
    """One exact reference action; return c=Q0 e*, not a computed inverse."""
    n = len(row)
    if (n < 3 or len(raw_core) != n or any(len(r) != n for r in raw_core)
            or (reference.nrows(), reference.ncols()) != (n-1, n-1)
            or not isinstance(reference, fmpq_mat) or reference != reference.transpose()
            or not isinstance(gamma, fmpq) or not 0 < gamma <= 1):
        raise ValueError('Invalid exact midpoint trial dimensions or reference.')
    e = fmpq_mat([[midpoint(p, bits) for p in row]])
    raw = fmpq_mat([[midpoint(p, bits) for p in r] for r in raw_core])
    core0 = (raw+raw.transpose())/2
    delta0 = -core0[0, 0]
    if not delta0 > 0:
        raise ValueError('No negative midpoint pivot.')
    t0 = fmpq_mat([[-core0[0, j+1]/core0[0, 0] for j in range(n-1)]])
    T0 = fmpq_mat([[t0[0, j] if i == 0 else fmpq(int(i == j+1))
                   for j in range(n-1)] for i in range(n)])
    z = e*T0
    action = gamma*reference*z.transpose()
    covector = T0*action
    covector[0, 0] -= e[0, 0]/delta0
    qinv0 = (e*covector)[0, 0]
    return dict(e=e, c=covector, t0=t0, delta0=delta0,
                k0=-e[0, 0]/delta0, qinv0=qinv0)


def centered_row_action(row, center, covector):
    """Keep the exact row-midpoint action separate from its interval error."""
    n = covector.nrows()
    if ((row.nrows(), row.ncols(), center.nrows(), center.ncols(), covector.ncols())
            != (1, n, 1, n, 1) or not isinstance(center, fmpq_mat)
            or not isinstance(covector, fmpq_mat)):
        raise ValueError('Incompatible exact anchored row action.')
    require_finite(list(row.entries()))
    return arb((center*covector)[0, 0])+((row-arb_mat(center))*arb_mat(covector))[0, 0]


def dual_profile(trial, sigma, rho, beta, rho0, beta0, alpha0, b, g, a):
    """Actual-Q dual norm with anchored core variation; no inverse/solve."""
    c, t0 = trial['c'], trial['t0']
    n = c.nrows()
    if ((sigma.nrows(), sigma.ncols(), t0.nrows(), t0.ncols()) != (n, n, 1, n-1)
            or not trial['delta0'] > 0 or not b > 0 or not a > 0):
        raise ValueError('Missing common-core and plane domains.')
    require_finite(list(sigma.entries())+[arb(x) for x in (alpha0, b, g, a)])
    delta = -sigma[0, 0]
    if not delta > 0:
        raise ValueError('No actual negative pivot.')
    t = arb_mat([[-sigma[0, j+1]/sigma[0, 0] for j in range(n-1)]])
    ctail = fmpq_mat([[c[j+1, 0]] for j in range(n-1)])
    dk = -((t-arb_mat(t0))*arb_mat(ctail))[0, 0]
    k0 = arb(trial['k0'])
    qinv = arb(trial['qinv0'])-(delta-arb(trial['delta0']))*k0**2 \
        -delta*(2*k0*dk+dk**2)
    rc = centered_row_action(rho, rho0, c)
    bc = centered_row_action(beta, beta0, c)
    dy = arb(alpha0)*rc+bc
    dv = rc+(g/b)*dy
    denominator = qinv+dy**2/b-dv**2/a
    result = dict(core_dual_norm=qinv, timelike_covector_action=dy,
        positive_covector_action=dv, projected_dual_norm=denominator,
        first_row_anchor_shift=dk)
    require_finite(list(result.values()))
    return result


def dual_residual_lower(cross, denominator_upper):
    """PSD residual Cauchy; a positive denominator LOWER is unnecessary."""
    require_finite([cross, arb(denominator_upper)])
    if not isinstance(denominator_upper, fmpq) or denominator_upper < 0:
        raise ValueError('Require a nonnegative exact dual-norm upper.')
    distance = max(fmpq(0), freeze_lower(abs(cross).lower()))
    if denominator_upper == 0:
        if distance > 0:
            raise ArithmeticError('Zero dual norm conflicts with a nonzero residual cross.')
        return fmpq(0)
    return distance**2/denominator_upper


def consume(out):
    started, source = time.perf_counter(), Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p:hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different fixed common-trial input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    inherited = dict(PINS)
    for packet in packets.values():
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256',
                      'cached_producer_source_sha256'):
            for name, digest in packet.get(field, {}).items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting common-trial binding: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        p = HERE/name
        if p in guards and guards[p] != digest:
            raise ValueError('Different invoked common-trial dependency: '+name)
        guards[p] = digest
    verify_guards(guards)
    get = lambda name: packets[name]
    row = get('full-gamma-joint-far-row.json.zlib')
    saved = get('full-gamma-residual-envelope.json.zlib')
    cached = get('full-gamma-cached-source-response.json.zlib')
    cauchy = get('full-gamma-signed-mu-direct-sine.json.cauchy.zlib')
    mu = get('full-gamma-signed-mu-direct-sine.json.zlib')
    seed = get('tail-moments.json.zlib')
    lift = get('full-gamma-lift-projection.json.zlib')
    signed = get('full-gamma-signed-correction.json.zlib')
    old = get('full-gamma-negative-channel-capacity.json.zlib')
    hankel = get('full-gamma-hankel-tail-energy.json.zlib')
    core = get('joint-gamma-6765-rho32-v2-schur.json')
    cert = get('full-gamma-directional-metric.json.zlib')
    projection = get('full-gamma-projection-upper.json.zlib')
    actual = get('full-gamma-actual-channel-projection.json.zlib')
    plane = get('full-gamma-joint-plane-paid-floor.json.zlib')
    target = get('structured-joint-2048-target.json.zlib')
    far = get('far-coupling.json.zlib')
    primal = get('full-gamma-primal-tail.json')
    rational = lambda p,key: fmpq(*map(int,p['exact_rationals'][key]))
    ball = lambda p,key: interval_ball(p['scalar_bounds'][key],BITS)
    upper = lambda p,key: fmpq(int(p['scalar_bounds'][key][1]),1 << BITS)
    expected = {name:row['packet_sha256'][name] for name in FULL_ROW_PINS}
    if (expected != FULL_ROW_PINS or any(p['grid_bits'] != BITS for p in packets.values())
            or row['schema'] != 'ACTUAL_C27_JOINT_FAR_ROW_V1'
            or saved['schema'] != 'ACTUAL_C27_RESIDUAL_ENVELOPE_V1'
            or any((p['source_bits'],p['source_row'],p['source_integer_sha256']) !=
                   (384,0,row['source_integer_sha256'])
                   for p in (saved,cached,cauchy,mu,signed,old,hankel,plane))
            or saved['packet_sha256']['full-gamma-joint-far-row.json.zlib'] !=
                PINS['full-gamma-joint-far-row.json.zlib']
            or any(saved['packet_sha256'][n] != PINS[n] for n in saved['packet_sha256'])
            or all(v['residual_only_refinement_cannot_certify_positive']
                   for v in saved['recipes'].values())
            or projection['exact_source_integer_sha256'][0] != row['source_integer_sha256']
            or projection['exact_source_bits'][0] != 384
            or any(p['remaining_modes'] != lift['remaining_modes']
                   or p['input_sha256'] != lift['input_sha256']
                   for p in (core,cert,projection,actual,old))
            or len(lift['remaining_modes']) != SIZE
            or (lift['N'],lift['M'],lift['last_mode']) != (1280,M,2048)
            or lift['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or cached['recovered_symbols_are_actual'] or not cached['inherited_symbol_error_paid']
            or seed['analytic_grid_bits'] != mu['symbol_metadata_grid_bits']
            or (seed['c'],seed['M']) != (27,M)
            or (cached['finite_first_mode'],cached['finite_last_mode'],
                cached['analytic_first_mode'],cached['normalization_scale']) != (M+1,2*M,2*M+1,M)
            or primal['metadata']['target_sha256'] != PINS['structured-joint-2048-target.json.zlib']
            or primal['metadata']['input_sha256'] != lift['input_sha256']
            or target['metadata']['input_sha256'] != lift['input_sha256']
            or (far['c'],far['N'],far['M'],far['jet_order'],far['moment_diagonals'],
                far['bound_grid_bits']) != (27,1280,M,R,R,628)):
        raise ValueError('Different SAME-source/form, cached error or jet contracts.')
    shapes = ((row['actual_far_row'],1,SIZE), (lift['far_features_on_lift'],112,SIZE),
        (lift['lift_gram'],SIZE,SIZE), (core['reduced_schur_intervals'],SIZE,SIZE),
        (cert['reference_inverse'],SIZE-1,SIZE-1), (target['kernel'],112,112),
        (cached['physical_left_kernel_action'],1,112),
        ([projection['matrices']['source_on_lift'][0]],1,SIZE),
        (actual['actual_finite_row_on_lift'],1,SIZE))
    width_ok = lambda pair: len(pair) == 2 and 0 <= int(pair[1])-int(pair[0]) <= 2
    if (any(len(rows) != n or any(len(r) != m for r in rows) for rows,n,m in shapes)
            or len(cached['recovered_approximate_symbols']) != M
            or len(cauchy['finite_cauchy_source_outputs']) != M
            or not 1 <= len(cached['interpolated_source_coefficients']) <= 13
            or len(seed['normalized_polynomial']) != 64
            or any(str(q) not in mu['real_hurwitz_kernels'] for q in range(2,188))
            or len(far['moment_pairs']) != R or len(far['zeta']) != R
            or {r['U_power'] for r in far['moment_pairs']} != set(range(2,2*R+1,2))
            or any(r['T_power'] != r['U_power']+1 for r in far['moment_pairs'])
            or {r['p'] for r in far['zeta']} != set(range(4,2*R+3,2))
            or any(not width_ok(r[k]) for r in far['moment_pairs'] for k in ('U','T'))
            or any(not width_ok(r['interval']) for r in far['zeta'])
            or not width_ok(far['alpha_pair'])):
        raise ValueError('Incomplete common-trial dimensions or retained moments.')
    kappa, weight = rational(cached,'kappa'), rational(cached,'original_far_weight')
    gamma, alpha0 = rational(plane,'gamma'), rational(plane,'alpha0')
    eps_far = rational(signed,'original_far_matrix_error')
    kphys = fmpq(int(target['metadata']['trial_operator_norm_sq_upper'][1]),1 << 512)
    eta = upper(old,'neumann_eta_upper')
    if (kappa != fmpq(33,32) or kappa != rational(signed,'kappa')
            or kappa != rational(mu,'kappa') or not weight > 0
            or weight != rational(signed,'original_far_weight')
            or weight != rational(mu,'original_far_weight')
            or weight != fmpq(int(far['D_numerator']),1 << far['D_grid_bits'])
            or not 0 <= eta < 1 or gamma != 1/(1+eta)
            or gamma != rational(saved,'gamma') or alpha0 != rational(old,'fixed_alpha')
            or eps_far != fmpq(int(target['metadata']['far_matrix_error'][1]),1 << 512)
            or kphys != rational(signed,'physical_operator_norm_sq_upper') or not kphys > 0
            or rational(signed,'moment_gram_norm_upper') !=
                fmpq(int(far['upper_bounds']['moment_gram_norm']),1 << 628)):
        raise ValueError('Different original weight, comparison or uniform error units.')
    reference = fmpq_mat([[fmpq(v) for v in r] for r in cert['reference_inverse']])
    trial = fixed_trial(row['actual_far_row'][0],core['reduced_schur_intervals'],reference,gamma)
    sigma = matrix_in(core['reduced_schur_intervals'],BITS)
    sigma = (sigma+sigma.transpose())/2
    rrows, brows = [projection['matrices']['source_on_lift'][0]], actual['actual_finite_row_on_lift']
    rho, beta = matrix_in(rrows,BITS), matrix_in(brows,BITS)
    exact_row = lambda rows: fmpq_mat([[midpoint(p) for p in r] for r in rows])
    b,g,a = ball(plane,'timelike_energy'),ball(plane,'rho_timelike_pairing'),ball(cached,'fixed_positive_plane_energy')
    scalars = dual_profile(trial,sigma,rho,beta,exact_row(rrows),exact_row(brows),alpha0,b,g,a)
    dupper = freeze_upper(scalars['projected_dual_norm'].upper())
    if dupper < 0:
        raise ArithmeticError('Dual-norm upper contradicts its structural nonnegativity.')
    exact = dict(projected_dual_norm_upper=dupper, midpoint_dual_norm=trial['qinv0'],
                 midpoint_pivot_magnitude=trial['delta0'], midpoint_first_anchor=trial['k0'])
    counters = dict(new_aggregated_feature_columns=0,new_masked_physical_quadratics=0,
        new_directional_norm_applications=0,new_finite_functionals=0,
        new_cached_tail_applications=0,new_mixed_dots=0)
    cross = arb(0)
    if dupper > 0:
        c = trial['c']
        hc = matrix_in(lift['far_features_on_lift'],BITS)*arb_mat(c)
        features = [hc[i,0] for i in range(112)]
        lcap = sum((abs(c[j,0])*freeze_upper(interval_ball(lift['lift_gram'][j][j],BITS).upper().sqrt())
                    for j in range(SIZE)),fmpq(0))
        old_lift, old_forward = rational(signed,'lift_norm_upper'),rational(signed,'sharp_forward_error_upper')
        forward = scaled_forward(old_forward,old_lift,lcap)
        raw_mixed_error = scaled_forward(rational(signed,'mixed_cross_error_upper'),old_lift,lcap)
        coefficient_norm = freeze_upper(arb(kphys).sqrt()*arb(lcap))
        quadratic,gram_error,broad2,_ = mixed_norm(hc,matrix_in(target['kernel'],BITS),eps_far,coefficient_norm)
        anchor_row = next(r for r in far['moment_pairs'] if r['U_power'] == 2*R)
        anchor_zeta = next(r for r in far['zeta'] if r['p'] == 2*R+2)
        anchor = [interval_ball(p,BITS) for p in (anchor_row['U'],anchor_row['T'],anchor_zeta['interval'])]
        norm_detail = directional_norm(features,freeze_upper(quadratic.upper()),
            freeze_upper(interval_ball(primal['bounds']['actual_symbol_budget'],BITS).upper()),
            fmpq(int(far['alpha_pair'][1]),1 << BITS),
            fmpq(int(far['alpha_pair'][1])-int(far['alpha_pair'][0]),1 << (BITS+1)),
            weight,fmpq(2*R,1 << BITS),rational(signed,'moment_gram_norm_upper'),
            forward,R,anchor,fmpq(M,M+1))
        norm2 = min(broad2,norm_detail['rebuilt_representer_energy_upper'])
        ncap = freeze_upper(arb(norm2).sqrt())
        errors = coordinate_errors(features[:R],
            interval_ball(seed['uniform_symbol_error'],seed['analytic_grid_bits']).upper(),weight,
            rational(cached,'source_tail_norm_error_upper'),ncap,forward,
            upper(hankel,'principal_weighted_norm_upper'),upper(hankel,'smooth_remainder_weighted_norm_upper'),
            raw_mixed_error,kappa)
        samples = [interval_ball(p,BITS) for p in cauchy['finite_cauchy_source_outputs']]
        symbols = [interval_ball(p,BITS) for p in cached['recovered_approximate_symbols']]
        finite = finite_pairing(samples,symbols,features,M,kappa,weight)
        polynomial = arb_poly([interval_ball(p,BITS) for p in cached['interpolated_source_coefficients']])
        smooth = arb_poly([interval_ball(p,BITS) for p in seed['normalized_polynomial']])
        kernels = {q:interval_ball(p,BITS) for q,p in mu['real_hurwitz_kernels'].items()}
        tail,tail_cap,oscillatory = cached_tail(polynomial,features,smooth,kernels,
            upper(mu,'inverse_phase_gap_weight_cap'),M,kappa,weight)
        mixed = (arb(kappa)*matrix_in(cached['physical_left_kernel_action'],BITS)*hc)[0,0]
        allowance = errors['nonoscillatory_error_upper']+oscillatory
        fc = finite+tail+mixed+arb(0,arb(allowance))
        cross = fc+ball(cached,'fixed_transported_theta')*scalars['timelike_covector_action']/b \
            -ball(saved,'positive_plane_response')*scalars['positive_covector_action']/a
        scalars.update(aggregate_finite=finite,aggregate_tail=tail,aggregate_mixed=mixed,
            aggregate_tail_coefficient_cap=tail_cap,aggregate_actual_pairing=fc,
            aggregate_masked_quadratic=quadratic)
        exact.update(norm_detail,**errors,lift_norm_upper=lcap,sharp_forward_error_upper=forward,
            physical_coefficient_norm_upper=coefficient_norm,physical_gram_error_upper=gram_error,
            broad_physical_norm2_upper=broad2,selected_physical_norm2_upper=norm2,
            selected_physical_norm_upper=ncap,oscillatory_error_upper=oscillatory,
            complete_error_upper=allowance)
        counters = {key:1 for key in counters}
    new_lower = dual_residual_lower(cross,dupper)
    previous = rational(saved,'selected_residual_lower')
    retained = max(previous,new_lower)
    if not 0 <= previous <= retained <= rational(saved,'selected_residual_upper'):
        raise ArithmeticError('Common-trial residual conflicts with saved SAME residual envelope.')
    exact.update(new_residual_lower=new_lower,previous_residual_lower=previous,
        retained_residual_lower=retained,residual_upper=rational(saved,'selected_residual_upper'))
    scalars['projected_residual_cross'] = cross
    recipes = {}
    for recipe in ('fixed','energy_endpoint'):
        same_F = rational(cached,recipe+'_all_alpha_floor_lower')
        prior = rational(row,recipe+'_retained_lower')
        if same_F != rational(saved,recipe+'_same_F_lower') or prior != rational(saved,recipe+'_previous_retained_lower'):
            raise ValueError('Different SAME-F lower certificate.')
        augmented = same_F+retained
        kept = max(prior,augmented)
        exact.update({recipe+'_same_F_lower':same_F,recipe+'_previous_retained_lower':prior,
                      recipe+'_residual_augmented_lower':augmented,recipe+'_retained_lower':kept})
        recipes[recipe] = dict(strictly_improved=bool(kept > prior),positive_paid_lower_certified=bool(kept > 0))
    require_finite(list(scalars.values()))
    exported = {k:envelope(v,BITS) for k,v in scalars.items()}
    if any(p is None or int(p[0]) > int(p[1]) for p in exported.values()):
        raise ArithmeticError('Invalid common-trial export; no retry.')
    exact_out = lambda x: [str(x.p),str(x.q)]
    result = dict(schema='ACTUAL_C27_COMMON_TRIAL_RESIDUAL_V1',packet_sha256=PINS,
        source_sha256=guards[source],invocation_source_sha256={p.name:h for p,h in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name:h for p,h in guards.items() if p not in invoked},
        source_integer_sha256=row['source_integer_sha256'],source_bits=384,source_row=0,
        grid_bits=BITS,working_bits=ctx.prec,remaining_modes=lift['remaining_modes'],
        exact_trial_row=[[exact_out(x) for x in trial['e'].entries()]],
        exact_trial_covector=[[exact_out(x)] for x in trial['c'].entries()],
        scalar_bounds=exported,exact_rationals={k:exact_out(v) for k,v in exact.items()},recipes=recipes,
        **counters,new_exact_reference_actions=1,new_core_variation_profiles=1,
        original_coordinate_queries=0,new_inverse_or_solve=0,
        old_scientific_producers_or_fixtures_replayed=0,source_direction_parameter_precision_cutoff_scan=False,
        strictly_improved_residual_lower=bool(retained > previous),
        actual_paid_energy_upper_claim=False,actual_center_sign_claim=False,
        actual_full_border_sign_claim=False,actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'Reuse SAME real fixed18-row f, original Q=gamma*T*Rref*T*-e0e0*/delta with exact SPD Rref, y0 timelike Q(y0,y0)=-b and v=rho+(g/b)y0 orthogonal to y0 with Q(v,v)=a>0. All inherited common-source and physical-error premises remain.',
            'Unique e is the exact dyadic midpoint of the saved full row. Exact symmetric core midpoint defines T0,delta0,Q0; one reference action gives c=Q0*e*. No inverse or solve is computed, and no alternative direction is queried.',
            'Actual-Q dual norm c*Q^-1*c=gamma*z*Rref*z*-delta*(c0-t*ctail)^2 uses ctail=gamma*Rref*z*, z=e*T0. Anchor at qinv0=e*c and k0=-e0/delta0 pays actual core variation explicitly. Finite rho,beta actions use exact midpoint anchors and retain all original row interval errors.',
            'The abstract Riesz representer Q^-1*c is projected off the SAME y0,v plane. Its squared norm D=c*Q^-1*c+(y0*c)^2/b-(v*c)^2/a is structurally nonnegative. Cross L=f*c+theta*(y0*c)/b-s*(v*c)/a uses actual s=gamma*zeta+(coordinate_ratio+g/b)*theta from the accepted SAME-source envelope.',
            'PSD residual Cauchy gives Rperp>=dist(0,I_L)^2/U_D for U_D>0. A positive numerical lower for D is unnecessary; Dupper<0 is inconsistent and Dupper=0 cannot supply a gain. Independent interval endpoints need not be jointly attainable.',
            'Aggregate the saved 112x18 feature map on this exact c before bounding errors. Original weighted physical kernel needs no extra weight or kappa. Pay broad Gram or complete directional-jet masking/moments/rounding/actual-symbol sharp forward, plus source/symbol/smooth/raw-mixed/Abel errors. Raw mixed allowance already includes kappa.',
            'Uniform forward and raw mixed errors scale by the SAME safe lift-norm upper. The full jet moments depend on the infinite high support, not whether physical low support ends at1280 or2048. Approximate recovered symbols remain approximate and their original error is paid.',
            'Retain max(old,new) residual lower, never sum. Add retained residual only to each saved SAME-F lower, then max with its previous retained paid certificate. Nonpositive lower bounds give no actual sign conclusion.',
            'Positive paid lower, if any, concerns only the two fixed-source paid comparison recipes under inherited premises. No actual-center/full-border/cofinal positivity, Robin, RH, generic novelty or Lean proof.',
        ])
    verify_guards(guards)
    digest = write_packet(out,result)
    print('saved common trial residual',digest,'seconds',result['elapsed_seconds'])
    for key in ('projected_dual_norm_upper','new_residual_lower','retained_residual_lower',
                'fixed_retained_lower','energy_endpoint_retained_lower'):
        print(key,arb(exact[key]).str(30))
    print('projected cross',cross.str(30))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; fixed1024/640, one trial, no retry.')
    ctx.prec = 1024
    consume(args.out)
