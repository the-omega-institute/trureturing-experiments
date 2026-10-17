"""Price the unknown far inverse energy using accepted same-source packets.

Only new saved-data contractions run: one exact17-reference action on the
complete lifted row, two bilinears, and the two-dimensional projection floor.
No source, far, target, old solve, inverse, factorization or scan is replayed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd', 'full-gamma-directional-metric.json.zlib': 'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2', 'full-gamma-negative-channel-capacity.json.zlib': '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e', 'full-gamma-projection-upper.json.zlib': 'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3', 'full-gamma-actual-channel-projection.json.zlib': '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745', 'full-gamma-surplus-budget.json': '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f', 'full-gamma-band-surplus-budget.json.zlib': '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1', 'full-gamma-joint-far-row.json.zlib': '00b9e4bdd8a99da4d6fdafc696d313bd43e4060e1ba7d8cec316f7cc5bcb7707', 'full-gamma-finite-good-block-metric.json.zlib': '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072', 'full-gamma-good-block-far-cross.json.zlib': '39ddee410fa83bf4d35cb97f33bc9e63568dff8d9404dee331b0b3925336a858', 'full-gamma-lift-projection.json.zlib': 'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873'}
CORE, CERT, CAPACITY, PROJECTION, FINITE, SURPLUS, BAND, FAR_ROW, METRIC, CROSS, LIFT = PINS


def finite(values):
    if any(not value.is_finite() for value in values):
        raise ArithmeticError('Nonfinite saved-data budget input or result.')


def projection_floor(gram, mixed):
    """Bessel floor m*G^-1*m, using a stable explicit real2x2 Schur form."""
    if (gram.nrows(), gram.ncols(), mixed.nrows(), mixed.ncols()) != (2, 2, 1, 2):
        raise ValueError('Require a real2x2 joint metric and two common mixed values.')
    finite(list(gram.entries())+list(mixed.entries()))
    a, c, d = gram[0, 0], gram[0, 1], gram[1, 1]
    if not a > 0:
        raise ArithmeticError('Missing positive first metric pivot.')
    schur = d-c*c/a
    if not schur > 0:
        raise ArithmeticError('Missing positive joint metric Schur pivot.')
    first = mixed[0, 0]**2/a
    remainder = (mixed[0, 1]-c*mixed[0, 0]/a)**2/schur
    floor = first+remainder
    finite([schur, floor])
    return dict(metric_schur=schur, projection_first=first,
                projection_second=remainder, far_energy_projection_floor=floor)


def core_terms(sigma, reference, rho, beta, saved_rr, eta):
    """Reuse a certified reference inverse; a single NEW complete-row action.

For A=Pplus and R=reference, the accepted relative eta implies
(1+eta)^-1 R <= A^-1 <= (1-eta)^-1 R. Therefore every bilinear
perturbation is bounded by eta/(1-eta) times its reference Cauchy cap.
The negative pivot term uses the actual enclosed first row of Sigma0.
"""
    n = sigma.nrows()
    if (n < 2 or sigma.ncols() != n
            or (reference.nrows(), reference.ncols()) != (n-1, n-1)
            or (rho.nrows(), rho.ncols(), beta.nrows(), beta.ncols()) != (1, n, 1, n)):
        raise ValueError('Incompatible actual core, reference and complete source rows.')
    finite(list(sigma.entries())+list(reference.entries())+list(rho.entries())
           +list(beta.entries())+[saved_rr, arb(eta)])
    pivot = sigma[0, 0]
    if not pivot < 0 or not 0 <= eta < 1 or not saved_rr.upper() >= 0:
        raise ArithmeticError('Missing negative pivot or certified relative reference bound.')
    t = arb_mat([[-sigma[0, j+1]/pivot if i == 0 else arb(int(i == j+1))
                  for j in range(n-1)] for i in range(n)])
    rt, bt = rho*t, beta*t
    action = reference*bt.transpose()
    rb, bb = (rt*action)[0, 0], (bt*action)[0, 0]
    if not bb.upper() >= 0:
        raise ArithmeticError('Complete positive-row form conflicts with reference positivity.')
    relative = arb(eta)/(1-arb(eta))
    rr_cap, bb_cap = max(arb(0), saved_rr.upper()), max(arb(0), bb.upper())
    er = (relative*rr_cap).upper()
    eb = (relative*bb_cap).upper()
    ec = (relative*(rr_cap*bb_cap).sqrt()).upper()
    rr_core = rho[0, 0]**2/pivot+saved_rr+arb(0, er)
    rb_core = rho[0, 0]*beta[0, 0]/pivot+rb+arb(0, ec)
    bb_core = beta[0, 0]**2/pivot+bb+arb(0, eb)
    result = dict(reference_complete_cross=rb, reference_complete_energy=bb,
        core_rr=rr_core, core_rb=rb_core, core_bb=bb_core,
        reference_rr_error=arb(er), reference_rb_error=arb(ec), reference_bb_error=arb(eb))
    finite(list(result.values()))
    return result


def budget(gram, mixed, core, t_plus, charge):
    if not t_plus > 0 or not charge > 0:
        raise ArithmeticError('Require positive comparison coefficient and fixed valid charge.')
    result = projection_floor(gram, mixed)
    a = 1/arb(t_plus)+gram[0, 0]+core['core_rr']
    off = gram[0, 1]+mixed[0, 0]+core['core_rb']
    d0 = arb(charge)+gram[1, 1]+2*mixed[0, 1]+core['core_bb']
    base = a*d0-off**2
    result.update(border_a=a, border_cross=off, border_d_without_far_energy=d0,
        determinant_constant=base, determinant_far_energy_coefficient=a,
        determinant_at_projection_floor=base+a*result['far_energy_projection_floor'])
    if a > 0:
        critical = off**2/a-d0
        delta = critical-result['far_energy_projection_floor']
        result.update(far_energy_critical=critical, remaining_energy_budget=delta)
        verdict = ('fixed_comparison_excluded_by_projection_floor' if delta.upper() <= 0
                   else 'positive_remaining_budget_only' if delta.lower() > 0
                   else 'inconclusive_remaining_budget_box')
    elif a < 0:
        verdict = 'negative_first_entry_of_optimistic_comparison_only'
    else:
        verdict = 'inconclusive_first_entry_zero_box'
    finite(list(result.values()))
    return result, verdict


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different same-source budget input: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    core,cert,capacity,projection,actual,surplus,band,far,metric,cross,lift = packets.values()
    modes = metric['remaining_modes']
    if (len(modes) != 18 or len(metric['positive_modes']) != 1263
            or sorted(modes+metric['positive_modes']) != list(range(1281))
            or any(p['remaining_modes'] != modes for p in (core,cert,capacity,projection,actual,far,lift))
            or any(p['input_sha256'] != core['input_sha256'] for p in (cert,capacity,projection,actual,lift))
            or any(p['grid_bits'] != BITS for p in (core,cert,capacity,projection,actual,band,far,metric,cross,lift))
            or cross['packet_sha256'][METRIC] != PINS[METRIC]
            or metric['packet_sha256']['joint-gamma-6765-rho32-v2-target.json.zlib'] != core['target_packet_sha256']
            or lift['base_core_sha256'] != PINS[CORE]
            or projection['packet_sha256']['core'] != PINS[CORE]
            or actual['packet_sha256']['projection'] != PINS[LIFT]
            or cert['packet_sha256']['core'] != PINS[CORE]
            or capacity['packet_sha256'][CERT] != PINS[CERT]
            or far['packet_sha256'][CERT] != PINS[CERT]
            or metric['joint_row_order'] != ['complete_unweighted_r','finite_bJ']
            or cert['verdict'] != 'directional-center-enclosed'
            or not cert['exact_reference_identities_checked'] or cert['positive_dimension'] != 17
            or cross['new_profiles'] != ['P^-1*r_g*','P^-1*bJ_g*']
            or cross['complete_joint_metric_supplied'] or metric['complete_joint_metric_supplied']
            or any((p['source_bits'],p['source_row'],p['source_integer_sha256']) != (384,0,SOURCE_ID)
                   for p in (capacity,surplus,band,far,metric,cross))
            or projection['exact_source_bits'][0] != 384
            or projection['independent_source_indices'][0] != 0
            or projection['exact_source_integer_sha256'][0] != SOURCE_ID
            or (core['N'],core['M'],core['last_mode']) != (1280,32768,2048)
            or metric['basis'] != core['basis'] or cross['basis'] != core['basis']
            or metric['positive_modes'] != cross['positive_modes']):
        raise ValueError('Different actual source, old core, same joint metric or reference roles.')
    raw = core['reduced_schur_intervals']
    if (len(raw) != 18 or any(len(row) != 18 for row in raw)
            or len(cert['reference_inverse']) != 17 or any(len(row) != 17 for row in cert['reference_inverse'])
            or len(actual['actual_finite_row_on_lift']) != 1 or len(actual['actual_finite_row_on_lift'][0]) != 18
            or len(far['actual_far_row']) != 1 or len(far['actual_far_row'][0]) != 18
            or len(projection['matrices']['source_on_lift'][0]) != 18
            or metric['joint_inverse_metric'][0][1] != metric['joint_inverse_metric'][1][0]):
        raise ValueError('Incomplete saved core, source rows or symmetric metric.')
    rational = lambda p,k: fmpq(*map(int,p['exact_rationals'][k]))
    cap = rational(band,'supplied_parameter_cap')
    h_lower, h_upper = [rational(surplus,k) for k in ('h_lower','h_upper')]
    if not 0 < cap < 1 or not 0 < h_lower <= h_upper:
        raise ArithmeticError('Missing positive supplied cap and h endpoints.')
    t_plus = cap/h_lower
    charge = fmpq(int(actual['scalar_bounds']['refined_common_charge_upper'][1]),1 << BITS)
    eta = fmpq(int(cert['scalar_bounds']['neumann_eta_upper'][1]),1 << BITS)
    sigma = symmetric(matrix_in(raw,BITS))
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    rho = matrix_in([projection['matrices']['source_on_lift'][0]],BITS)
    beta = matrix_in(actual['actual_finite_row_on_lift'],BITS)+matrix_in(far['actual_far_row'],BITS)
    saved_rr = interval_ball(capacity['scalar_bounds']['steering_reference_rRr'],BITS)
    terms = core_terms(sigma,reference,rho,beta,saved_rr,eta)
    gram = matrix_in(metric['joint_inverse_metric'],BITS)
    mixed = matrix_in(cross['actual_far_cross'],BITS)
    scalars,verdict = budget(gram,mixed,terms,t_plus,charge)
    result = dict(schema='ACTUAL_C27_FAR_INVERSE_ENERGY_BUDGET_V1',packet_sha256=PINS,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID,source_bits=384,source_row=0,
        N=1280,M=32768,remaining_modes=modes,grid_bits=BITS,working_bits=ctx.prec,
        complete_lifted_row=matrix_out(beta,BITS),
        scalar_bounds={k:envelope(v,BITS) for k,v in dict(**terms,**scalars).items()},
        exact_rationals={k:[str(v.p),str(v.q)] for k,v in dict(comparison_t_plus=t_plus,
            fixed_charge=charge,neumann_eta_upper=eta).items()},verdict=verdict,
        comparison_coefficient_eligibility_claim=False,new_reference_actions=1,
        new_reference_bilinears=2,reused_reference_rRr=True,new_two_dimensional_projection_form=1,
        new_source_far_target_inverse_solve_queries=0,old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,complete_joint_metric_supplied=False,
        actual_center_sign_claim=False,actual_full_border_sign_claim=False,actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'All rows and inverse metrics refer to ONE original real384-bit row0 source, the SAME1263-positive principal block and old18-core. The full old18far row is consumed without another source query.',
            'The accepted exact17reference inverse and uniform eta<1 enclose positive-block inverse bilinears. The saved rRr is reused. The sole new reference action is on the complete betaJ+far row.',
            'The true far inverse energy z is at least m*G_JJ^-1m. This projection lower is not the missing upper. Shared containing intervals do not assert independent endpoint attainability.',
            't_plus=cap/h_lower is an optimistic upper comparison, not asserted admitted. D is the fixed refined common charge upper. Any exclusion applies only to 0<t<=t_plus and D>=that fixed charge on the SAME lower target.',
            'When a>0, detN2=a*(z-zcrit). Its strict negative sign requires z<zcrit. A positive remaining budget is not a sign certificate; a box crossing zero is inconclusive without retries.',
            'All inherited common form/domain, residual-map, source/lift, complete-far, normalization and comparison premises remain. Classical projection/Schur/Neumann mechanisms only; no generic novelty or Lean proof. Actual-center/full-border/cofinal positivity, Robin and RH remain open.',
        ])
    if any(v is None for v in result['scalar_bounds'].values()):
        raise ArithmeticError('Nonfinite budget export; no retry.')
    print('saved far inverse-energy budget',write_packet(out,result),'seconds',result['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; no overwrite or retry.')
    ctx.prec = 1024
    consume(args.out)
