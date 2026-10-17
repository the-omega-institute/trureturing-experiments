"""Apply the accepted source functional to two NEW good-block profiles.

Reuse the saved1263x2 inverse actions. Aggregate the original112-dimensional
far map before source/physical/tail error bounds. No old lift profile,
source recovery, target reconstruction, solve or general inverse runs.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, arb_poly, ctx, fmpq
from diagonal_floor_witness import interval_ball
from full_gamma_cached_source_response import finite_pairing, cached_tail
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_common_trial_residual import PINS as ACCEPTED_PINS
from full_gamma_directional_jet_norm import directional_norm
from full_gamma_mixed_norm_reprice import mixed_norm
from full_gamma_signed_coordinate import coordinate_errors
from full_gamma_signed_mu import require_finite
from full_gamma_unweighted_retained_row import SOURCE_ID
from full_gamma_zero_channel_jet_refine import scaled_forward
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, R = 640, 32768, 56
METRIC = 'full-gamma-finite-good-block-metric.json.zlib'
ROW = 'full-gamma-unweighted-retained-row.json.zlib'
TARGET = 'structured-joint-2048-target.json.zlib'
CACHED = 'full-gamma-cached-source-response.json.zlib'
CAUCHY = 'full-gamma-signed-mu-direct-sine.json.cauchy.zlib'
MU = 'full-gamma-signed-mu-direct-sine.json.zlib'
SEED = 'tail-moments.json.zlib'
SIGNED = 'full-gamma-signed-correction.json.zlib'
HANKEL = 'full-gamma-hankel-tail-energy.json.zlib'
FAR = 'far-coupling.json.zlib'
PRIMAL = 'full-gamma-primal-tail.json'
LIFT = 'full-gamma-lift-projection.json.zlib'
PINS = {name: ACCEPTED_PINS[name] for name in
        (TARGET,CACHED,CAUCHY,MU,SEED,SIGNED,HANKEL,FAR,PRIMAL,LIFT)}
PINS.update({METRIC: '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072',
             ROW: '566cb68530f31dd6e26af076055bb368e2b2e82fdcc0857abb7b2b475a036f48'})


def aggregate_good_profiles(features, positive, actions):
    """Embed the SAME two good-block actions; bound their full retained norms."""
    if (actions.ncols() != 2 or actions.nrows() != len(positive)
            or not positive or len(set(positive)) != len(positive)
            or any(type(n) is not int or not 0 <= n < features.nrows() for n in positive)
            or features.ncols() < 1):
        raise ValueError('Require two complete actions on distinct retained coordinates.')
    require_finite(list(features.entries())+list(actions.entries()))
    selected = arb_mat([[features[n,j] for j in range(features.ncols())] for n in positive])
    projected = selected.transpose()*actions
    caps = [freeze_upper(sum((abs(actions[i,j]).upper()**2
                for i in range(actions.nrows())),arb(0)).sqrt()) for j in range(2)]
    require_finite(list(projected.entries())+[arb(x) for x in caps])
    return projected, caps


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in PINS.items():
        packets[name], current = read_json(HERE/name)
        if current != digest:
            raise ValueError('Different accepted source/target/action input: '+name)
    target,cached,cauchy,mu,seed,signed,hankel,far,primal,lift = [
        packets[n] for n in (TARGET,CACHED,CAUCHY,MU,SEED,SIGNED,HANKEL,FAR,PRIMAL,LIFT)]
    metric,row = packets[METRIC],packets[ROW]
    if [p['schema'] for p in (target,cached,cauchy,mu,seed,signed,hankel,far,primal,lift,metric,row)] != [
            'ACTUAL_C27_STRUCTURED_TARGET_V1','ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1',
            'ACTUAL_C27_SIGNED_MU_CAUCHY_BLOCK_V1','ACTUAL_C27_SAME_SOURCE_SIGNED_MU_ABEL_V1',
            'ACTUAL_WEIL_NORMALIZED_TAIL_MOMENTS_V1','ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1','ACTUAL_C27_FAR_COUPLING_GRAM_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1','ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1','ACTUAL_C27_COMPLETE_UNWEIGHTED_RETAINED_ROW_V1']:
        raise ValueError('Different complete-source/profile/physical roles.')
    positive,remaining = metric['positive_modes'],metric['remaining_modes']
    if (any(p['grid_bits'] != BITS for p in packets.values())
            or any((p['source_bits'],p['source_row'],p['source_integer_sha256']) != (384,0,SOURCE_ID)
                for p in (cached,cauchy,mu,signed,hankel,metric,row))
            or len(positive) != 1263 or len(remaining) != 18
            or sorted(positive+remaining) != list(range(1281))
            or metric['joint_row_order'] != ['complete_unweighted_r','finite_bJ']
            or metric['packet_sha256'][ROW] != PINS[ROW]
            or metric['packet_sha256']['joint-gamma-6765-rho32-v2-target.json.zlib'] != row['base_target_sha256']
            or row['original_target_sha256'] != PINS[TARGET]
            or any(p['positive_modes'] != positive or p['remaining_modes'] != remaining for p in (row,lift))
            or lift['base_target_sha256'] != row['base_target_sha256']
            or lift['base_lift_sha256'] != row['base_lift_sha256']
            or lift['original_target_sha256'] != PINS[TARGET]
            or lift['input_sha256'] != row['input_sha256']
            or target['metadata']['input_sha256'] != row['input_sha256']
            or signed['packet_sha256'][TARGET] != PINS[TARGET]
            or (target['metadata']['N'],target['metadata']['M'],target['metadata']['last_mode']) != (1280,M,2048)
            or any(p['basis'] != target['metadata']['basis'] for p in (metric,row,lift))
            or row['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or (row['first_source_mode'],row['last_source_mode']) != (2049,M)
            or (cached['finite_first_mode'],cached['finite_last_mode'],cached['analytic_first_mode'],
                cached['normalization_scale']) != (M+1,2*M,2*M+1,M)
            or (cauchy['finite_first_mode'],cauchy['finite_last_mode'],cauchy['normalization_scale']) != (M+1,2*M,M)
            or cached['recovered_symbols_are_actual'] or not cached['inherited_symbol_error_paid']
            or (seed['c'],seed['M'],seed['analytic_grid_bits']) != (27,M,mu['symbol_metadata_grid_bits'])
            or (far['c'],far['N'],far['M'],far['jet_order'],far['moment_diagonals']) != (27,1280,M,R,R)
            or metric['complete_b_row_supplied'] or metric['complete_joint_metric_supplied']):
        raise ValueError('Different same-source, positive block, basis, weight or completion scope.')
    rational = lambda p,k: fmpq(*map(int,p['exact_rationals'][k]))
    upper = lambda p,k: fmpq(int(p['scalar_bounds'][k][1]),1 << BITS)
    kappa,weight = rational(cached,'kappa'),rational(cached,'original_far_weight')
    metadata_upper = lambda k: fmpq(int(target['metadata'][k][1]),1 << 512)
    kphys = rational(signed,'physical_operator_norm_sq_upper')
    eps = fmpq(lift['exact_far_allowance'])
    if (kappa != fmpq(33,32) or any(rational(p,'kappa') != kappa for p in (mu,signed))
            or any(rational(p,'original_far_weight') != weight for p in (mu,signed))
            or weight != fmpq(int(far['D_numerator']),1 << far['D_grid_bits'])
            or not weight > 0 or not kphys > 0
            or kphys != metadata_upper('trial_operator_norm_sq_upper')
            or rational(signed,'original_far_matrix_error') != metadata_upper('far_matrix_error')
            or lift['original_metadata_grid_bits'] != 512
            or eps < metadata_upper('complete_far_operator_allowance')):
        raise ValueError('Different original physical weight, amplification or complete allowance.')
    actions = matrix_in(metric['inverse_actions'],BITS)
    features = matrix_in([g[4:4+2*R] for g in target['generators']],BITS)
    kernel = matrix_in(target['kernel'],BITS)
    if ((actions.nrows(),actions.ncols(),features.nrows(),features.ncols(),kernel.nrows(),kernel.ncols())
            != (1263,2,1281,2*R,2*R,2*R)
            or len(cached['recovered_approximate_symbols']) != M
            or len(cauchy['finite_cauchy_source_outputs']) != M
            or not 1 <= len(cached['interpolated_source_coefficients']) <= 13
            or len(seed['normalized_polynomial']) != 64
            or len(cached['physical_left_kernel_action']) != 1
            or len(cached['physical_left_kernel_action'][0]) != 2*R
            or any(str(q) not in mu['real_hurwitz_kernels'] for q in range(2,188))):
        raise ValueError('Incomplete physical profile, source cache, kernel or saved inverse action.')
    projected,norm_caps = aggregate_good_profiles(features,positive,actions)
    mixed = arb(kappa)*matrix_in(cached['physical_left_kernel_action'],BITS)*projected
    samples = [interval_ball(x,BITS) for x in cauchy['finite_cauchy_source_outputs']]
    symbols = [interval_ball(x,BITS) for x in cached['recovered_approximate_symbols']]
    polynomial = arb_poly([interval_ball(x,BITS) for x in cached['interpolated_source_coefficients']])
    smooth = arb_poly([interval_ball(x,BITS) for x in seed['normalized_polynomial']])
    kernels = {q:interval_ball(v,BITS) for q,v in mu['real_hurwitz_kernels'].items()}
    anchor_row = next(v for v in far['moment_pairs'] if v['U_power'] == 2*R)
    anchor_zeta = next(v for v in far['zeta'] if v['p'] == 2*R+2)
    anchor = [interval_ball(v,BITS) for v in (anchor_row['U'],anchor_row['T'],anchor_zeta['interval'])]
    details,actual = [],[]
    for j in range(2):
        h = arb_mat([[projected[i,j]] for i in range(2*R)])
        values = list(h.entries())
        lcap = norm_caps[j]
        forward = scaled_forward(rational(signed,'sharp_forward_error_upper'),
            rational(signed,'lift_norm_upper'),lcap)
        mixed_error = scaled_forward(rational(signed,'mixed_cross_error_upper'),
            rational(signed,'lift_norm_upper'),lcap)
        quadratic,gram_error,broad2,_ = mixed_norm(h,kernel,eps,lcap)
        directional = directional_norm(values,freeze_upper(quadratic.upper()),
            freeze_upper(interval_ball(primal['bounds']['actual_symbol_budget'],BITS).upper()),
            fmpq(int(far['alpha_pair'][1]),1 << BITS),
            fmpq(int(far['alpha_pair'][1])-int(far['alpha_pair'][0]),1 << (BITS+1)),
            weight,fmpq(2*R,1 << BITS),rational(signed,'moment_gram_norm_upper'),
            forward,R,anchor,fmpq(M,M+1))
        ncap2 = min(broad2,directional['rebuilt_representer_energy_upper'])
        ncap = freeze_upper(arb(ncap2).sqrt())
        errors = coordinate_errors(values[:R],
            interval_ball(seed['uniform_symbol_error'],seed['analytic_grid_bits']).upper(),
            weight,rational(cached,'source_tail_norm_error_upper'),ncap,forward,
            upper(hankel,'principal_weighted_norm_upper'),upper(hankel,'smooth_remainder_weighted_norm_upper'),
            mixed_error,kappa)
        finite = finite_pairing(samples,symbols,values,M,kappa,weight)
        tail,tail_cap,oscillatory = cached_tail(polynomial,values,smooth,kernels,
            upper(mu,'inverse_phase_gap_weight_cap'),M,kappa,weight)
        allowance = errors['nonoscillatory_error_upper']+oscillatory
        value = finite+tail+mixed[0,j]+arb(0,arb(allowance))
        require_finite([value])
        actual.append(value)
        exact = dict(**errors,**directional,full_retained_norm_upper=lcap,
            physical_gram_error_upper=gram_error,broad_physical_norm2_upper=broad2,
            selected_physical_norm2_upper=ncap2,sharp_forward_error_upper=forward,
            oscillatory_error_upper=oscillatory,tail_coefficient_cap=freeze_upper(tail_cap),
            complete_error_upper=allowance)
        details.append(dict(profile=j,finite=envelope(finite,BITS),tail=envelope(tail,BITS),
            mixed=envelope(mixed[0,j],BITS),masked_physical_quadratic=envelope(quadratic,BITS),
            actual_pairing=envelope(value,BITS),
            exact_rationals={k:[str(v.p),str(v.q)] for k,v in exact.items()}))
    gram = matrix_in(metric['joint_inverse_metric'],BITS)
    known = arb_mat([[gram[0,0],gram[0,1]+actual[0]],
                     [gram[0,1]+actual[0],gram[1,1]+2*actual[1]]])
    result = dict(schema='ACTUAL_C27_GOOD_BLOCK_FAR_CROSS_V1',packet_sha256=PINS,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        N=1280,M=M,basis=row['basis'],positive_modes=positive,remaining_modes=remaining,
        source_integer_sha256=SOURCE_ID,source_bits=384,source_row=0,grid_bits=BITS,working_bits=ctx.prec,
        new_profiles=['P^-1*r_g*','P^-1*bJ_g*'],projected_far_features=matrix_out(projected,BITS),
        actual_far_cross=matrix_out(arb_mat([actual]),BITS),profile_details=details,
        complete_metric_without_unknown_far_energy=matrix_out(known,BITS),
        unknown_far_energy='z=f_g P^-1 f_g^* >= 0; no upper supplied',
        new_aggregate_profile_columns=2,new_masked_physical_quadratics=2,
        new_directional_norm_applications=2,new_finite_functionals=2,new_cached_tail_applications=2,
        new_mixed_batch_products=1,new_source_recovery_target_reconstruction_inverse_solve=0,
        original18_lift_queries=0,old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        complete_joint_metric_supplied=False,actual_center_sign_claim=False,
        actual_full_border_sign_claim=False,actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'Two saved inverse actions belong to the SAME1263-positive principal block of the rho32-v2 lower target; no general inverse or whole-center positivity is assumed.',
            'Aggregate the original112-dimensional far features before source-functionals and error bounds. The embedded retained profiles have zero remaining coordinates and their full Euclidean norms are bounded directly.',
            'Both the complete physical Gram allowance and accepted directional jet norm bound apply to each SAME profile. The signed kernel is not assumed PSD; no original W or kappa factor is duplicated.',
            'The uniform sharp forward and raw mixed allowances scale by the full retained norm. Original symbol/source/smooth/forward/mixed and Abel payments are all retained. Saved approximate symbols are not the actual symbol.',
            'The resulting pair is a joint containing box for two functionals of ONE actual source, not independently attainable endpoints. It encloses f_g P^-1 r_g* and f_g P^-1 bJ_g* with the existing kappa33/32 raw-row normalization.',
            'The complete joint metric equals the known matrix plus diag(0,z), z=f_g P^-1 f_g* >=0. The far inverse energy upper is missing; the known matrix is not claimed PSD or the complete metric. The old18-coordinate far row remains reused, not recomputed.',
            'All inherited common-form/domain/residual/source/lift and comparison premises remain. Existing methods only; no generic novelty, Lean proof, actual-center/full-border/cofinal positivity, Robin or RH claim.',
        ])
    print('new good-block far cross',write_packet(out,result),'seconds',result['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; no overwrite or retry.')
    ctx.prec = 1024
    consume(args.out)
