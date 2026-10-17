"""Condition the saved full far lift on two paid inverse-action observations.

Build only the NEW two-by-two and two-by-eighteen physical Gram blocks.
Reuse the old eighteen-by-eighteen block and both saved feature projections.
No source functional, target, inverse action, old producer or scan runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_correlated_sign import symmetric
from full_gamma_far_energy_budget import finite, projection_floor
from full_gamma_negative_index_joint import freeze_lower
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'structured-joint-2048-target.json.zlib': 'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
    'full-gamma-lift-projection.json.zlib': 'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'joint-gamma-6765-rho32-v2-actions.json.zlib': 'ff68d8350005be78daaa78a7e3f73ff7c09437dc9b1069825c5cc98b74c61962',
    'full-gamma-finite-good-block-metric.json.zlib': '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072',
    'full-gamma-good-block-far-cross.json.zlib': '39ddee410fa83bf4d35cb97f33bc9e63568dff8d9404dee331b0b3925336a858',
    'full-gamma-coherent-split-capacity.json.zlib': '6e4232d6079a4279511a67c13c6501d51a3b4fdaf410ccfee957adce29e1d7bc',
    'full-gamma-joint-far-row.json.zlib': '00b9e4bdd8a99da4d6fdafc696d313bd43e4060e1ba7d8cec316f7cc5bcb7707',
}
TARGET, LIFT, ACTIONS, METRIC, CROSS, ENERGY, FAR_ROW = PINS


def join_blocks(mm, mv, vv):
    n = vv.nrows()
    if (mm.nrows(), mm.ncols(), mv.nrows(), mv.ncols(), vv.ncols()) != (2, 2, 2, n, n):
        raise ValueError('Require complete two-observation and lifted blocks.')
    return arb_mat([[mm[i,j] if i < 2 and j < 2 else
                     mv[i,j-2] if i < 2 else
                     mv[j,i-2] if j < 2 else vv[i-2,j-2]
                     for j in range(n+2)] for i in range(n+2)])


def extend_physical_gram(kernel, hy, hv, yy, yv, vv, old_far, eps):
    """Retain the signed kernel, full allowance and old physical block."""
    n, rank = hv.ncols(), kernel.nrows()
    if (n < 1 or rank < 1 or kernel.ncols() != rank
            or (hy.nrows(),hy.ncols(),hv.nrows()) != (rank,2,rank)
            or (yy.nrows(),yy.ncols(),yv.nrows(),yv.ncols()) != (2,2,2,n)
            or (vv.nrows(),vv.ncols(),old_far.nrows(),old_far.ncols()) != (n,n,n,n)
            or not eps > 0):
        raise ValueError('Incompatible signed physical profiles or complete allowance.')
    finite([arb(eps)]+[x for m in (kernel,hy,hv,yy,yv,vv,old_far) for x in m.entries()])
    kernel, yy, vv, old_far = map(symmetric, (kernel,yy,vv,old_far))
    new_mm = symmetric(hy.transpose()*kernel*hy+arb(eps)*yy)
    new_mv = hy.transpose()*kernel*hv+arb(eps)*yv
    physical = join_blocks(new_mm,new_mv,old_far)
    euclidean = join_blocks(yy,yv,vv)
    finite(list(physical.entries())+list(euclidean.entries()))
    return physical, euclidean


def condition_profiles(physical, measured, gamma, old_row, raw_old_pairs=None):
    """Outer conditional supports, under the inherited exact PSD majorant.

The exact response column xi obeys xi*xi^T <= gamma*F. F is the
physical upper Gram, not the P-inverse metric. Only the observed2x2
pivots are numerically required. No full entry box is asserted PSD.
"""
    n = physical.nrows()-2
    if (n < 1 or physical.ncols() != n+2
            or (measured.nrows(),measured.ncols(),old_row.nrows(),old_row.ncols()) != (1,2,1,n)
            or not gamma > 0):
        raise ValueError('Incomplete common-source conditioning domains.')
    finite(list(physical.entries())+list(measured.entries())+list(old_row.entries())+[arb(gamma)])
    mm = arb_mat([[physical[i,j] for j in range(2)] for i in range(2)])
    mv = arb_mat([[physical[i,j+2] for j in range(n)] for i in range(2)])
    vv = arb_mat([[physical[i+2,j+2] for j in range(n)] for i in range(n)])
    floor = projection_floor(mm,measured)
    a, c, s = mm[0,0], mm[0,1], floor['metric_schur']
    inverse = arb_mat([[1/a+c*c/(a*a*s),-c/(a*s)],[-c/(a*s),1/s]])
    regression = inverse*mv
    schur = symmetric(vv-mv.transpose()*regression)
    means = measured*regression
    spent = floor['far_energy_projection_floor']
    spent_lower = max(fmpq(0),freeze_lower(spent.lower()))
    if spent_lower > gamma:
        raise ArithmeticError('Measured response box conflicts with common physical energy.')
    remaining = gamma-spent_lower
    rounded_pairs = matrix_out(old_row,BITS)[0]
    old_pairs = rounded_pairs if raw_old_pairs is None else raw_old_pairs
    if (len(old_pairs) != n or any(len(pair) != 2 for pair in old_pairs)
            or any(not int(broad[0]) <= int(pair[0]) <= int(pair[1]) <= int(broad[1])
                   for pair,broad in zip(old_pairs,rounded_pairs))):
        raise ValueError('Raw old grid intervals must lie within the enclosed old row.')
    support,intersections,details = [],[],[]
    shrinking = []
    for j in range(n):
        if schur[j,j].upper() < 0 or vv[j,j].upper() < 0:
            raise ArithmeticError('Conditional diagonal has no PSD physical realization.')
        variance = max(fmpq(0),freeze_upper(schur[j,j].upper()))
        radius = freeze_upper((arb(remaining)*arb(variance)).sqrt())
        interval = means[0,j]+arb(0,arb(radius))
        pair = envelope(interval,BITS)
        lo,hi = max(int(pair[0]),int(old_pairs[j][0])),min(int(pair[1]),int(old_pairs[j][1]))
        if lo > hi:
            raise ArithmeticError('Conditional support contradicts the saved same-source row.')
        intersections.append([str(lo),str(hi)])
        support.append(pair)
        old_width = int(old_pairs[j][1])-int(old_pairs[j][0])
        if hi-lo < old_width:
            shrinking.append(j)
        detail = dict(coordinate=j,conditional_radius_upper=[str(radius.p),str(radius.q)],
            old_width_units=str(old_width),intersection_width_units=str(hi-lo),
            conditional_variance=envelope(schur[j,j],BITS))
        if vv[j,j] > 0:
            detail['conditional_to_unconditional_variance_ratio'] = envelope(schur[j,j]/vv[j,j],BITS)
        details.append(detail)
    finite(list(regression.entries())+list(schur.entries())+list(means.entries())+[spent])
    return dict(observed_inverse=matrix_out(inverse,BITS),
        regression=matrix_out(regression,BITS),conditional_gram=matrix_out(schur,BITS),
        conditional_center=matrix_out(means,BITS),conditional_coordinate_support=[support],
        intersected_far_row=[intersections],coordinate_details=details,
        shrinking_original_coordinates=shrinking,
        scalar_bounds={'observed_first_pivot':envelope(a,BITS),
            'observed_schur_pivot':envelope(s,BITS),
            'observed_energy':envelope(spent,BITS),
            'observed_energy_fraction':envelope(spent/arb(gamma),BITS)},
        exact_rationals={'gamma':[str(gamma.p),str(gamma.q)],
            'observed_energy_lower':[str(spent_lower.p),str(spent_lower.q)],
            'remaining_energy_upper':[str(remaining.p),str(remaining.q)]})


def consume(out):
    started = time.perf_counter()
    guards = {Path(__file__).resolve():hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    for module in tuple(sys.modules.values()):
        p = getattr(module,'__file__',None)
        if p and Path(p).resolve().parent == HERE:
            p = Path(p).resolve();guards[p] = hashlib.sha256(p.read_bytes()).hexdigest()
    packets = {}
    for name,digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted profile input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    target,lift,actions,metric,cross,energy,far = packets.values()
    if [p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_STRUCTURED_TARGET_V1','ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_SCHUR_LIFT_V1','ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1',
            'ACTUAL_C27_GOOD_BLOCK_FAR_CROSS_V1','ACTUAL_C27_COHERENT_SPLIT_CAPACITY_V1',
            'ACTUAL_C27_JOINT_FAR_ROW_V1']:
        raise ValueError('Different common-source, physical Gram or inverse-action roles.')
    positive,remaining = metric['positive_modes'],metric['remaining_modes']
    if (len(positive) != 1263 or len(remaining) != 18 or sorted(positive+remaining) != list(range(1281))
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['remaining_modes'] != remaining for p in (lift,actions,cross,energy,far))
            or any(p['positive_modes'] != positive for p in (lift,actions,cross))
            or any((p['source_bits'],p['source_row'],p['source_integer_sha256']) != (384,0,SOURCE_ID)
                   for p in (metric,cross,energy,far))
            or lift['base_lift_sha256'] != PINS[ACTIONS]
            or lift['original_target_sha256'] != PINS[TARGET]
            or lift['base_target_sha256'] != actions['target_packet_sha256']
            or metric['packet_sha256']['joint-gamma-6765-rho32-v2-target.json.zlib'] != lift['base_target_sha256']
            or cross['packet_sha256'][METRIC] != PINS[METRIC]
            or cross['packet_sha256'][LIFT] != PINS[LIFT]
            or cross['packet_sha256'][TARGET] != PINS[TARGET]
            or energy['input_sha256'] != lift['input_sha256']
            or target['metadata']['input_sha256'] != lift['input_sha256']
            or (target['metadata']['N'],target['metadata']['M'],target['metadata']['last_mode']) != (1280,32768,2048)
            or any(p['basis'] != lift['basis'] for p in (metric,cross))
            or lift['original_metadata_grid_bits'] != 512
            or metric['joint_row_order'] != ['complete_unweighted_r','finite_bJ']
            or cross['new_profiles'] != ['P^-1*r_g*','P^-1*bJ_g*']
            or metric['complete_joint_metric_supplied'] or cross['complete_joint_metric_supplied']):
        raise ValueError('Different SAME source, old positive block, lift, basis or observation.')
    eps = fmpq(lift['exact_far_allowance'])
    original_eps = fmpq(int(target['metadata']['complete_far_operator_allowance'][1]),1 << 512)
    rational = lambda p,k:fmpq(*map(int,p['exact_rationals'][k]))
    kappa,source_energy = rational(energy,'kappa'),rational(energy,'actual_source_energy_cap')
    if kappa != fmpq(33,32) or not source_energy > 0 or not eps >= original_eps > 0:
        raise ValueError('Missing original allowance or SAME physical source-energy normalization.')
    y = matrix_in(metric['inverse_actions'],BITS)
    old_action = matrix_in(actions['inverse_action'],BITS)
    if (y.nrows(),y.ncols(),old_action.nrows(),old_action.ncols()) != (1263,2,1263,18):
        raise ValueError('Incomplete saved exact inverse-action enclosures.')
    yy = y.transpose()*y
    yv = -(y.transpose()*old_action)
    hy,hv = matrix_in(cross['projected_far_features'],BITS),matrix_in(lift['far_features_on_lift'],BITS)
    kernel = matrix_in(target['kernel'],BITS)
    vv,old_far = matrix_in(lift['lift_gram'],BITS),matrix_in(lift['far_gram_upper'],BITS)
    if (kernel.nrows(),hy.nrows(),hv.nrows(),hv.ncols()) != (112,112,112,18):
        raise ValueError('Incomplete original112-feature physical profiles.')
    physical,euclidean = extend_physical_gram(kernel,hy,hv,yy,yv,vv,old_far,eps)
    measured,old_row = matrix_in(cross['actual_far_cross'],BITS),matrix_in(far['actual_far_row'],BITS)
    try:
        conditional = condition_profiles(physical,measured,kappa*kappa*source_energy,
                                         old_row,far['actual_far_row'][0])
        conditional['conditioning_certified'] = True
    except ArithmeticError as error:
        if str(error) not in {'Missing positive first metric pivot.',
                              'Missing positive joint metric Schur pivot.'}:
            raise
        conditional = dict(conditioning_certified=False,
            conditioning_diagnostic=str(error),shrinking_original_coordinates=[],
            scalar_bounds={'observed_first_pivot':envelope(physical[0,0],BITS)})
    result = dict(schema='ACTUAL_C27_CONDITIONED_TWENTY_PROFILE_GRAM_V1',packet_sha256=PINS,
        source_sha256=guards[Path(__file__).resolve()],
        inherited_dependency_sha256={p.name:h for p,h in guards.items() if p.name not in PINS and p != Path(__file__).resolve()},
        source_integer_sha256=SOURCE_ID,source_bits=384,source_row=0,N=1280,M=32768,
        remaining_modes=remaining,positive_modes=positive,basis=lift['basis'],
        grid_bits=BITS,working_bits=ctx.prec,original_metadata_grid_bits=512,
        exact_far_allowance=str(eps),source_energy_upper=str(source_energy),kappa=str(kappa),
        profile_order=['iota_g P^-1 r_g*','iota_g P^-1 bJ_g*']+['V column '+str(j) for j in range(18)],
        joint_physical_gram=matrix_out(physical,BITS),joint_euclidean_gram=matrix_out(euclidean,BITS),
        measured_pair=matrix_out(measured,BITS),**conditional,
        verdict=('inconclusive_observed_physical_gram_pivots' if not conditional['conditioning_certified'] else
                 'conditional_coordinate_support_improved' if conditional['shrinking_original_coordinates']
                 else 'no_original_coordinate_box_improvement_joint_gram_retained'),
        new_euclidean_YY=1,new_euclidean_YV=1,new_physical_mm=1,new_physical_mv=1,
        reused_old_physical_vv=True,reused_old_feature_projections=True,
        new_two_observation_conditioning=1,new_source_far_target_inverse_solve_queries=0,
        old_scientific_producers_or_fixtures_replayed=0,source_direction_parameter_precision_cutoff_scan=False,
        full_entry_box_asserted_PSD=False,unseen_far_inverse_energy_upper_supplied=False,
        actual_center_sign_claim=False,actual_full_border_sign_claim=False,actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'One original384-bit row0 source, SAME1263-positive block, exact enclosed old18-lift and two paid inverse-action observations. No old functional or inverse-action producer is replayed.',
            'A=[iota_g Y,V]. Reuse H*Y,H*V and the old V*Fphys*V block; build only new observed and mixed physical blocks, including ALL epsilon*A*A blocks. Retain the complete original512-metadata allowance.',
            'The signed112-kernel is not assumed PSD. Under the inherited physical majorant R_F*W_F R_F<=H K H*+eps I, the exact joint physical Gram F is PSD. Its full interval entry box is not asserted PSD.',
            'For the same physical source norm upper E, the column xi=(m,fV)^T obeys xi*xi^T<=kappa^2 E F. Original W is already included in F and E; kappa33/32 appears once in the response energy.',
            'Positive observed2x2 pivots justify its explicit Schur regression. For the actual measured m, the centered y residual is bounded by (gamma-m*Fmm^-1*m) times the conditional Gram. A lower enclosure of the spent energy gives a safe remaining-energy upper.',
            'Coordinate supports and intersections use the original saved640-grid row intervals for ONE common realization, not independently attainable endpoints. No width improvement is credited merely against a round-trip enlargement. Failure to shrink an old coordinate box does not exclude improvement of a correlated core quadratic.',
            'Uncertified observed2x2positive pivots retain the new joint Gram with an explicit inconclusive diagnostic; they do not trigger another contraction or higher-precision retry.',
            'This supplies a common-source observation geometry, not a certificate for the indefinite core or determinant. The unseen far P-inverse energy z remains uncontrolled above; a sufficient direct joint comparison remains an alternative.',
            'All inherited common-form, domain, residual-map, physical-majorant and source/lift premises remain. Existing Schur/ellipsoid methods only; no generic novelty or Lean claim. Actual-center/full-border/cofinal positivity, Robin and RH remain open.',
        ])
    for path,digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError('Input/source changed during new profile application: '+str(path))
    print('saved conditioned20-profile Gram',write_packet(out,result),'seconds',result['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; no overwrite or precision retry.')
    ctx.prec=1024
    consume(args.out)
