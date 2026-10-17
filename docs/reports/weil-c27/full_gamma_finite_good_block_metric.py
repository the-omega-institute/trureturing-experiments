"""Enclose the missing joint good-block metric on the SAME finite rows.

The two NEW right sides are complete r and finite bJ on the1263-positive
block. Reuse the old target representation and positive-pivot evidence.
No old18-column lift solve, scientific supplier or generic inverse runs.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb_mat, ctx
from gamma_boundary_update import read_json
from full_gamma_finite_retained_row import export
from full_gamma_unweighted_retained_row import BITS, SOURCE_ID
from joint_gamma_boundary import matrix_in, write_packet
from structured_residual_schur_finish import load_target, reconstruct_target

HERE=Path(__file__).resolve().parent
TARGET='joint-gamma-6765-rho32-v2-target.json.zlib'
PIVOTS='joint-gamma-6765-rho32-v2-pivots.json'
ROW='full-gamma-unweighted-retained-row.json.zlib'
FINITE='full-gamma-finite-retained-row.json.zlib'


def symmetric_joint_energy(rows, solution):
    """Real-row Gram; both cross evaluations enclose the same actual entry."""
    if rows.nrows()!=2 or (solution.nrows(),solution.ncols())!=(rows.ncols(),2):
        raise ValueError('Require two joint rows and matching inverse actions.')
    gram=rows*solution
    cross=(gram[0,1]+gram[1,0])/2
    return arb_mat([[gram[0,0],cross],[cross,gram[1,1]]])


def consume(out):
    started=time.perf_counter()
    target,positive,remaining,target_hash,pivot_raw=load_target(HERE/TARGET,HERE/PIVOTS)
    r,rhash=read_json(HERE/ROW)
    finite,bhash=read_json(HERE/FINITE)
    if (r['schema'],finite['schema']) != (
            'ACTUAL_C27_COMPLETE_UNWEIGHTED_RETAINED_ROW_V1','ACTUAL_C27_FULL_FINITE_ROW_V1'):
        raise ValueError('Different unweighted/finite-row roles.')
    if (len(positive)!=1263 or len(remaining)!=18
            or r['base_target_sha256']!=target_hash
            or target['metadata']['original_target_sha256']!=r['original_target_sha256']
            or target['metadata']['input_sha256']!=r['input_sha256']
            or finite['input_sha256']!=r['input_sha256']
            or finite['packet_sha256']['joint-gamma-6765-rho32-v2-actions.json.zlib']!=r['base_lift_sha256']
            or any(v['positive_modes']!=positive or v['remaining_modes']!=remaining
                or v['retained_modes']!=list(range(1281))
                or (v['source_bits'],v['source_row'],v['source_integer_sha256'],v['grid_bits'])
                    !=(384,0,SOURCE_ID,BITS) for v in (r,finite))
            or r['basis']!=finite['basis'] or r['basis']!=target['metadata']['basis']
            or not r['complete_for_original_finite_source'] or r['complete_b_row_supplied']
            or r['actual_retained_sign_claim'] or r['actual_full_border_sign_claim'] or r['actual_rh_proof']):
        raise ValueError('Different actual source, base, coordinate partition or finite scope.')
    unweighted=matrix_in(r['complete_unweighted_row'],BITS)
    finite_row=matrix_in(finite['finite_common_row'],BITS)
    if any((v.nrows(),v.ncols())!=(1,1281) for v in (unweighted,finite_row)):
        raise ValueError('Incomplete joint retained rows.')
    rows=arb_mat([[v[0,j] for j in positive] for v in (unweighted,finite_row)])
    matrix=reconstruct_target(target)
    P=arb_mat([[matrix[i,j] for j in positive] for i in positive])
    reconstructed=time.perf_counter()
    print('same positive block reconstructed for NEW two-RHS query',flush=True)
    solution=P.solve(rows.transpose(),algorithm='precond')
    solved=time.perf_counter()
    residual=P*solution-rows.transpose()
    if any(not x.is_finite() or not x.contains(0) for x in residual.entries()):
        raise ArithmeticError('New inverse-action enclosure fails its zero-residual check; no retry.')
    gram=symmetric_joint_energy(rows,solution)
    determinant=gram[0,0]*gram[1,1]-gram[0,1]**2
    if any(not x.is_finite() for x in gram.entries()) or any(
            x.upper()<0 for x in (gram[0,0],gram[1,1],determinant)):
        raise ArithmeticError('New joint metric contradicts its inherited PSD identity.')
    result=dict(schema='ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1',
        packet_sha256={TARGET:target_hash,PIVOTS:hashlib.sha256(pivot_raw).hexdigest(),ROW:rhash,FINITE:bhash},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        N=1280,M=32768,basis=r['basis'],positive_modes=positive,remaining_modes=remaining,
        source_bits=384,source_row=0,source_integer_sha256=SOURCE_ID,
        grid_bits=BITS,working_bits=ctx.prec,joint_row_order=['complete_unweighted_r','finite_bJ'],
        joint_inverse_metric=export(gram),inverse_actions=export(solution),
        new_positive_block_reconstruction_queries=1,new_inverse_action_solves=1,new_rhs_columns=2,
        old18_rhs_lift_solve_replayed=False,new_general_inverse_calls=0,new_source_queries=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        complete_b_row_supplied=False,complete_joint_metric_supplied=False,
        actual_retained_sign_claim=False,actual_full_border_sign_claim=False,actual_rh_proof=False,
        reconstruction_seconds=reconstructed-started,solve_seconds=solved-reconstructed,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'P is the SAME1263-positive principal block of the rho32-v2 lower target, with the inherited accepted positive-pivot certificate. It is not a new proof of actual full-center positivity.',
            'The newly solved columns are P^-1*r_g* and P^-1*bJ_g*, not the archived1263x18 lift action. Reconstructing the saved target is a necessary query input to this new action, not an old target-producer or Schur/lift-solve replay.',
            'The joint metric retains the cross term for all real combinations of these two SAME rows. r is complete for the original finite w, whereas bJ excludes the actual far mixed contribution. The complete joint metric still needs that contribution and its correlations.',
            'Future complete-row estimates may reuse the saved two inverse actions. A finite bJ energy cannot be substituted for complete b energy. No retained/full-border/cofinal positivity, Robin, RH, generic novelty or Lean proof is claimed.',
        ])
    digest=write_packet(out,result)
    print('finite joint good-block metric',digest,'seconds',result['elapsed_seconds'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists() or args.out.suffix!='.zlib':
        parser.error('Require a fresh .zlib output; no overwrite or retry.')
    ctx.prec=1024
    consume(args.out)
