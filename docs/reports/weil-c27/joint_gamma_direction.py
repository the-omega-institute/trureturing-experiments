"""Distinguish a joint Young-bound obstruction from the ideal updated target.

Consume the saved small-block lift; no dense solve or row production is
repeated. On an exact full dyadic direction, bound the ideal Gamma-updated
inverse-residual target above and below using its actual far covariance.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, N, M
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, original_target
from structured_inverse_residual_interval import envelope
from structured_residual_schur_finish import load_target, reconstruct_target


def run(args):
    started=time.monotonic()
    old,oh,metadata_bits=original_target(args)
    finite,fh=read_json(args.finite)
    updated,prefix,remaining,uh,pivot_raw=load_target(args.lift_target,args.lift_pivots)
    actions,ah=read_json(args.actions)
    if (actions['schema'],actions['target_packet_sha256'],actions['pivot_result_sha256'],
            actions['positive_modes'],actions['remaining_modes']) != (
            'ACTUAL_C27_SCHUR_LIFT_V1',uh,hashlib.sha256(pivot_raw).hexdigest(),prefix,remaining):
        raise ValueError('Require the exact saved small-block lift.')
    if (updated['metadata']['original_target_sha256'],updated['metadata']['joint_boundary_action_sha256'],
            finite['target_packet_sha256'],args.mode in remaining) != (oh,fh,oh,True):
        raise ValueError('The old and updated target must share their actual source and thin action.')
    lift=matrix_in(actions['inverse_action'],actions['grid_bits'])
    if (lift.nrows(),lift.ncols()) != (len(prefix),len(remaining)):
        raise ValueError('Incomplete saved lift.')
    integers=[0]*(N+1)
    integers[args.mode]=1<<args.direction_bits
    column=remaining.index(args.mode)
    for i,n in enumerate(prefix):
        integers[n]=trial_integer(-lift[i,column],args.direction_bits)
    vq=[fmpq(x,1<<args.direction_bits) for x in integers]
    v=arb_mat([[arb(x)] for x in vq])
    norm=sum((x*x for x in vq),fmpq(0))
    old_target=reconstruct_target(old)
    qold=(v.transpose()*old_target*v)[0,0]
    oldG=matrix_in(old['generators'],old['grid_bits'])
    K=matrix_in(old['kernel'],old['grid_bits'])
    P=arb_mat([[oldG[n,4+k] for n in range(N+1)] for k in range(K.nrows())])
    pv=P*v
    allowance=interval_ball(old['metadata']['complete_far_operator_allowance'],metadata_bits).upper()
    far_upper=(pv.transpose()*K*pv)[0,0]+allowance*arb(norm)
    if not far_upper > 0:
        raise ArithmeticError('Require a positive upper enclosure for the actual far residual energy.')
    beta=matrix_in(finite['finite_beta'],finite['grid_bits'])*v
    A=matrix_in(finite['finite_gram'],finite['grid_bits'])
    delta=interval_ball(finite['far_gram_trace_bound'],finite['grid_bits']).upper()
    S=arb_mat(A)
    for k in range(S.nrows()):
        S[k,k] += 1
    solution=S.solve(beta,algorithm='precond')
    zintegers=[trial_integer(solution[k,0],args.direction_bits) for k in range(S.nrows())]
    zq=[fmpq(x,1<<args.direction_bits) for x in zintegers]
    z=arb_mat([[arb(x)] for x in zq])
    znorm=arb(sum((x*x for x in zq),fmpq(0)))
    error=(delta*far_upper.upper()).sqrt()
    beta_norm=sum((beta[k,0]**2 for k in range(beta.nrows())),arb(0)).sqrt()
    gain_lower=(2*z.transpose()*beta-z.transpose()*S*z)[0,0]-2*znorm.sqrt()*error-delta*znorm
    gain_upper=(beta.transpose()*solution)[0,0]+2*beta_norm*error+error**2
    # L_old differs from the ideal old target only by the paid far enclosure.
    # Its actual discrepancy lies between0 and2*a*||v||².
    margin_lower=qold+gain_lower
    margin_upper=qold+2*allowance*arb(norm)+gain_upper
    rho=fmpq(updated['metadata']['common_young_parameter'])
    Srho=arb_mat(S)
    for k in range(Srho.nrows()):
        Srho[k,k] += (1+1/arb(rho))*delta
    joint_quadratic=qold+(beta.transpose()*Srho.solve(beta,algorithm='precond'))[0,0]-arb(rho)*far_upper
    verdict=('ideal-updated-fixed-trial-rejected' if margin_upper < 0 else
             'this-direction-repaired' if margin_lower > 0 else 'inconclusive')
    result=dict(schema='ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1',N=N,M=M,
                original_target_sha256=oh,updated_target_sha256=uh,finite_action_sha256=fh,
                lift_actions_sha256=ah,pivot_result_sha256=hashlib.sha256(pivot_raw).hexdigest(),
                selected_remaining_mode=args.mode,basis=old['metadata']['basis'],
                direction_dyadic_bits=args.direction_bits,direction_dyadic_integers=[str(x) for x in integers],
                direction_norm_sq=str(norm),gain_trial_dyadic_bits=args.direction_bits,
                gain_trial_dyadic_integers=[str(x) for x in zintegers],
                working_bits=args.precision,grid_bits=args.output_bits,
                original_metadata_grid_bits=metadata_bits,
                old_lower_target_quadratic=envelope(qold,args.output_bits),
                actual_far_residual_energy_upper=envelope(far_upper,args.output_bits),
                gamma_gain_lower=envelope(gain_lower,args.output_bits),gamma_gain_upper=envelope(gain_upper,args.output_bits),
                ideal_updated_margin_lower=envelope(margin_lower,args.output_bits),
                ideal_updated_margin_upper=envelope(margin_upper,args.output_bits),
                common_young_target_quadratic=envelope(joint_quadratic,args.output_bits),
                verdict=verdict,elapsed_seconds=time.monotonic()-started,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                conditions=['Same actual common-form/domain and complete exact trial Y=E F^-1 B.',
                            'Existing Gamma boundary/diagonal inverse comparisons and actual far Gram enclosure.',
                            'Sound reused full target, exact saved lift and new complete finite thin action.'],
                scope='Upper and lower estimates concern the ideal inverse-residual target with this fixed Y and W_G, '
                      'not an upper bound for the actual retained Schur form. A negative upper rejects this fixed trial target, '
                      'not other trials on E, other majorants, actual Weil positivity or RH. No production or dense solve rerun.')
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(verdict,'Young target',joint_quadratic.str(24),'ideal margin bounds',margin_lower.str(24),margin_upper.str(24),
          'seconds',time.monotonic()-started,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--target',type=Path,default=HERE/'structured-joint-2048-target.json.zlib')
    p.add_argument('--pivots',type=Path,default=HERE/'structured-joint-2048-pivots.json')
    p.add_argument('--finite',type=Path,default=HERE/'joint-gamma-6765-finite.json.zlib')
    p.add_argument('--lift-target',type=Path,required=True)
    p.add_argument('--lift-pivots',type=Path,required=True)
    p.add_argument('--actions',type=Path,required=True)
    p.add_argument('--mode',type=int,default=0)
    p.add_argument('--direction-bits',type=int,default=256)
    p.add_argument('--precision',type=int,default=1024)
    p.add_argument('--output-bits',type=int,default=640)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists() or a.precision < max(640,a.output_bits,a.direction_bits)+192:
        p.error('Use a fresh output and at least192 working guard bits.')
    ctx.prec=a.precision
    run(a)


if __name__ == '__main__':
    main()
