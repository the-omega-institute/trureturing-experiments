"""Consume the saved nineteen channels with the accepted uniform surplus.

One fixed whole-space Young inverse gain, followed by a necessary saved
full-direction screen. Only a passing screen permits a NEW19-RHS query
on the saved positive block and its fully paid small Schur complement.
No accepted scientific producer, old lift solve or fixture is replayed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq
from diagonal_floor_witness import interval_ball
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope
from structured_residual_schur_finish import load_target, reconstruct_target

HERE = Path(__file__).resolve().parent
BITS = 640
ORIGINAL = 'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22'
BASIS = '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
TARGET = 'joint-gamma-6765-rho32-v2-target.json.zlib'
PIVOTS = 'joint-gamma-6765-rho32-v2-pivots.json'
FINITE = 'joint-gamma-6765-finite.json.zlib'
LIFT = 'joint-gamma-6765-rho32-v2-actions.json.zlib'
CORE = 'joint-gamma-6765-rho32-v2-schur.json'
DIRECTION = 'joint-gamma-6765-direction.json'
PILOT = 'inverse-residual-pilot.json.zlib'
ROW = 'propagated-prime-row.json.zlib'
PINS = {
    TARGET: '7546a915368b4a40e0e41ea0ea693403c239598c7a6e96b7b8b07e67e22e98c6',
    PIVOTS: 'd441b9f9fc49338a121430a501043b76afcf7dffa9ec09c8d672aa9a2358b3f0',
    FINITE: '94aec57dacacad6bae1698a7e643514e3dace976fddbcc0008b2b7fd5aab447c',
    LIFT: 'ff68d8350005be78daaa78a7e3f73ff7c09437dc9b1069825c5cc98b74c61962',
    CORE: '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    DIRECTION: 'ccdb3360923f4fd07086871c5d02512a71301143376db8fd54aa23a3797152cc',
    PILOT: '6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b',
    ROW: '935b2c0707c3187fa652730b1e40d5937a4e071158b7cb9f1f4a6a78f8459ce6',
}


def rational(q):
    return arb(fmpq(q.numerator, q.denominator))


def pair(q):
    return [str(q.numerator), str(q.denominator)]


def symmetric(a):
    if a.nrows() != a.ncols():
        raise ValueError('A symmetric enclosure must be square.')
    return (a+a.transpose())/2


def require_shape(a, rows, columns):
    if (a.nrows(), a.ncols()) != (rows, columns):
        raise ValueError('Incomplete saved thin matrix.')
    if any(not x.is_finite() for x in a.entries()):
        raise ValueError('Nonfinite matrix; stop without a retry.')


def leading_positive(a):
    """Use FLINT determinants and the existing degree-dependent export grid."""
    require_shape(a, a.nrows(), a.nrows())
    minors = []
    for k in range(1, a.nrows()+1):
        block = arb_mat([[a[i,j] for j in range(k)] for i in range(k)])
        bits = BITS*k+64
        value = envelope(block.det(), bits)
        sign = (1 if value is not None and int(value[0]) > 0 else
                -1 if value is not None and int(value[1]) < 0 else 0)
        minors.append(dict(order=k, grid_bits=bits, interval=value, sign=sign))
        if sign != 1:
            return ('negative-leading-minor' if sign == -1 else 'inconclusive-leading-minor'), minors
    return 'positive-definite', minors


def new_denominator(a, s, ell, gap):
    """A CA+(ell/g)CA A CA; no new unweighted-channel Gram is produced."""
    if s <= 0 or ell <= 0 or gap <= 0:
        raise ValueError('Require positive fixed charge,weight bound and surplus.')
    a = symmetric(a)
    ca = a/rational(s)
    for k in range(ca.nrows()):
        ca[k,k] += 1
    return symmetric(a*ca+rational(ell/gap)*(ca*a*ca))


def classify_screen(old, gain):
    value = old+gain
    if not value.is_finite():
        return 'inconclusive-screen', value
    return ('fixed-nineteen-target-negative' if value.upper() < 0 else
            'necessary-direction-positive-only' if value.lower() > 0 else
            'inconclusive-screen'), value


def checked_solve(a, rhs):
    answer = a.solve(rhs, algorithm='precond')
    residual = a*answer-rhs
    if any(not x.is_finite() or not x.contains(0) for x in residual.entries()):
        raise ArithmeticError('New action failed the zero-residual check; no retry.')
    return answer


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    target, positive, remaining, target_hash, pivot_raw = load_target(HERE/TARGET, HERE/PIVOTS)
    finite, lift, core, direction, pilot, row = (packets[k] for k in
                                               [FINITE,LIFT,CORE,DIRECTION,PILOT,ROW])
    meta = target['metadata']
    if (target_hash != PINS[TARGET] or hashlib.sha256(pivot_raw).hexdigest() != PINS[PIVOTS]
            or len(positive) != 1263 or len(remaining) != 18
            or meta['original_target_sha256'] != ORIGINAL
            or meta['joint_boundary_action_sha256'] != PINS[FINITE]
            or meta['common_young_parameter'] != '1/32'
            or meta['positive_boundary_bins'] != finite['bins'] or len(finite['bins']) != 19
            or target['grid_bits'] != BITS or target['producer_working_bits'] != 1024
            or (meta['N'],meta['M'],meta['last_mode'],meta['basis']) != (1280,32768,2048,BASIS)
            or finite['schema'] != 'ACTUAL_C27_JOINT_GAMMA_ACTION_V1'
            or finite['target_packet_sha256'] != ORIGINAL
            or (finite['N'],finite['M'],finite['last_mode'],finite['basis']) != (1280,32768,2048,BASIS)
            or finite['grid_bits'] != BITS or finite['working_bits'] != 1024
            or finite['input_sha256'] != meta['input_sha256']
            or meta['input_sha256'][PILOT] != PINS[PILOT]
            or direction['schema'] != 'ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1'
            or direction['original_target_sha256'] != ORIGINAL
            or direction['updated_target_sha256'] != target_hash
            or direction['finite_action_sha256'] != PINS[FINITE]
            or direction['pivot_result_sha256'] != PINS[PIVOTS]
            or (direction['N'],direction['M'],direction['basis']) != (1280,32768,BASIS)
            or direction['grid_bits'] != BITS or direction['direction_dyadic_bits'] != 256
            or direction['selected_remaining_mode'] != 0
            or len(direction['direction_dyadic_integers']) != 1281
            or any(type(n) is not str for n in direction['direction_dyadic_integers'])
            or pilot['schema'] != 'ACTUAL_C27_INVERSE_RESIDUAL_PILOT_V1'
            or (pilot['N'],pilot['M'],pilot['basis'],pilot['middle_weight_bits']) != (1280,32768,BASIS,256)
            or row['schema'] != 'ACTUAL_C27_CANONICAL_PROPAGATED_PRIME_ROW_V1'
            or row['canonical_weight'] != 'q1=C_p*(1+(29/50)*(x/L)^2)'
            or row['propagation_applications'] != 1
            or row['verdict'] != 'actual_slack_and_new_eligible_coefficient'
            or row['actual_retained_or_cofinal_positivity_claim'] is not False):
        raise ValueError('Different actual target,channel,charge,direction,basis or uniform surplus.')
    for packet, schema in [(lift,'ACTUAL_C27_SCHUR_LIFT_V1'),
                           (core,'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1')]:
        if (packet['schema'] != schema or packet['target_packet_sha256'] != target_hash
                or packet['pivot_result_sha256'] != PINS[PIVOTS]
                or packet['grid_bits'] != BITS or packet['remaining_modes'] != remaining):
            raise ValueError('Different saved lift or small complement.')
    if (lift['positive_modes'] != positive or core['reused_positive_coordinate_count'] != 1263
            or core['input_sha256'] != meta['input_sha256'] or core['basis'] != BASIS):
        raise ValueError('Different good-block coordinates or actual inputs.')
    gap = F(*map(int,row['exact_rationals']['actual_surplus_floor_lower']))
    if gap != F(16641320634704717921,46116860184273879040):
        raise ValueError('Use the fixed accepted surplus,not a parameter selection.')
    nums = pilot['middle_weight_numerators']
    if len(nums) != 31488 or any(type(n) is not str for n in nums):
        raise ValueError('Need all1281..32768 weights,including selectedE.')
    weights = [F(int(n),1 << pilot['middle_weight_bits']) for n in nums]
    if min(weights) <= 0:
        raise ValueError('Positive weights required before inversion.')
    ell = 1/min(weights)
    # Match EXACTLY the original Young decoder and its working precision.
    delta = interval_ball(finite['far_gram_trace_bound'],BITS).upper()
    if not delta.is_finite() or not delta.is_exact() or delta < 0:
        raise ValueError('No exact original far charge; stop.')
    delta_q = F(str(delta.fmpq()))
    s = 1+33*delta_q
    old_s = 1+(1+1/arb(fmpq(1,32)))*delta
    if not old_s.is_exact() or not (old_s-rational(s)).is_zero():
        raise ValueError('Original Young charge cannot be preserved exactly.')
    a = matrix_in(finite['finite_gram'],BITS)
    beta = matrix_in(finite['finite_beta'],BITS)
    require_shape(a,19,19)
    require_shape(beta,19,1281)
    d = new_denominator(a,s,ell,gap)
    result = dict(schema='ACTUAL_C27_NINETEEN_UNIFORM_JOINT_GAIN_V1',
        packet_sha256=PINS,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        original_target_sha256=ORIGINAL,base_target_sha256=target_hash,
        N=1280,M=32768,last_mode=2048,basis=BASIS,grid_bits=BITS,working_bits=ctx.prec,
        common_young_parameter='1/32',channels=19,positive_modes=positive,remaining_modes=remaining,
        exact_rationals=dict(g=pair(gap),ell_all_middle=pair(ell),delta_ch=pair(delta_q),s=pair(s)),
        denominator=matrix_out(d,BITS),fixed_joint_gain_applications=1,
        new_positive_block_reconstruction_queries=0,new_good_block_solves=0,new_good_block_rhs=0,
        old_scientific_producers_or_fixtures_replayed=0,old18_rhs_lift_solve_replayed=False,
        norm_weight_source_direction_parameter_cutoff_precision_scan=False,
        actual_retained_or_cofinal_positivity_claim=False,actual_rh_proof=False,
        conditional_fixed_lower_target_positive=False,
        conditions=[
            'Same actual common form/domain,original exact trialY,retained orthonormal basis and sound-cache premises remain; no Lean verification.',
            'Original19-channel A uses ALL middle rows1281..32768 includingE,whereas beta can omitE since R_E=0. Original rho=1/32 and original delta_ch decoder at1024 are retained.',
            'Full-space Young V>=T^-1 includes nonzeroE vectors. H>=T+gI>=V^-1+gI gives Kg=V-(V^-1+gI)^-1>0. No false operator monotonicity is assumed for t^-1-(t+g)^-1.',
            'Original signed far core and complete allowance remain paid: Lstar=Q-R*V R-Efar,Efar>=0. Thus S_actual>=Lstar+R*KgR>=Lstar+beta*Dg^-1 beta when Dg>0.',
            'Dg=A CA+(ell/g)CA A CA bounds the dual observation Gram using U*U<=ell A. No unweighted-channel Gram or arbitrary missing residual energy is invented.',
            'Negative saved-direction upper rejects ONLY this fixed sufficient19-channel lower target. Positive lower passes ONLY one necessary screen and must be followed by the paid good-block query.',
            'If queried,Omega=beta_G P^-1 beta_G* is included in Dg+Omega. SavedX supplies betaV but cannot replace the new19-RHS inverse action.',
            'No additive combination with separately paid row gains,no retry of stopped norm-contraction or rank-one endpoint families,no actual retained/cofinal Weil,Robin orRH sign conclusion.',
        ])

    def finish(verdict):
        for name,digest in PINS.items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Input changed during fixed consumption: '+name)
        for name,digest in meta['input_sha256'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Inherited input changed: '+name)
        result.update(verdict=verdict,elapsed_seconds=time.perf_counter()-started)
        digest = write_packet(out,result)
        print(verdict,digest,'seconds',result['elapsed_seconds'],flush=True)

    status, minors = leading_positive(d)
    result['denominator_status'],result['denominator_minors'] = status,minors
    if status != 'positive-definite':
        finish('stop-denominator-'+status)
        return
    v = arb_mat([[arb(fmpq(int(n),1 << 256))] for n in direction['direction_dyadic_integers']])
    u = beta*v
    action = checked_solve(d,u)
    old = interval_ball(direction['common_young_target_quadratic'],BITS)
    gain = (u.transpose()*action)[0,0]
    status, value = classify_screen(old,gain)
    result['screen'] = dict(status=status,direction_bits=256,direction_packet_sha256=PINS[DIRECTION],
        direction_coordinates=1281,old_target=envelope(old,BITS),gain=envelope(gain,BITS),
        updated_target=envelope(value,BITS),channel_projection=matrix_out(u,BITS),
        denominator_action=matrix_out(action,BITS))
    if status != 'necessary-direction-positive-only':
        finish(status)
        return
    print('necessary direction passes; proceed to NEW19-RHS good-block query',flush=True)
    matrix = reconstruct_target(target)
    p = arb_mat([[matrix[i,j] for j in positive] for i in positive])
    del matrix
    bg = arb_mat([[beta[i,j] for j in positive] for i in range(19)])
    bb = arb_mat([[beta[i,j] for j in remaining] for i in range(19)])
    x = matrix_in(lift['inverse_action'],BITS)
    sigma = matrix_in(core['reduced_schur_intervals'],BITS)
    require_shape(x,1263,18)
    require_shape(sigma,18,18)
    result['new_positive_block_reconstruction_queries'] = 1
    query_started = time.perf_counter()
    inverse_action = checked_solve(p,bg.transpose())
    result.update(new_good_block_solves=1,new_good_block_rhs=19,
                  new_good_block_solve_seconds=time.perf_counter()-query_started)
    omega = symmetric(bg*inverse_action)
    bv = bb-bg*x
    paid_d = symmetric(d+omega)
    pd_status,pd_minors = leading_positive(paid_d)
    result['good_block'] = dict(inverse_actions=matrix_out(inverse_action,BITS),
        inverse_metric=matrix_out(omega,BITS),lifted_beta=matrix_out(bv,BITS),
        paid_denominator=matrix_out(paid_d,BITS),paid_denominator_status=pd_status,
        paid_denominator_minors=pd_minors)
    if pd_status != 'positive-definite':
        finish('stop-paid-denominator-'+pd_status)
        return
    small_gain = bv.transpose()*checked_solve(paid_d,bv)
    updated_core = symmetric(sigma+small_gain)
    status,minors = leading_positive(updated_core)
    result['updated_complement'] = matrix_out(updated_core,BITS)
    result['updated_complement_status'],result['updated_complement_minors'] = status,minors
    result['conditional_fixed_lower_target_positive'] = status == 'positive-definite'
    finish('conditional-fixed1281-target-positive' if status == 'positive-definite' else
           'stop-updated-complement-'+status)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; no overwrite or retry.')
    ctx.prec = 1024
    consume(args.out)
