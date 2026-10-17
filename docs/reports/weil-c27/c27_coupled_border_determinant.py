"""Classify the saved joint lower target with one uneliminated border.

Reuse the saved19-RHS action, paid denominator, lifted beta and old
one-negative core. One NEW37x37 FLINT determinant at unchanged1024 bits
avoids the explicit inverse that produced an inconclusive small core.
No scientific supplier, old solve, determinant or fixture is replayed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import zlib

import flint
from flint import arb_mat, ctx
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
JOINT = 'nineteen-uniform-joint-gain.json.zlib'
CORE = 'joint-gamma-6765-rho32-v2-schur.json'
INERTIA = 'joint-gamma-6765-rho32-v2-rank-floor.json'
BASE = '7546a915368b4a40e0e41ea0ea693403c239598c7a6e96b7b8b07e67e22e98c6'
PIVOT = 'd441b9f9fc49338a121430a501043b76afcf7dffa9ec09c8d672aa9a2358b3f0'
PINS = {JOINT:'bc0c6f4cda73f7ace24a7c07956a8323feb529e74b140ffc96b7b4b739658455',
        CORE:'13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
        INERTIA:'f4ceb4b526de202b540089f7532ac3aee71e7013df98a7dfe0ab3975345cdf89'}
HELPERS = {'joint_gamma_boundary.py':'10021a4faf837683e8430a7124241f3145e1bcfdce42aca3f817cdd4a491ea07',
           'structured_inverse_residual_interval.py':'eae66659370f74c8a4c7adbf8df59755c8aa25d532a9652a3c16e7f2da10567c',
           'diagonal_floor_witness.py':'5942cdf0779e13d4ba0763122c311e4fec730ed1660e9fe01ff89eb8db16c6a6'}


def certified_sign(pair):
    if (not isinstance(pair,list) or len(pair)!=2 or
            any(type(x) is not str for x in pair)):
        return 0
    lo,hi = map(int,pair)
    if lo > hi:
        raise ValueError('Reversed determinant enclosure.')
    return 1 if lo > 0 else -1 if hi < 0 else 0


def require_certificate(rows, signs, pair_key, sign_key):
    """Read saved certificates, without recomputing their determinants."""
    if len(rows) != len(signs):
        raise ValueError('Incomplete accepted leading-minor certificate.')
    for k,(row,sign) in enumerate(zip(rows,signs),1):
        if (row['order'] != k or row['grid_bits'] != BITS*k+64
                or row[sign_key] != sign or certified_sign(row[pair_key]) != sign):
            raise ValueError('Different accepted positive/inertia certificate.')


def coupled_border(d, b, sigma):
    """Retain the shared B in BOTH off-diagonal blocks before elimination."""
    n,m = d.nrows(),sigma.nrows()
    if (n <= 0 or m <= 0 or d.ncols()!=n or sigma.ncols()!=m
            or (b.nrows(),b.ncols())!=(n,m)):
        raise ValueError('Incomplete coupled-block coordinates.')
    if any(not x.is_finite() for a in (d,b,sigma) for x in a.entries()):
        raise ValueError('Nonfinite saved enclosure; no alternate query.')
    d,sigma = (d+d.transpose())/2,(sigma+sigma.transpose())/2
    z = arb_mat(n+m,n+m)
    for i in range(n):
        for j in range(n):
            z[i,j] = d[i,j]
        for j in range(m):
            z[i,n+j] = z[n+j,i] = b[i,j]
    for i in range(m):
        for j in range(m):
            z[n+i,n+j] = -sigma[i,j]
    return z


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name,digest in {**PINS,**HELPERS}.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input/helper: '+name)
        if name in PINS:
            packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    joint,core,inertia = (packets[k] for k in (JOINT,CORE,INERTIA))
    modes = joint['remaining_modes']
    if (joint['schema'] != 'ACTUAL_C27_NINETEEN_UNIFORM_JOINT_GAIN_V1'
            or joint['source_sha256'] != '029415413d767a41fff8349386d6d484c6085bac212e3f0a3473c2d59cd32cda'
            or joint['packet_sha256'][CORE] != PINS[CORE]
            or joint['base_target_sha256'] != BASE
            or joint['packet_sha256']['joint-gamma-6765-rho32-v2-pivots.json'] != PIVOT
            or joint['grid_bits'] != BITS or joint['working_bits'] != 1024
            or joint['channels'] != 19 or joint['common_young_parameter'] != '1/32'
            or joint['exact_rationals']['g'] != ['16641320634704717921','46116860184273879040']
            or (joint['N'],joint['M'],joint['last_mode']) != (1280,32768,2048)
            or joint['verdict'] != 'stop-updated-complement-inconclusive-leading-minor'
            or joint['denominator_status'] != 'positive-definite'
            or joint['good_block']['paid_denominator_status'] != 'positive-definite'
            or (joint['new_good_block_solves'],joint['new_good_block_rhs']) != (1,19)
            or joint['fixed_joint_gain_applications'] != 1
            or joint['old_scientific_producers_or_fixtures_replayed'] != 0
            or joint['norm_weight_source_direction_parameter_cutoff_precision_scan']
            or joint['conditional_fixed_lower_target_positive']
            or joint['actual_retained_or_cofinal_positivity_claim'] or joint['actual_rh_proof']
            or len(modes) != 18 or len(joint['positive_modes']) != 1263
            or sorted(modes+joint['positive_modes']) != list(range(1281))
            or core['schema'] != 'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1'
            or core['target_packet_sha256'] != BASE or core['pivot_result_sha256'] != PIVOT
            or core['remaining_modes'] != modes or core['grid_bits'] != BITS
            or core['reused_positive_coordinate_count'] != 1263 or core['basis'] != joint['basis']
            or inertia['schema'] != 'ACTUAL_C27_SMALL_SCHUR_RANK_FLOOR_V1'
            or inertia['schur_result_sha256'] != PINS[CORE]
            or inertia['target_packet_sha256'] != BASE or inertia['pivot_result_sha256'] != PIVOT
            or inertia['remaining_modes'] != modes or inertia['reused_positive_coordinate_count'] != 1263
            or inertia['verdict'] != 'full-complement-inertia-counted'
            or inertia['negative_count_lower'] != 1):
        raise ValueError('Different SAME lower target,old core,payment or coordinates.')
    require_certificate(joint['denominator_minors'],[1]*19,'interval','sign')
    require_certificate(joint['good_block']['paid_denominator_minors'],[1]*19,'interval','sign')
    require_certificate(inertia['leading_minors'],[-1]*18,'determinant_interval','determinant_sign')
    d = matrix_in(joint['good_block']['paid_denominator'],BITS)
    b = matrix_in(joint['good_block']['lifted_beta'],BITS)
    sigma = matrix_in(core['reduced_schur_intervals'],BITS)
    if ((d.nrows(),d.ncols()),(b.nrows(),b.ncols()),(sigma.nrows(),sigma.ncols())) != (
            (19,19),(19,18),(18,18)):
        raise ValueError('Require the complete19+18 coupled border.')
    border = coupled_border(d,b,sigma)
    # Exactly ONE new library determinant, not old minors or an inverse retry.
    determinant = border.det()
    determinant_bits = 37*BITS+64
    endpoints = envelope(determinant,determinant_bits)
    sign = certified_sign(endpoints)
    verdict = {1:'conditional-fixed1281-lower-target-positive',
               -1:'fixed1281-lower-target-one-negative',
               0:'inconclusive-coupled-border-determinant'}[sign]
    for name,digest in {**PINS,**HELPERS}.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Input/helper changed during the single application: '+name)
    result = dict(schema='ACTUAL_C27_COUPLED_BORDER_DETERMINANT_V1',
        packet_sha256=PINS,helper_sha256=HELPERS,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        base_target_sha256=BASE,N=1280,M=32768,last_mode=2048,basis=joint['basis'],
        positive_modes=joint['positive_modes'],remaining_modes=modes,
        original_g=joint['exact_rationals']['g'],original_rho='1/32',
        working_bits=ctx.prec,matrix_entry_grid_bits=BITS,border_order=37,
        border_layout='[paidD19,betaV19x18;betaV_transpose,-oldSigma18]',
        old_core_inertia=[17,1,0],accepted_paid_denominator_positive=True,
        determinant_grid_bits=determinant_bits,determinant_interval=endpoints,
        determinant_sign=sign,native_determinant_calls=1,
        verdict=verdict,conditional_fixed_lower_target_positive=sign==1,
        conditional_fixed_lower_target_one_negative=sign==-1,
        new_inverse_solve_source_Gram_gain_or_scientific_producer_calls=0,
        old_determinant_or_inertia_fixture_replays=0,
        parameter_precision_basis_direction_weight_cutoff_scan=False,
        actual_retained_or_cofinal_positivity_claim=False,actual_rh_proof=False,
        runtime=dict(python=platform.python_version(),python_flint=flint.__version__,
                     machine=platform.machine(),integer_string_digit_limit=sys.get_int_max_str_digits()),
        elapsed_seconds=time.perf_counter()-started,
        conditions=joint['conditions']+[
            'The saved D=Dg+Omega and betaV enclose the SAME fully paid update; oldSigma is the SAME18 complement,not an upper threshold or its midpoint.',
            'Accepted D>0 and oldSigma inertia17+,1-,0 are reused; certificate bounds are read without rerunning their determinant/inertia computations.',
            'Schur gives detZ=detD*det(-updatedSigma)=detD*det(updatedSigma),because18 is even. The PSD gain gives updatedSigma at least17strict positives.',
            'Strictpositive detZ therefore classifies the fixed1281 lower target positive; strictnegative gives exactlyonenegative of this lower target. A crossing/nonfinite interval implies neither sign.',
            'The coupled37 border is evaluated before explicit inversion. No entrywise endpoint matrix is used as a Loewner bound; all arrays are common actual-object enclosures.',
            'One new37x37 FLINT determinant at the unchanged1024working bits;23744 is its degree-aware EXPORTgrid,not an increase in computation precision. No basis balancing or alternative evaluation/precision retry.',
            'Classical Schur,inertia and determinant methods only. No generic novelty,Lean,actual/cofinal Weil,Robin orRH conclusion.'
        ])
    digest = write_packet(out,result)
    print(verdict,digest,'seconds',result['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh.zlib output; no overwrite or retry.')
    # Exact23744-bit endpoint strings can exceed Python's default4300 digits.
    sys.set_int_max_str_digits(20000)
    ctx.prec = 1024
    consume(args.out)
