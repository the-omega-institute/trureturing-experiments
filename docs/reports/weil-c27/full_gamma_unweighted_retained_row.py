"""Supply the missing complete unweighted row r=w*R on retained modes.

Use the SAME exact finite source and existing residual generators. This is
one new unweighted Cauchy product, not the archived weighted common column.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from gamma_channel_resolution import cauchy_projection
from full_gamma_finite_retained_row import export, restrict_to_lift
from joint_gamma_boundary import matrix_in, write_packet

HERE = Path(__file__).resolve().parent
BITS, ORDER = 640, 512
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
FILES = (
    'full-gamma-dual-32768-finite.json.zlib',
    'joint-residual-generators-2048.json.zlib',
    'structured-joint-2048-target.json.zlib',
    'joint-gamma-6765-rho32-v2-actions.json.zlib',
    'full-gamma-projection-upper.json.zlib',
)


def unweighted_coefficients(integers, bits, first, last, retained_last):
    """Place the original exact source; no W, M, kappa or source rescaling."""
    if (type(bits) is not int or bits <= 0
            or not retained_last < first <= last
            or len(integers) != last-first+1
            or any(type(x) is not str for x in integers)):
        raise ValueError('Incomplete exact finite source support.')
    coefficients = arb_mat(last-retained_last, 1)
    for m, value in enumerate(integers, first):
        coefficients[m-retained_last-1, 0] = arb(fmpq(int(value), 1 << bits))
    return coefficients


def consume(out):
    started = time.perf_counter()
    values, pins = {}, {}
    for name in FILES:
        values[name], pins[name] = read_json(HERE/name)
    finite, generators, target, actions, projection = [values[n] for n in FILES]
    if [v['schema'] for v in values.values()] != [
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', 'ACTUAL_C27_RESIDUAL_GENERATOR_V1',
            'ACTUAL_C27_STRUCTURED_TARGET_V1', 'ACTUAL_C27_SCHUR_LIFT_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1']:
        raise ValueError('Different source/generator/lift/projection roles.')
    metadata = finite['metadata']
    if (metadata != dict(target_sha256=pins[FILES[2]],
            support_sha256=generators['support_action_sha256'],
            generators_sha256=pins[FILES[1]], direction_sha256=metadata['direction_sha256'],
            input_sha256=target['metadata']['input_sha256'], N=1280, M=32768,
            last_mode=2048, basis=target['metadata']['basis'])
            or generators['metadata']['target_packet_sha256'] != pins[FILES[2]]
            or generators['metadata']['input_sha256'] != metadata['input_sha256']
            or (generators['metadata']['N'],generators['metadata']['M'],
                generators['metadata']['last_mode'],generators['metadata']['basis'])
                != (1280,32768,2048,metadata['basis'])
            or actions['N'] != 1280 or actions['M'] != 32768
            or projection['packet_sha256']['old'] != pins[FILES[2]]
            or projection['packet_sha256']['generators'] != pins[FILES[1]]
            or projection['packet_sha256']['finite'] != pins[FILES[0]]
            or projection['packet_sha256']['actions'] != pins[FILES[3]]
            or projection['packet_sha256']['base'] != actions['target_packet_sha256']
            or projection['remaining_modes'] != actions['remaining_modes']
            or projection['exact_source_bits'][0] != 384
            or projection['exact_source_integer_sha256'][0] != SOURCE_ID
            or (finite['dual_first'], finite['dual_last'], finite['dual_bits']) != (2049,32768,384)
            or any(v['grid_bits'] != BITS for v in values.values())
            or metadata['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different original unphased finite source or exact lift.')
    coefficients = unweighted_coefficients(finite['dual_integers'],384,2049,32768,1280)
    GO, GP = matrix_in(generators['GO'],BITS), matrix_in([r[:2] for r in target['generators']],BITS)
    if (GO.nrows(),GO.ncols(),GP.nrows(),GP.ncols()) != (30720,2,1281,2):
        raise ValueError('Incomplete residual-generator map.')
    row, geometry = cauchy_projection(coefficients, lambda m,j:GO[m-2049,j], GP,
        list(range(2049,32769)),list(range(1281)),interval_ball(target['kappa'],BITS),1280,ORDER)
    lifted = restrict_to_lift(row,actions['positive_modes'],actions['remaining_modes'],
        matrix_in(actions['inverse_action'],BITS))
    saved = matrix_in([projection['matrices']['source_on_lift'][0]],BITS)
    if (saved.nrows(),saved.ncols()) != (1,18) or any(
            not lifted[0,j].overlaps(saved[0,j]) for j in range(18)):
        raise ArithmeticError('New complete unweighted row disagrees with its SAME archived lift.')
    result = dict(schema='ACTUAL_C27_COMPLETE_UNWEIGHTED_RETAINED_ROW_V1',
        packet_sha256=pins, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID,source_bits=384,source_row=0,
        N=1280,M=32768,first_source_mode=2049,last_source_mode=32768,
        retained_modes=list(range(1281)),positive_modes=actions['positive_modes'],
        remaining_modes=actions['remaining_modes'],basis=metadata['basis'],
        input_sha256=metadata['input_sha256'],original_target_sha256=pins[FILES[2]],
        base_target_sha256=actions['target_packet_sha256'],
        base_lift_sha256=pins[FILES[3]],grid_bits=BITS,working_bits=ctx.prec,
        complete_unweighted_row=export(row),same_lift_projection_check=export(lifted),
        projection_geometry=geometry,new_unweighted_source_projections=1,
        new_retained_coefficients=1281,new_inverse_or_solve_calls=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        complete_for_original_finite_source=True,complete_b_row_supplied=False,
        actual_retained_sign_claim=False,actual_full_border_sign_claim=False,actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,conditions=[
            'The original exact384-bit row0 w has finite support2049..32768. r=w*R therefore has no extra w-source tail. This does not remove the far mixed part in b=u*Wbar R-w*R.',
            'Existing GO/GP,kappa and their geometric absolute tail represent the SAME full retained residual map and unphased basis. No W or Mw enters the new coefficients; kappa is applied only inside the existing projection.',
            'The exact old lift identifies rV with the saved source_on_lift. Overlap checks implementation, not source equality or positivity.',
            'All inherited actual form, residual map, domain, source and exact trial premises remain. No generic novelty, Lean proof, actual retained/full-border/cofinal positivity, Robin or RH result.',
        ])
    digest=write_packet(out,result)
    print('complete unweighted retained row',digest,'seconds',result['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; no overwrite or retry.')
    ctx.prec=1024
    consume(args.out)
