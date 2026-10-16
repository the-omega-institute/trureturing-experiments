"""Bound the common charge using the saved scalar defect and channel inverse.

Weighted Cauchy and the reverse triangle inequality contract one new
nineteen-channel scalar. Existing inverse actions and joint signed
bounds are reused; no solve, inverse or scientific supplier is run.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, N
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    packets, hashes = {}, {}
    for name in ('shared', 'joint', 'cross', 'boundary', 'base', 'defect', 'direction'):
        packets[name], hashes[name] = read_json(getattr(args, name))
    shared, joint, cross, boundary, base, defect, direction = (
        packets[name] for name in ('shared', 'joint', 'cross', 'boundary', 'base', 'defect', 'direction'))
    schemas = dict(shared='ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1',
                   joint='ACTUAL_C27_FULL_GAMMA_JOINT_SIGN_V1',
                   cross='ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1',
                   boundary='ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
                   defect='ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
                   direction='ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1')
    if any(packets[name]['schema'] != schema for name, schema in schemas.items()):
        raise ValueError('Require the accepted old-base, scalar-defect and joint-sign packets.')
    parent = shared['packet_sha256']
    if any(parent[name] != hashes[name] for name in ('cross', 'boundary', 'base', 'defect', 'direction')):
        raise ValueError('Different common row, channel inverse or scalar defect.')
    if (joint['packet_sha256']['shared'], joint['packet_sha256']['cross'],
            cross['scalar_defect_sha256'], cross['exact_direction_sha256'],
            direction['updated_target_sha256'], base['metadata']['common_young_parameter']) != (
            hashes['shared'], hashes['cross'], hashes['defect'], hashes['direction'], hashes['base'], '1/32'):
        raise ValueError('Different joint sign or common-source old nineteen-channel base.')
    expected, modes = shared['input_sha256'], shared['remaining_modes']
    if any(packets[name]['input_sha256'] != expected for name in ('joint', 'cross', 'boundary')):
        raise ValueError('Different scientific provenance.')
    if (joint['remaining_modes'], cross['remaining_modes']) != (modes, modes):
        raise ValueError('Different retained lift coordinates.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    beta = matrix_in(boundary['finite_beta'], boundary['grid_bits'])
    kernel = matrix_in(base['kernel'], base['grid_bits'])
    channels, values = len(boundary['bins']), direction['direction_dyadic_integers']
    features = kernel.nrows()-channels
    if (channels, beta.nrows(), beta.ncols(), len(values), kernel.ncols()) != (
            19, channels, N+1, N+1, features+channels) or features < 0:
        raise ValueError('Incomplete nineteen-channel inverse or exact retained direction.')
    inverse = symmetric(arb_mat([[-kernel[features+i, features+j] for j in range(channels)]
                                 for i in range(channels)]))
    v0 = arb_mat([[arb(fmpq(int(value), 1 << direction['direction_dyadic_bits']))]
                  for value in values])
    theta = arb(fmpq(int(defect['theta_integer']), 1 << defect['theta_bits']))
    if not abs(theta) > 0:
        raise ArithmeticError('Require the nonzero exact scalar-defect direction scale.')
    Ee = interval_ball(defect['weighted_residual_defect'], defect['grid_bits']).upper()
    if not Ee >= 0:
        raise ArithmeticError('Require a nonnegative accepted defect ceiling.')
    F0star = (beta*v0)/theta
    inverse_F0star = inverse*F0star
    energy = (F0star.transpose()*inverse_F0star)[0, 0]
    if not energy.is_finite() or not energy.upper() >= 0:
        raise ArithmeticError('Require a finite nonnegative center-energy upper.')
    defect_radius = (Ee/theta**2).sqrt().upper()
    center_norm_lower = max(arb(0), energy.lower()).sqrt().lower()
    center_norm_upper = energy.upper().sqrt().upper()
    channel_lower = (max(arb(0), (center_norm_lower-defect_radius).lower())**2).lower()
    channel_upper = ((center_norm_upper+defect_radius)**2).upper()
    D = interval_ball(cross['common_denominator_upper'], cross['grid_bits'])
    charge_upper = (D.upper()-channel_lower).upper()
    if not charge_upper >= 0:
        raise ArithmeticError('Inconsistent with the inherited nonnegative actual common charge.')
    qup = interval_ball(joint['scalar_bounds']['signed_inverse_metric_upper'], joint['grid_bits']).upper()
    budget = (-charge_upper-qup).lower()
    scalar = dict(channel_center_energy=energy, defect_channel_norm_upper=defect_radius,
                  channel_energy_lower=channel_lower, channel_energy_upper=channel_upper,
                  old_charge_upper=D.upper(), refined_common_charge_upper=charge_upper,
                  reused_signed_inverse_metric_upper=qup,
                  sufficient_good_metric_budget_lower=budget)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_CHANNEL_CHARGE_V1',
                  packet_sha256=hashes, input_sha256=expected, remaining_modes=modes,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  channel_center_column=matrix_out(F0star, args.output_bits),
                  channel_inverse_center_column=matrix_out(inverse_F0star, args.output_bits),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()},
                  verdict=('positive-good-metric-budget' if budget > 0
                           else 'separate-charge-budget-inconclusive'),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=joint['conditions']+[
                      'The SAME actual A_rho inverse is the saved negative channel block of the old kernel.',
                      'Fphi=F0-eta with eta=e*W_J U_J/theta and e*W_J e<=Ee for the SAME defect and original weights.',
                      'A_J=U_J*W_J U_J<=A0<A_rho, so eta A_rho^-1 eta*<=Ee/theta^2 by the operator norm, without a channel-count factor.',
                      'The actual common charge is nonnegative and Kcal<=Dstar-Fphi A_rho^-1 Fphi*.'],
                  scope='A new scalar contraction for the existing fixed common gain. '
                    'The separate charge ceiling and reused q upper supply only a sufficient budget; '
                    'a nonpositive budget neither excludes a correlated joint estimate nor establishes '
                    'actual Schur negativity. Gamma remains missing. No generic novelty, Lean, '
                    'full retained/cofinal positivity or RH conclusion.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    for name, value in scalar.items():
        print(name, value.str(16), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(shared='full-gamma-correlated-sign.json',
                             joint='full-gamma-joint-sign.json',
                             cross='full-gamma-defect-cross.json',
                             boundary='joint-gamma-6765-finite.json.zlib',
                             base='joint-gamma-6765-rho32-v2-target.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib',
                             direction='joint-gamma-6765-direction.json').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh output and the inherited640-bit grid with192 guard bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
