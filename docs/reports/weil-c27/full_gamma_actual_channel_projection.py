"""Project the cached common effective action onto the original19channels.

Supply the missing Fphi=u_J*W_J U_J and finite common row. Evaluate only
the original channel shapes needed by this new functional, one row at a
time; their full values were not archived. The effective action, old
channel Gram/map, exact lift and all inverse actions are reused.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import acb, arb, arb_mat, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_inverse_residual_obstruction import HERE, M, N
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def produce(args):
    started = time.monotonic()
    packets, hashes = {}, {}
    for name in ('charge', 'shared', 'projection', 'defect', 'common', 'pilot', 'boundary', 'base', 'cross'):
        packets[name], hashes[name] = read_json(getattr(args, name))
    charge, shared, projection, defect, common, pilot, boundary, base, cross = (
        packets[name] for name in ('charge', 'shared', 'projection', 'defect', 'common', 'pilot', 'boundary', 'base', 'cross'))
    schemas = dict(charge='ACTUAL_C27_FULL_GAMMA_CHANNEL_CHARGE_V1',
                   shared='ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1',
                   projection='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                   defect='ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
                   common='ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
                   boundary='ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
                   cross='ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1')
    if any(packets[name]['schema'] != value for name, value in schemas.items()):
        raise ValueError('Require the same accepted original19channel common-source packets.')
    for name in ('shared', 'defect', 'boundary', 'base', 'cross'):
        if charge['packet_sha256'][name] != hashes[name]:
            raise ValueError('Different accepted channel-charge provenance: '+name)
    if any(shared['packet_sha256'][name] != hashes[name] for name in (
            'projection', 'common', 'boundary', 'base', 'defect', 'cross')):
        raise ValueError('Different common source or old base/lift.')
    if (common['finite_defect_sha256'], common['pilot_sha256'], common['metadata'],
            common['metadata']['N'], common['metadata']['M'], common['metadata']['last_mode'],
            projection['boundary_sha256'], projection['base_target_sha256'],
            base['metadata']['common_young_parameter'],
            boundary['boundary_source_sha256']) != (
            hashes['defect'], hashes['pilot'], defect['metadata'], N, M, 2048,
            hashes['boundary'], hashes['base'], '1/32',
            hashlib.sha256((HERE/'gamma_boundary_update.py').read_bytes()).hexdigest()):
        raise ValueError('Different actual action, original weights or channel-definition source.')
    expected, modes = charge['input_sha256'], charge['remaining_modes']
    if any(packet['input_sha256'] != expected for packet in (shared, projection, boundary, cross)):
        raise ValueError('Different scientific provenance.')
    if any(packet['remaining_modes'] != modes for packet in (shared, projection, cross)):
        raise ValueError('Different full retained lift.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    bins = boundary['bins']
    if len(bins) != 19 or bins[0] != [0, 1] or bins[-1][1] != 6765 or any(
            not isinstance(a, int) or not isinstance(b, int) or a >= b
            for a, b in bins):
        raise ValueError('Require the unchanged original19Fibonacci channel bins.')
    mw = matrix_in(defect['effective_dual_action'], defect['grid_bits'])
    weights = [arb(fmpq(int(value), 1 << pilot['middle_weight_bits']))
               for value in pilot['middle_weight_numerators']]
    if (mw.nrows(), mw.ncols(), len(weights)) != (M-2048, 1, M-N) or not all(
            weight > 0 for weight in weights):
        raise ValueError('Incomplete cached action or original positive weights.')
    arithmetic = Arithmetic(27, args.output_bits)
    factors = [((1-(-arb(fmpq(4*a+1, 2))*arithmetic.d).exp())
                /((b-a)*arithmetic.d)).sqrt() for a, b in bins]
    endpoints = sorted({value for pair in bins for value in pair})
    Fphi = arb_mat(1, len(bins))
    for i, m in enumerate(range(2049, M+1)):
        im = arithmetic.pi*m/arithmetic.d
        values = {value: acb(arb(value)+arb(fmpq(1, 4)), im).digamma().real
                  for value in endpoints}
        weighted = (-1)**m*weights[m-N-1]*mw[i, 0]
        for j, (a, b) in enumerate(bins):
            Fphi[0, j] += weighted*factors[j]*(values[b]-values[a])
        if (i+1) % 4096 == 0:
            print('new actual19channel projection through', m,
                  'seconds', time.monotonic()-started, flush=True)
    kernel = matrix_in(base['kernel'], base['grid_bits'])
    rank = len(projection['far_features_on_lift'])
    if (kernel.nrows(), kernel.ncols()) != (rank+19, rank+19):
        raise ValueError('Incomplete saved original19channel inverse.')
    inverse = symmetric(arb_mat([[-kernel[rank+i, rank+j] for j in range(19)] for i in range(19)]))
    F0star = matrix_in(charge['channel_center_column'], charge['grid_bits'])
    if (F0star.nrows(), F0star.ncols()) != (19, 1):
        raise ValueError('Require the accepted cached scalar-defect channel center.')
    eta = F0star.transpose()-Fphi
    inverse_phi = inverse*Fphi.transpose()
    energy = (Fphi*inverse_phi)[0, 0]
    eta_energy = (eta*inverse*eta.transpose())[0, 0]
    radius = interval_ball(charge['scalar_bounds']['defect_channel_norm_upper'], charge['grid_bits'])
    if not eta_energy.is_finite() or eta_energy.lower() > radius.upper()**2:
        raise ArithmeticError('Actual projected defect contradicts the accepted channel-error envelope.')
    c = matrix_in(projection['common_column_on_lift'], projection['grid_bits'])
    betaV = matrix_in(projection['channel_map_on_lift'], projection['grid_bits'])
    if (c.nrows(), c.ncols(), betaV.nrows(), betaV.ncols()) != (1, len(modes), 19, len(modes)):
        raise ValueError('Incomplete original finite column or channel map on the same lift.')
    finite_row = c-inverse_phi.transpose()*betaV
    Z = inverse*betaV
    old_center = arb_mat([[interval_ball(row['finite_center'], cross['grid_bits'])
                           for row in cross['rows']]])
    transport_check = old_center+eta*Z
    if any(not finite_row[0, j].overlaps(transport_check[0, j]) for j in range(len(modes))):
        raise ArithmeticError('Direct channel projection disagrees with the cached same-source transport.')
    Dstar = interval_ball(charge['scalar_bounds']['old_charge_upper'], charge['grid_bits'])
    charge_upper = (Dstar.upper()-energy.lower()).upper()
    scalar = dict(actual_channel_energy=energy, actual_projected_defect_energy=eta_energy,
                  inherited_projected_defect_energy_upper=(radius.upper()**2).upper(),
                  refined_common_charge_upper=charge_upper)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
                  packet_sha256=hashes, input_sha256=expected, remaining_modes=modes,
                  N=N, M=M, first_mode=2049, bins=bins,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  actual_channel_row=matrix_out(Fphi, args.output_bits),
                  actual_channel_error_row=matrix_out(eta, args.output_bits),
                  channel_inverse_actual_column=matrix_out(inverse_phi, args.output_bits),
                  actual_finite_row_on_lift=matrix_out(finite_row, args.output_bits),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()},
                  transport_coordinates_checked=len(modes), elapsed_seconds=time.monotonic()-started,
                  conditions=charge['conditions']+[
                      'The original19channel shapes are evaluated by the unchanged boundary-definition formula, with no pole/contrast rescaling.',
                      'The same exact effective action u_J=Mw and original weights are enclosed by the accepted cached arrays.',
                      'Only Fphi=u_J*W_J U_J is newly projected; old Mw, residuals, channel Grams/maps and inverse actions are not regenerated.'],
                  scope='Actual finite-source channel information for the fixed common gain. '
                    'The finite common row now uses the actual Fphi; the far mixed row and '
                    'good-block inverse metric remain missing. No sign certificate, generic '
                    'novelty, Lean, full retained/cofinal positivity, Robin-tail or RH claim.')
    write_packet(args.out, result)
    for name, value in scalar.items():
        print(name, value.str(20), flush=True)
    print('actual finite-channel source enclosed; seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(charge='full-gamma-channel-charge.json',
                             shared='full-gamma-correlated-sign.json',
                             projection='full-gamma-lift-projection.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib',
                             common='full-gamma-common-column.json.zlib',
                             pilot='inverse-residual-pilot.json.zlib',
                             boundary='joint-gamma-6765-finite.json.zlib',
                             base='joint-gamma-6765-rho32-v2-target.json.zlib',
                             cross='full-gamma-defect-cross.json').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh .zlib output and inherited640-bit grid with192guardbits.')
    ctx.prec = args.precision
    produce(args)


if __name__ == '__main__':
    main()
