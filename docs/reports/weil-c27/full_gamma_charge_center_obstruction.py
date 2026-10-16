"""Test the zero-error point in the stated common-charge relaxation.

Consume accepted center scalars only, using exact integer endpoint
arithmetic. No inverse, contraction, supplier or prior consumer runs.
The witness belongs to an error-envelope relaxation; it is not claimed
to come from the actual common source.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib


HERE = Path(__file__).resolve().parent


def read(path):
    data = path.read_bytes()
    return json.loads(zlib.decompress(data) if path.suffix == '.zlib' else data), hashlib.sha256(data).hexdigest()


def endpoints(packet, field, section='scalar_bounds'):
    values = (packet[section] if section is not None else packet)[field]
    if len(values) != 2 or any(not isinstance(value, str) for value in values):
        raise ValueError('Require two exact integer endpoint strings: '+field)
    lower, upper = map(int, values)
    if lower > upper:
        raise ValueError('Reversed exact endpoint interval: '+field)
    return lower, upper


def consume(args):
    started = time.monotonic()
    charge, charge_hash = read(args.charge)
    joint, joint_hash = read(args.joint)
    common, common_hash = read(args.common)
    defect, defect_hash = read(args.defect)
    cross, cross_hash = read(args.cross)
    if (charge['schema'], joint['schema'], charge['grid_bits'], joint['grid_bits'],
            charge['packet_sha256']['joint'], charge['packet_sha256']['shared'],
            charge['packet_sha256']['cross']) != (
            'ACTUAL_C27_FULL_GAMMA_CHANNEL_CHARGE_V1',
            'ACTUAL_C27_FULL_GAMMA_JOINT_SIGN_V1', 640, 640, joint_hash,
            joint['packet_sha256']['shared'], joint['packet_sha256']['cross']):
        raise ValueError('Require accepted common-source charge and signed-center packets.')
    if (charge['input_sha256'], charge['remaining_modes']) != (
            joint['input_sha256'], joint['remaining_modes']):
        raise ValueError('Different scientific source or retained lift.')
    if (common['schema'], defect['schema'], cross['schema'],
            common['grid_bits'], defect['grid_bits'],
            charge['packet_sha256']['cross'], charge['packet_sha256']['defect'],
            cross['common_column_sha256'], cross['scalar_defect_sha256'],
            common['finite_defect_sha256'], common['metadata']) != (
            'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1', 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
            'ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1', 640, 640,
            cross_hash, defect_hash, common_hash, defect_hash, defect_hash, defect['metadata']):
        raise ValueError('Different common finite energy, dual denominator or scalar defect.')
    D = endpoints(charge, 'old_charge_upper')
    E0 = endpoints(charge, 'channel_center_energy')
    qc = endpoints(joint, 'center_signed_energy')
    center_charge = (D[0]-E0[1], D[1]-E0[0])
    center_sign = (center_charge[0]+qc[0], center_charge[1]+qc[1])
    middle = endpoints(common, 'middle_action_energy', 'bounds')
    denominator = endpoints(defect, 'denominator', None)
    finite_floor = (middle[0]-denominator[1], middle[1]-denominator[0])
    floor_charge = (finite_floor[0]-E0[1], finite_floor[1]-E0[0])
    floor_sign = (floor_charge[0]+qc[0], floor_charge[1]+qc[1])
    if not (E0[0] >= 0 and center_charge[0] > 0 and qc[1] < 0 and center_sign[0] > 0
            and floor_charge[0] > 0 and floor_sign[0] > 0):
        raise ArithmeticError('The cached center does not certify this positive obstruction.')
    intervals = dict(old_charge_upper=D, channel_center_energy=E0,
                     signed_center_energy=qc, relaxed_zero_error_charge=center_charge,
                     relaxed_zero_error_joint_scalar=center_sign,
                     middle_action_energy=middle, dual_denominator=denominator,
                     finite_common_charge_floor=finite_floor,
                     relaxed_zero_error_charge_floor=floor_charge,
                     relaxed_zero_error_joint_floor=floor_sign)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_CHARGE_CENTER_OBSTRUCTION_V1',
                  packet_sha256=dict(charge=charge_hash, joint=joint_hash,
                                     common=common_hash, defect=defect_hash, cross=cross_hash),
                  input_sha256=charge['input_sha256'], remaining_modes=charge['remaining_modes'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  grid_bits=640,
                  scalar_intervals={name: list(map(str, values)) for name, values in intervals.items()},
                  relaxed_point=dict(finite_channel_error='zero row in19channels',
                                     far_lift_error='zero row in18retained coordinates',
                                     good_metric='any Gamma>=0',
                                     far_action_energy='any s>=0, including the relaxed value zero'),
                  verdict='zero-error-point-obstructs-uniform-relaxed-upper-sign',
                  elapsed_seconds=time.monotonic()-started,
                  conditions=charge['conditions']+[
                      'The relaxation constrains eta and t only by the accepted PSD error upper envelopes; it imposes no positive lower size or further same-source realization condition.',
                      'The tested denominator is Dstar-(F0-eta) A_rho^-1 (F0-eta)*, not the exact actual common charge.',
                      'The same-source charge identity is Kcal=middle_action_energy-h+(1+rho)*s-(F0-eta) A_rho^-1 (F0-eta)*, with s=||u_far||_W^2>=0 and rho=1/32.',
                      'The accepted old core has one negative direction and a positive good block; the same inherited rank-one criterion is used only with the positive relaxed denominator.'],
                  scope='At eta=t=0, both error upper envelopes hold and the positive relaxed '
                    'denominator plus q is strictly positive even before adding Gamma>=0. '
                    'Even replacing the far action energy by its ideal relaxed value zero '
                    'leaves a strictly positive lower bound for the exact-charge scalar '
                    'Kcal+Gamma+q at eta=t=0. Thus a uniform negative upper on this stated '
                    'relaxation cannot follow from these envelopes or a sharper far-energy upper alone. '
                    'This is not an actual-source realization, an exclusion of sharper charge '
                    'data or other lower constructions, actual Schur negativity, or an RH result. '
                    'Existing feasibility reasoning and exact cached arithmetic; no generic novelty or Lean claim.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    for name, (lower, upper) in intervals.items():
        print(name, f'({lower/(1 << 640):.14f}, {upper/(1 << 640):.14f})')
    print(result['verdict'], 'seconds', result['elapsed_seconds'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--charge', type=Path, default=HERE/'full-gamma-channel-charge.json')
    parser.add_argument('--joint', type=Path, default=HERE/'full-gamma-joint-sign.json')
    parser.add_argument('--common', type=Path, default=HERE/'full-gamma-common-column.json.zlib')
    parser.add_argument('--defect', type=Path, default=HERE/'full-gamma-dual-32768-defect.json.zlib')
    parser.add_argument('--cross', type=Path, default=HERE/'full-gamma-defect-cross.json')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    consume(args)


if __name__ == '__main__':
    main()
