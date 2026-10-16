"""Assemble the finite common row on all1281 retained coordinates.

The common column,19channel map and inverse applied to the actual
channel source are cached. Reuse their full arrays and the old implicit
Schur lift; no source projection, inverse, solve or supplier is replayed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb_mat, ctx

from joint_gamma_boundary import matrix_in, matrix_out

HERE = Path(__file__).resolve().parent
BITS = 640
MODES = [0, 1, 2, 3, 4, 5, 6, 9, 10, 12, 14, 15, 18, 19, 24, 29, 33, 39]
PINS = {
    'full-gamma-common-column.json.zlib':
        '314e341189f72516aa9f251c6b873587f1bfae54d902c35f6c87c9b03130dfe7',
    'joint-gamma-6765-finite.json.zlib':
        '94aec57dacacad6bae1698a7e643514e3dace976fddbcc0008b2b7fd5aab447c',
    'full-gamma-actual-channel-projection.json.zlib':
        '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'joint-gamma-6765-rho32-v2-actions.json.zlib':
        'ff68d8350005be78daaa78a7e3f73ff7c09437dc9b1069825c5cc98b74c61962',
    'full-gamma-lift-projection.json.zlib':
        'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
}


def assemble(column, inverse_channel_column, channel_map):
    """Real application: bJ=c*-(A_rho^-1 Fphi*)* Beta."""
    n, channels = column.nrows(), inverse_channel_column.nrows()
    if (column.ncols() != 1 or inverse_channel_column.ncols() != 1
            or channel_map.nrows() != channels or channel_map.ncols() != n):
        raise ValueError('Incompatible complete column/channel dimensions.')
    return column.transpose()-inverse_channel_column.transpose()*channel_map


def restrict_to_lift(row, positive, remaining, inverse_action):
    """Use the enclosed exact X=P^-1 C, preserving its coordinate order."""
    n = row.ncols()
    if (row.nrows() != 1 or len(positive)+len(remaining) != n
            or sorted(positive+remaining) != list(range(n))
            or (inverse_action.nrows(), inverse_action.ncols()) != (len(positive), len(remaining))):
        raise ValueError('Incomplete lift partition or inverse-action dimensions.')
    good = arb_mat([[row[0, i] for i in positive]])
    exceptional = arb_mat([[row[0, i] for i in remaining]])
    return exceptional-good*inverse_action


def export(matrix):
    if any(not value.is_finite() for value in matrix.entries()):
        raise ArithmeticError('Nonfinite retained-row export; stop without retry.')
    result = matrix_out(matrix, BITS)
    if any(pair is None or len(pair) != 2 or int(pair[0]) > int(pair[1])
           for row in result for pair in row):
        raise ArithmeticError('Invalid ordered retained-row export.')
    return result


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    common, boundary, actual, actions, lift, budget = packets
    if [p['schema'] for p in packets] != [
            'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
            'ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_SCHUR_LIFT_V1',
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1']:
        raise ValueError('Different accepted packet roles.')
    if any(p['grid_bits'] != BITS for p in packets[:-1]):
        raise ValueError('Different retained-row/lift interval grid.')
    for role, name in dict(common='full-gamma-common-column.json.zlib',
                           boundary='joint-gamma-6765-finite.json.zlib',
                           projection='full-gamma-lift-projection.json.zlib').items():
        if actual['packet_sha256'][role] != PINS[name]:
            raise ValueError('Different actual-channel dependency: '+role)
    if (lift['common_column_sha256'] != PINS['full-gamma-common-column.json.zlib']
            or lift['boundary_sha256'] != PINS['joint-gamma-6765-finite.json.zlib']
            or lift['base_lift_sha256'] != PINS['joint-gamma-6765-rho32-v2-actions.json.zlib']
            or lift['base_target_sha256'] != actions['target_packet_sha256']
            or lift['base_pivots_sha256'] != actions['pivot_result_sha256']
            or lift['original_target_sha256'] != common['target_sha256']
            or common['target_sha256'] != boundary['target_packet_sha256']
            or budget['packet_sha256']['full-gamma-common-column.json.zlib'] != PINS['full-gamma-common-column.json.zlib']):
        raise ValueError('Different common source, target or exact implicit lift.')
    expected = common['metadata']['input_sha256']
    if any(p['input_sha256'] != expected for p in (boundary, actual, lift)):
        raise ValueError('Different actual scientific source inputs.')
    if (any((p['N'], p['M']) != (1280, 32768) for p in
            (boundary, actual, actions, lift, budget, common['metadata']))
            or common['metadata']['last_mode'] != 2048 or actual['first_mode'] != 2049
            or boundary['last_mode'] != 2048 or lift['last_mode'] != 2048
            or boundary['bins'] != actual['bins'] or len(boundary['bins']) != 19
            or not boundary['basis'] == lift['basis'] == common['metadata']['basis']
            or any(p['remaining_modes'] != MODES for p in (actual, actions, lift))
            or actions['positive_modes'] != lift['positive_modes']
            or budget['source_bits'] != 384 or budget['source_row'] != 0
            or budget['source_integer_sha256'] !=
            '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'):
        raise ValueError('Different source identity, phase/basis or retained/channel order.')
    positive = actions['positive_modes']
    if len(positive) != 1263 or sorted(positive+MODES) != list(range(1281)):
        raise ValueError('Incomplete retained-coordinate partition.')
    column = matrix_in(common['column'], BITS)
    inverse_channel = matrix_in(actual['channel_inverse_actual_column'], BITS)
    channel_map = matrix_in(boundary['finite_beta'], BITS)
    inverse_action = matrix_in(actions['inverse_action'], BITS)
    archived = matrix_in(actual['actual_finite_row_on_lift'], BITS)
    if ((column.nrows(), column.ncols(), inverse_channel.nrows(), inverse_channel.ncols(),
         channel_map.nrows(), channel_map.ncols(), inverse_action.nrows(), inverse_action.ncols(),
         archived.nrows(), archived.ncols()) != (1281, 1, 19, 1, 19, 1281, 1263, 18, 1, 18)):
        raise ValueError('Incomplete full retained row inputs or old lift check.')
    full = assemble(column, inverse_channel, channel_map)
    projected = restrict_to_lift(full, positive, MODES, inverse_action)
    if any(not projected[0, j].overlaps(archived[0, j]) for j in range(18)):
        raise ArithmeticError('Full finite row disagrees with the SAME archived lift projection.')
    full_export, check_export = export(full), export(projected)
    result = dict(schema='ACTUAL_C27_FULL_FINITE_ROW_V1', packet_sha256=PINS,
        source_sha256=guards[source],
        invocation_source_sha256={p.name: digest for p, digest in guards.items()},
        N=1280, M=32768, first_source_mode=2049, last_source_mode=32768,
        retained_modes=list(range(1281)), positive_modes=positive, remaining_modes=MODES,
        bins=boundary['bins'], basis=boundary['basis'], grid_bits=BITS, working_bits=1024,
        source_bits=384, source_row=0, source_integer_sha256=budget['source_integer_sha256'],
        input_sha256=expected, finite_common_row=full_export,
        same_lift_projection_check=check_export,
        archived_lift_row_source='actual.actual_finite_row_on_lift',
        verdict='conditional-complete-finite-row-assembled-and-lift-consistent',
        new_channel_row_products=1, new_lift_consistency_products=1,
        lift_consistency_coordinates=18, new_retained_coefficients=1281,
        old_producers_replayed=0, new_source_or_adjoint_projections=0,
        new_inverse_or_solve_calls=0, parameter_source_precision_scan=False,
        conditions=actual['conditions']+lift['conditions']+[
            'The real unphased common column, finite19channel map and cached A_rho^-1 Fphi* describe the SAME complete finite source, weights and retained residual map.',
            'The saved1263x18 inverse_action encloses the SAME exact X=P^-1 C and is used only for V=[-X;I], not as an inverse for arbitrary new good-block rows.',
            'Projection overlap is a numerical consistency check; source identity and the exact distributive formula establish the common mathematical row.',
            'All1281 retained coordinates of bJ are assembled, but bJ omits the actual far mixed row. No full external-source B, unweighted r row or good-block inverse metric is supplied.'],
        scope='Finite-source common row on all retained coordinates, for later same-source '
            'good-block metric work. The old18lift row is reused as a check. No actual '
            'sufficient score, retained/cofinal positivity, maximal capacity, Robin-tail, '
            'RH, generic novelty or Lean result.', elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed: '+name)
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed: '+path.name)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':')).encode(), 9))
    print('complete finite retained row', len(full_export[0]), flush=True)
    print('same-lift consistency coordinates', len(check_export[0]), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh .zlib output; accepted archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
