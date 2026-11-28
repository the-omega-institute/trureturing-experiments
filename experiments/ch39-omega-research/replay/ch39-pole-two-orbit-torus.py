"""Certify two-orbit triple erasures on the common two-clock torus."""

from argparse import ArgumentParser
from collections import Counter
from ctypes import CDLL, c_size_t
from itertools import combinations, permutations, product
import json
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
import numpy as np
from numpy.ctypeslib import ndpointer


SCALE = 2**30
CELLS = 256
LEFT_SCALE = 256
LEFT_CAP = 100_000
DEFECT_SHIFT = 2**15
TRANSFER_UNIT_STEP = 6 * 22
DEPTH = 13
FLOOR = '5e-7'
HERE = Path(__file__).resolve().parent


def cubic_controls(nodes):
    """Convert equally spaced interpolation nodes to Bernstein controls."""
    def axis(data, direction):
        lo = np.take(data, np.arange(0, data.shape[direction], 3), axis=direction)
        one = np.take(data, np.arange(1, data.shape[direction], 3), axis=direction)
        two = np.take(data, np.arange(2, data.shape[direction], 3), axis=direction)
        hi = np.roll(lo, -1, axis=direction)
        first = (18 * one - 9 * two - 5 * lo + 2 * hi) // 6
        second = (18 * two - 9 * one + 2 * lo - 5 * hi) // 6
        return np.stack((lo, first, second, hi), axis=-1)

    return axis(axis(nodes, 0), 1)


def split_cubic(data, axis):
    a, b, c, d = np.moveaxis(data, axis, 0)
    ab, bc, cd = (a + b) // 2, (b + c) // 2, (c + d) // 2
    ac, bd = (ab + bc) // 2, (bc + cd) // 2
    middle = (ac + bd) // 2
    return [np.stack(group, axis=axis) for group in ((a, ab, ac, middle), (middle, bd, cd, d))]


def physical_line(nets, weights, line, scale):
    cells = nets[0].shape[1]
    controls = np.empty((cells, 4, 4, 18, 15, 2), dtype=np.int32)
    for k, (net, weight) in enumerate(zip(nets, weights)):
        block = net[line].astype(np.int64)
        re = block[..., 0][..., None, :]
        im = block[..., 1][..., None, :]
        wr = weight[:, 0, None]
        wi = weight[:, 1, None]
        physical = controls[..., 5*k:5*(k+1), :].reshape(cells, 4, 4, 6, 3, 5, 2)
        real = (re * wr - im * wi + scale // 2) // scale
        imag = (re * wi + im * wr + scale // 2) // scale
        # This bound also keeps the sums in int32 during subdivision.
        assert -scale < real.min() and real.max() < scale
        assert -scale < imag.min() and imag.max() < scale
        physical[..., 0] = real
        physical[..., 1] = imag
    return controls


def verify_nets(balls, scale, threshold):
    weights = np.array([1, 3, 3, 1], dtype=np.float64) / 8
    center_parts = np.einsum('i,j,nij...->n...', weights, weights, balls) / scale
    center = center_parts[..., 0] + 1j * center_parts[..., 1]
    proposed = np.linalg.pinv(center)
    assert np.isfinite(proposed).all()
    parts = np.stack((proposed.real, proposed.imag), axis=-1)
    left = np.rint(LEFT_SCALE * np.clip(parts, -LEFT_CAP, LEFT_CAP)).astype(np.int64)
    assert np.abs(left).max() <= LEFT_SCALE * LEFT_CAP
    assert -scale < balls.min() and balls.max() < scale
    unit = scale * LEFT_SCALE
    max_left_square = 15 * 16 * 2 * (LEFT_SCALE * LEFT_CAP) ** 2
    max_product = 32 * LEFT_SCALE * LEFT_CAP * 3 * scale + unit
    assert max(max_left_square, max_product) < 2 ** 63
    norms = np.fromiter((isqrt(int(x)) + 1 for x in (left**2).sum(axis=(1,2,3))), count=len(left), dtype=np.int64)
    lr, li = left[..., 0], left[..., 1]
    eye = np.eye(15, dtype=np.int64) * unit
    max_defect = np.zeros(len(balls), dtype=np.int64)
    for j, k in product(range(4), repeat=2):
        r, i = balls[:, j, k, ..., 0], balls[:, j, k, ..., 1]
        re = np.minimum(abs(lr @ r - li @ i - eye), unit)
        im = np.minimum(abs(lr @ i + li @ r), unit)
        re = (re + DEFECT_SHIFT - 1)//DEFECT_SHIFT
        im = (im + DEFECT_SHIFT - 1)//DEFECT_SHIFT
        assert 450 * ((unit + DEFECT_SHIFT - 1)//DEFECT_SHIFT)**2 < 2**63
        maximum = (re**2 + im**2).sum(axis=(1,2))
        max_defect = np.maximum(max_defect, maximum)
    defect = np.fromiter((isqrt(int(x)) + 1 for x in max_defect), count=len(left), dtype=np.int64)
    left_margin = unit - defect * DEFECT_SHIFT
    assert (isqrt(max_left_square) + 1) * threshold < 2**63
    return (left_margin > 0) & (norms * threshold < left_margin)


def control_net(k, cells, scale, basis, full, pump, shell, storage, node_scale):
    # The seventh coordinate in character zero is the pole. Retain it in the
    # same interpolation and eigenline error estimates as the ordinary rows.
    dim = basis.ncols()
    assert node_scale >= scale and node_scale % scale == 0 and 64 * node_scale < 2**63
    divisor = node_scale // scale
    full['quantized'].__globals__['SCALE'] = node_scale
    net = np.lib.format.open_memmap(storage / f'character-{k}.npy', mode='w+',
                                   dtype=np.int32, shape=(cells, cells, 4, 4, dim, 5, 2))
    squared_fourth_derivatives = [arb(0), arb(0)]
    residuals = []
    degree = full['DEGREE']
    for col, charge in enumerate(full['CHARGES']):
        base, links = pump['clock_matrices'](shell, charge)
        adjoint = pump['adjoint'](basis)
        operators = [adjoint * matrix * basis for matrix in (base, *links)]
        trial, _, residual = full['verify_line'](pump, operators, charge)
        fourth_derivatives = [arb(0), arb(0)]
        for i, j in np.ndindex(trial.shape[:2]):
            norm = sum(abs(z)**2 for z in trial[i, j]).sqrt()
            for axis, index in enumerate((i, j)):
                fourth_derivatives[axis] += norm * (charge * (index - degree))**4
        for axis, bound in enumerate(fourth_derivatives):
            squared_fourth_derivatives[axis] += bound**2
        residuals.append(residual)
        for row in range(dim):
            grid = full['sample'](trial[:, :, row], 3 * cells, charge)
            for part in range(2):
                values = [[z.imag if part else z.real for z in line] for line in grid]
                controls = cubic_controls(full['quantized'](values))
                controls = (controls + divisor // 2) // divisor
                assert -scale < controls.min() and controls.max() < scale
                net[:, :, :, :, row, col, part] = controls
    net.flush()
    step = 2 * arb.pi() / cells
    lo, hi = sorted(bound.sqrt().upper() for bound in squared_fourth_derivatives)
    interpolation = (hi + arb(7) / 4 * lo) * step**4 / (81 * 24)
    residual = max(residuals) + 6 * arb('1e-8') * arb('1.5')
    relative = residual / arb('.4')
    assert relative < arb('.5')
    line_error = residual / (1 - relative)
    print('character', k, 'interpolation', interpolation.str(8), 'line', line_error.str(8), flush=True)
    return net, interpolation + arb(5).sqrt() * line_error


def pole_cases(sites, orbit_count):
    orbits = [sites[start:start + 3] for start in range(0, len(sites), 3)]
    assert len(orbits) == 6 and len(set(sites)) == 18
    rotation = {}
    for orbit in orbits:
        for j, (x, y, z) in enumerate(orbit):
            image = orbit[(j + 1) % 3]
            assert image == (y, z, x)
            rotation[orbit[j]] = image
    covered = set()
    cases = []
    if orbit_count == 2:
        masks = ((o[1], o[2], site) for o, p in permutations(orbits, 2) for site in p)
    else:
        masks = ((a[0], b[j], c[k]) for a, b, c in combinations(orbits, 3)
                 for j, k in product(range(3), repeat=2))
    for names in masks:
        rotated = {frozenset(names),
                   frozenset(rotation[s] for s in names),
                   frozenset(rotation[rotation[s]] for s in names)}
        assert len(rotated) == 3 and not covered.intersection(rotated)
        covered.update(rotated)
        keep = [i for i, s in enumerate(sites) if s not in names] + [18]
        assert len(keep) == 16
        cases.append((names, keep))
    count = {2: 90, 3: 180}[orbit_count]
    assert len(cases) == count and len(covered) == 3 * count
    return cases


def exact_defects(kernel, left, balls):
    assert left.shape == (len(balls), 15, 16, 2)
    assert balls.shape[1:] == (4, 4, 16, 15, 2)
    assert left.dtype == balls.dtype == np.int32
    assert max(int(left.max()), -int(left.min())) <= LEFT_SCALE * LEFT_CAP < 2**25
    assert max(int(balls.max()), -int(balls.min())) < SCALE
    defects = np.empty(len(balls), dtype=np.int64)
    kernel(np.ascontiguousarray(left), np.ascontiguousarray(balls), defects, len(balls))
    return defects


def native_check(kernel, balls, scale, threshold):
    assert scale == SCALE and 0 < threshold < scale
    midpoint = np.array([1, 3, 3, 1], dtype=np.float64) / 8
    parts = np.einsum("i,j,nij...->n...", midpoint, midpoint, balls) / scale
    approximate = np.linalg.pinv(parts[..., 0] + 1j * parts[..., 1])
    assert np.isfinite(approximate).all()
    pairs = np.stack((approximate.real, approximate.imag), axis=-1)
    left = np.rint(LEFT_SCALE * np.clip(pairs, -LEFT_CAP, LEFT_CAP)).astype(np.int32)
    defects_squared = exact_defects(kernel, left, balls)
    norm_squared = (left.astype(np.int64)**2).sum(axis=(1, 2, 3))
    norm = np.fromiter((isqrt(int(value))+1 for value in norm_squared), dtype=np.int64)
    defect = np.fromiter((isqrt(int(value))+1 for value in defects_squared), dtype=np.int64)
    margin = SCALE * LEFT_SCALE - defect * DEFECT_SHIFT
    assert (isqrt(15*16*2*(LEFT_SCALE*LEFT_CAP)**2)+1)*threshold < 2**63
    return (margin > 0) & (norm*threshold < margin)


def load_native_kernel(path):
    library = CDLL(str(path.resolve()))
    kernel = library.max_squared_defects
    kernel.argtypes = [
        ndpointer(dtype=np.int32, flags=("C_CONTIGUOUS", "ALIGNED")),
        ndpointer(dtype=np.int32, flags=("C_CONTIGUOUS", "ALIGNED")),
        ndpointer(dtype=np.int64, flags=("C_CONTIGUOUS", "ALIGNED")),
        c_size_t,
    ]
    kernel.restype = None
    return kernel


def generate(model, storage):
    full = run_path(str(HERE / 'ch39-five-charge-full-torus.py'))
    pump = run_path(str(model / 'ch39-native-light-projector-pump.py'))
    shell = run_path(str(model / 'ch39-three-ray-chiral-boundary.py'))
    pump['verify'](model)
    ctx.prec = 160
    full['DEGREE'] = 14
    full['verify_line'].__globals__['DEGREE'] = full['DEGREE']
    storage.mkdir(parents=False, exist_ok=False)
    errors = []
    for character, frame in enumerate(pump['rotation_bases'](shell)):
        _, error = control_net(character, CELLS, SCALE, frame, full, pump, shell, storage,
                                        node_scale=2**42)
        errors.append(error.upper())
    assert len(errors) == 3
    # Two tensor conversions amplify the high-scale nodal error by the
    # cubic l1 bound 17/3 and add one unit of truncation per axis.
    high_scale_error = arb(17)/3 * (arb(17)/3 + 1) + 1
    coarse_error = high_scale_error / 2**12 + arb('.5')
    assert coarse_error < arb('.51')
    # Each physical row uses one character coefficient. Rounding the
    # coefficient contributes at most one unit and rounding the product
    # contributes at most one half. Projection onto retained rows is contractive.
    assert (arb(2)/3).sqrt() * coarse_error + arb('1.5') < arb(2)
    assert arb(480).sqrt() < 22
    bound = max(errors) + arb(22 * 3) / SCALE
    # Store a decimal that is verified above the complete depth-zero bound.
    numerator = int(float(bound.upper()) * 10**18) + 2
    decimal = f'{numerator // 10**18}.{numerator % 10**18:018}'
    assert arb(decimal) > bound
    threshold = int(float(bound.upper()) * SCALE) + 2
    assert arb(threshold) / SCALE > arb(decimal)
    settings = {'transfer': decimal, 'threshold': threshold, 'accounted_splits': 0}
    (storage / 'transfer.json').write_text(json.dumps(settings) + '\n')
    print(json.dumps(settings), flush=True)


def verify_reflection(shell):
    sites = shell['geo']['sites']
    pole = shell['soft']['pole']
    tuned = set(shell['soft']['star']['orbit'])
    shift = [sites.index((y, z, x)) for x, y, z in sites]
    flip = [sites.index((y, x, z)) for x, y, z in sites]
    assert all(flip[flip[i]] == i for i in range(19))
    assert all(flip[shift[flip[i]]] == shift[shift[i]] for i in range(19))
    assert flip[pole] == pole
    reflected_tuned = {flip[i] for i in tuned}
    assert len(reflected_tuned) == 3 and {shift[i] for i in reflected_tuned} == reflected_tuned
    geo = shell['geo']
    edges = {(i, j): distance for i, j, distance, _ in geo['edges']}
    pole_site = sites[pole]
    dot, det = geo['dot'], geo['det']

    def phase_numerator(i, j):
        x, y = sites[i], sites[j]
        return geo['R'] * (geo['R'] + dot(pole_site, x) + dot(x, y) + dot(y, pole_site)), det(pole_site, x, y)

    assert len(edges) == 135
    clocks = Counter()
    for (i, j), distance in edges.items():
        a, b = flip[i], flip[j]
        assert edges[min(a, b), max(a, b)] == distance
        real, imag = phase_numerator(i, j)
        assert phase_numerator(a, b) == (real, -imag)
        if distance in (724, 804):
            assert a < b
            clocks[distance] += 1
    assert clocks == {724: 6, 804: 6}
    ordinary = set(range(19)) - {pole}
    orbit = {i: min(i, shift[i], shift[shift[i]]) for i in ordinary}
    masks = {frozenset(mask) for mask in combinations(ordinary, 3) if len({orbit[i] for i in mask}) == 2}
    assert len(masks) == 270 and {frozenset(flip[i] for i in mask) for mask in masks} == masks


def replay(model, storage, start, stop, kernel):
    limit = DEPTH
    pump = run_path(str(model / 'ch39-native-light-projector-pump.py'))
    shell = run_path(str(model / 'ch39-three-ray-chiral-boundary.py'))
    ctx.prec = 160
    verify_reflection(shell)
    sites = shell['geo']['sites']
    pole = shell['soft']['pole']
    cycle = [sites.index((y, z, x)) for x, y, z in sites]
    omega = acb(-1, -arb(3).sqrt()) / 2
    weights, sectors = [], []
    for k, frame in enumerate(pump['rotation_bases'](shell)):
        net = np.load(storage / f'character-{k}.npy', mmap_mode='r')
        assert net.shape == (CELLS, CELLS, 4, 4, frame.ncols(), 5, 2)
        assert net.dtype == np.int32
        ordinary = [j for j in range(frame.ncols()) if frame[pole, j].is_zero()]
        representatives = [min(i for i in range(19) if not frame[i, j].is_zero()) for j in ordinary]
        order = [sites[v] for i in representatives for v in (i, cycle[i], cycle[cycle[i]])]
        if k == 0:
            ordered = order
            trivial = net
            pole_column = next(j for j in range(frame.ncols()) if j not in ordinary)
            assert frame[pole, pole_column] == 1
        assert ordered == order
        sectors.append((net, ordinary))
        weight = []
        for r in range(3):
            coefficient = omega**(k*r) / arb(3).sqrt()
            parts = (coefficient.real, coefficient.imag)
            integers = [round(float(part * SCALE)) for part in parts]
            assert all(abs(part * SCALE - value) < arb('.5000000001')
                       for part, value in zip(parts, integers))
            weight.append(integers)
        weights.append(weight)
    weights = np.array(weights, dtype=np.int64)
    cases = pole_cases(ordered, 2)
    settings = json.loads((storage / 'transfer.json').read_text())
    transfer = arb(settings['transfer'])
    requested = arb(FLOOR)
    assert requested >= 0
    threshold = settings['threshold']
    if requested > 0:
        threshold = max(threshold, int(float((transfer + arb('1.5') * requested).upper()) * SCALE) + 2)
    accounted = settings['accounted_splits']
    floor = (arb(threshold) / SCALE - transfer) / arb('1.5')
    assert floor > requested

    def threshold_at(depth):
        return threshold + TRANSFER_UNIT_STEP * max(0, depth - accounted)

    def check(controls, depth):
        threshold = threshold_at(depth)
        if kernel is None:
            return verify_nets(controls, SCALE, threshold)
        return native_check(kernel, controls, SCALE, threshold)

    def visit(controls, origins, depth, counts, failures):
        # Coordinates are exact ticks on the finest dyadic mesh and are only
        # needed for locating unresolved leaves.
        passed = check(controls, depth)
        counts[depth] += int(passed.sum())
        if passed.all():
            return
        remaining = controls[~passed]
        origins = origins[~passed]
        if depth == limit:
            failures.extend(map(tuple, origins.tolist()))
            return
        stride = 1 << (limit - 1 - depth)
        for y, children in enumerate(split_cubic(remaining, 1)):
            for x, child in enumerate(split_cubic(children, 2)):
                visit(child, origins + np.array([y, x]) * stride, depth + 1, counts, failures)

    print(json.dumps({'start': start, 'stop': stop, 'depth': limit, 'floor': floor.str(12),
                      'cyclic_mask_classes': len(cases)}), flush=True)
    all_certified = True
    for y in range(start, stop):
        nets = [net[y:y+1][..., ordinary, :, :] for net, ordinary in sectors]
        readings = np.zeros((CELLS, 4, 4, 19, 15, 2), dtype=np.int32)
        readings[..., :18, :, :] = physical_line(nets, weights, 0, SCALE)
        readings[..., 18, :5, :] = trivial[y, ..., pole_column, :, :]
        origins = np.column_stack((np.full(CELLS, y), np.arange(CELLS))) << limit
        counts = [0] * (limit + 1)
        open_cases = []
        region_size = CELLS * 4**limit
        for names, keep in cases:
            local = [0] * (limit + 1)
            failures = []
            visit(readings[..., keep, :, :], origins, 0, local, failures)
            assert sum(n * 4**(limit-d) for d, n in enumerate(local)) + len(failures) == region_size
            counts = [a+b for a, b in zip(counts, local)]
            if failures:
                open_cases.append({'mask': names, 'count': len(failures), 'first_ticks': sorted(failures)[:10]})
        all_certified &= not open_cases
        print(json.dumps({'row': y, 'certified_by_depth': counts, 'unresolved': sum(c['count'] for c in open_cases),
                          'open_cases': open_cases}), flush=True)
    return all_certified


if __name__ == '__main__':
    assert __debug__
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('model', type=Path, help='pinned shell-model replay directory')
    parser.add_argument('storage', type=Path, help='directory for generated control nets')
    parser.add_argument('--generate', action='store_true', help='create the control nets before replay')
    parser.add_argument('--generate-only', action='store_true', help='create the control nets and exit')
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int, default=CELLS // 2)
    parser.add_argument('--native', type=Path, help='optional compiled exact-integer kernel')
    args = parser.parse_args()
    assert 0 <= args.start < args.stop <= CELLS // 2
    if args.generate or args.generate_only:
        generate(args.model, args.storage)
    if not args.generate_only:
        kernel = load_native_kernel(args.native) if args.native is not None else None
        if not replay(args.model, args.storage, args.start, args.stop, kernel):
            raise SystemExit('unresolved angular boxes remain')
