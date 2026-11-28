"""Certify full-torus coverage for two lost stations in different orbits."""

from itertools import combinations, product
from math import isqrt
from pathlib import Path
from runpy import run_path
from tempfile import TemporaryDirectory
import argparse

from flint import acb, arb, ctx
import numpy as np


LEFT_SCALE = 256
LEFT_CAP = 100_000
DEFECT_SHIFT = 2**15
CELLS = 1024
TRANSFER_UNITS = 3


def control_nets(base, verifier, cache):
    full = run_path(str(verifier / "ch39-five-charge-full-torus.py"))
    pump = run_path(str(base / "ch39-native-light-projector-pump.py"))
    shell = run_path(str(base / "ch39-three-ray-chiral-boundary.py"))
    pump["verify"](base)
    ctx.prec = 160
    frames = list(pump["rotation_bases"](shell))
    pole = shell["soft"]["pole"]
    omega = acb(-1, -arb(3).sqrt()) / 2
    normalization = arb(3).sqrt()
    scale = 2**30
    control_globals = full["fill_controls"].__globals__
    control_globals["SCALE"] = scale
    build_controls = control_globals["controls"]

    def checked_controls(nodes):
        controls = build_controls(nodes)
        assert -scale < controls.min() and controls.max() < scale
        return controls

    control_globals["controls"] = checked_controls
    weights = []
    error = []
    nets = []
    for k in range(3):
        net_path = cache / f"character-{k}.npy"
        shape = (CELLS, CELLS, 3, 3, 6, 5, 2)
        net = np.lib.format.open_memmap(net_path, mode="w+", dtype=np.int32, shape=shape)
        interpolation, off_line = full["fill_controls"](pump, shell, k, net)
        net.flush()
        error.append(interpolation + arb(5).sqrt() * off_line)
        nets.append(net)
        row_indices = [i for i in range(frames[k].ncols()) if frames[k][pole, i].is_zero()]
        rows = [[i for i in range(19) if not frames[k][i,j].is_zero()] for j in row_indices]
        sites = shell["geo"]["sites"]
        shift = [sites.index((y,z,x)) for x,y,z in sites]
        orbits = [(min(row), shift[min(row)], shift[shift[min(row)]]) for row in rows]
        if k == 0:
            orbit_sites = [tuple(sites[i]) for orbit in orbits for i in orbit]
        else:
            assert orbit_sites == [tuple(sites[i]) for orbit in orbits for i in orbit]
        phases = [omega**(k*j) / normalization for j in range(3)]
        integers = []
        for value in phases:
            parts = [int(round(float(c * scale))) for c in (value.real, value.imag)]
            assert all(abs(c * scale - p) < arb("0.5000000001") for c, p in zip((value.real, value.imag), parts))
            integers.append(parts)
        weights.append(integers)
        print("character", k, "interpolation", interpolation.str(9), "line error", off_line.str(9), flush=True)
    return full, orbit_sites, nets, np.array(weights, dtype=np.int64), max(item.upper() for item in error)


def physical_line(nets, weights, line, scale):
    controls = np.empty((CELLS, 3, 3, 18, 15, 2), dtype=np.int32)
    for k, (net, weight) in enumerate(zip(nets, weights)):
        block = net[line].astype(np.int64)
        re = block[..., 0][..., None, :]
        im = block[..., 1][..., None, :]
        wr = weight[:, 0, None]
        wi = weight[:, 1, None]
        physical = controls[..., 5*k:5*(k+1), :].reshape(CELLS, 3, 3, 6, 3, 5, 2)
        real = (re * wr - im * wi + scale // 2) // scale
        imag = (re * wi + im * wr + scale // 2) // scale
        # This bound also keeps the sums in int32 during subdivision.
        assert -scale < real.min() and real.max() < scale
        assert -scale < imag.min() and imag.max() < scale
        physical[..., 0] = real
        physical[..., 1] = imag
    return controls


def verify_nets(balls, scale, threshold):
    center = (balls[:, 1, 1, ..., 0] + 1j * balls[:, 1, 1, ..., 1]) / scale
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
    for j, k in product(range(3), repeat=2):
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


def unresolved(balls, origin, width, split, threshold, scale, limit, depth=0):
    fail = ~verify_nets(balls, scale, threshold)
    if not fail.any():
        return []
    if depth == limit:
        return [(float(y), float(x), width) for y, x in origin[fail]]
    branches = []
    for iy, v in enumerate(split(balls[fail], 1)):
        for ix, q in enumerate(split(v, 2)):
            offset = width * np.array([iy, ix]) / 2
            branches.extend(unresolved(q, origin[fail]+offset, width/2, split, threshold, scale, limit, depth+1))
    return branches


def orbit_representatives(sites):
    assert len(sites) == len(set(sites)) == 18
    orbits = [sites[start:start + 3] for start in range(0, len(sites), 3)]
    for orbit in orbits:
        for index, (x, y, z) in enumerate(orbit):
            assert orbit[(index + 1) % 3] == (y, z, x)
    rotation = {site: orbit[(j + 1) % 3] for orbit in orbits for j, site in enumerate(orbit)}
    hard_pairs = [((1, 5, 29), (13, 23, 13)), ((7, 17, 23), (7, 23, 17))]
    representatives = []
    covered = set()
    for first, second in combinations(orbits, 2):
        for offset in range(3):
            equivalent = {
                frozenset((first[j], second[(j + offset) % 3])) for j in range(3)
            }
            assert len(equivalent) == 3 and not covered.intersection(equivalent)
            covered.update(equivalent)
            representative = next(
                (pair for pair in hard_pairs if frozenset(pair) in equivalent),
                (first[0], second[offset]),
            )
            assert {
                frozenset(representative),
                frozenset(rotation[site] for site in representative),
                frozenset(rotation[rotation[site]] for site in representative),
            } == equivalent
            representatives.append(representative)
    assert len(representatives) == 45 and len(covered) == 135
    return representatives


def main(base, verifier, cache):
    full, sites, nets, weights, error = control_nets(base, verifier, cache)
    scale = full["fill_controls"].__globals__["SCALE"]
    limit = 8
    entry_rounding = (13 + 4*limit + TRANSFER_UNITS)
    error += arb(480).sqrt() * entry_rounding / scale
    threshold = 4290
    assert error < arb("3.993061e-6")
    normalized_floor = (arb(threshold) / scale - error) / arb("1.5")
    assert normalized_floor > arb("1.5e-9")
    print("uniform trial transfer", error.str(10), "test threshold", threshold, scale, flush=True)
    named_pairs = orbit_representatives(sites)
    keep = [[i for i, site in enumerate(sites) if site not in pair] for pair in named_pairs]
    assert all(len(k) == 16 for k in keep)
    h = float(2*np.pi/CELLS)
    failures = [[] for _ in keep]
    for y in range(CELLS):
        controls = physical_line(nets, weights, y, scale)
        origin = np.column_stack((np.full(CELLS, y*h), np.arange(CELLS)*h))
        for i, rows in enumerate(keep):
            v = controls[..., rows, :, :]
            failures[i].extend(unresolved(v, origin, h, full["split_axis"], threshold, scale, limit))
        if y % 64 == 63:
            print("rows", y+1, "unresolved", [len(v) for v in failures], flush=True)
    for names, failed in zip(named_pairs, failures):
        print(names, "coverage at integer threshold", "unresolved", len(failed), failed[:2], flush=True)
    if any(failures):
        raise SystemExit("uncovered boxes remain")
    print("normalized singular-value lower bound", normalized_floor.str(10), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path, help="replay directory at cb58c8e64ecea67184e0c1d68a690be098a7a384")
    parser.add_argument("verifier", type=Path, help="replay directory at 1c485355843e787bf64af78c6c21b0dffd0fbd78")
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without -O because the verification uses assertions")
    with TemporaryDirectory(prefix="station-control-nets-") as temporary:
        main(args.model.resolve(), args.verifier.resolve(), Path(temporary))
