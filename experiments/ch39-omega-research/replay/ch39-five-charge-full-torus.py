"""Enclose the five-charge synthesis on the full common clock torus."""

from itertools import product
from math import isqrt
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import acb, acb_mat, arb, ctx

CHARGES = (1, -4, 2, -3, 6)
DEGREE = 10
CELLS = 256
SCALE = 2**25
LEFT_SCALE = 1024


def matrix_array(matrix):
    return np.array([[complex(z) for z in row] for row in matrix.tolist()])


def candidates(base, a, b):
    """Only proposes dyadic Laurent coefficients. The subsequent tests verify them."""
    base, a, b = map(matrix_array, (base, a, b))
    roots = np.exp(2j * np.pi * np.arange(64) / 64)
    links = roots[:, None, None, None] * a + roots[None, :, None, None] * b
    energy, vectors = np.linalg.eigh(base + links + links.conj().swapaxes(-2, -1))
    v = vectors[..., 0]
    anchor = np.linalg.eigh(base)[1][:, 0]
    trial = v * np.sum(v.conj() * anchor, axis=-1)[..., None]
    indices = np.arange(-DEGREE, DEGREE + 1) % 64

    def coeff(field):
        table = np.fft.fft2(field, axes=(0, 1)) / 64**2
        return table[indices[:, None], indices]

    return coeff(trial), coeff(energy[..., 0])


def coefficients(array):
    return np.array(
        [acb(float(z.real), float(z.imag)) for z in array.flat], dtype=object
    ).reshape(array.shape)


def sample(table, count, charge=1):
    # The forward transform at negative frequencies gives unscaled positive phases.
    frequencies = [(-charge * i) % count for i in range(-DEGREE, DEGREE + 1)]
    assert len(set(frequencies)) == len(frequencies)
    zero = acb(0)

    def transform(values):
        v = [zero] * count
        for i, z in zip(frequencies, values, strict=True):
            v[i] = z
        return acb.dft(v)

    rows = [transform(row) for row in table]
    return [transform(col) for col in zip(*rows, strict=True)]  # [b][a]


def forward2(grid):
    rows = [acb.dft(row) for row in grid]
    return [acb.dft(col) for col in zip(*rows, strict=True)]


def verify_line(pump, operators, q):
    base, a, b = operators
    trial, energy = candidates(base, a, b)
    v, e = coefficients(trial), coefficients(energy)
    n = base.nrows()
    size = np.empty(e.shape, dtype=object)
    for i, j in np.ndindex(size.shape):
        size[i, j] = sum(abs(z) ** 2 for z in v[i, j]).sqrt()
    total = sum(size.flat)
    assert total < arb("1.5") and 2 * size[DEGREE, DEGREE] - total > arb(".4")
    # Rayleigh parameter is real by symmetrization, including for nonsymmetric proposals.
    e = (
        e
        + np.array([z.conjugate() for z in e[::-1, ::-1].flat], dtype=object).reshape(
            e.shape
        )
    ) / 2
    nodes = 64
    roots = [acb(0, 2 * arb.pi() * i / nodes).exp() for i in range(nodes)]
    adj = pump["adjoint"]
    ax = [a * z + adj(a * z) for z in roots]
    by = [b * z + adj(b * z) for z in roots]
    values = [sample(v[:, :, l], nodes) for l in range(n)]
    energies = sample(e, nodes)
    residual = [[[] for _ in roots] for _ in range(n)]
    for y in range(nodes):
        for x in range(nodes):
            column = acb_mat([[w[y][x]] for w in values])
            r = (base + ax[x] + by[y]) * column - energies[y][x] * column
            for l in range(n):
                residual[l][y].append(r[l, 0])
    # Bidegree <= 20, strictly below 64/2: the DFT recovers all residual coefficients.
    residual_coeff = [forward2(grid) for grid in residual]
    bound = (
        sum(
            sum(t[i][j].abs_upper() ** 2 for t in residual_coeff).sqrt()
            for i, j in product(range(nodes), repeat=2)
        )
        / nodes**2
    )
    assert bound < arb("5e-7"), (q, bound)
    lam = (base + ax[0] + by[0]).eig(algorithm="rump")
    ground = min(z.real.lower() for z in lam).union(min(z.real.upper() for z in lam))
    assert abs(energies[0][0] - ground) < arb(".1")
    third = sum(
        size[i, j] * (abs(q) * (abs(i - DEGREE) + abs(j - DEGREE))) ** 3
        for i, j in np.ndindex(size.shape)
    )
    return v, third, bound.upper()


def quantized(grid):
    # Verify containment in an integer box in units of 1/SCALE before discarding balls.
    out = np.empty((len(grid), len(grid)), dtype=np.int64)
    for y, row in enumerate(grid):
        for x, t in enumerate(row):
            v = t * SCALE
            m = round(float(v))
            if not (abs(v - m) < 1 and abs(m) < 2 * SCALE):
                raise ArithmeticError("nodal enclosure")
            out[y, x] = m
    return out


def controls(nodes):
    # Quadratic Bernstein controls from the values at 0, 1/2 and 1 on both axes.
    def axis(data, direction):
        lo = np.take(data, np.arange(0, 2 * CELLS, 2), axis=direction)
        mid = np.take(data, np.arange(1, 2 * CELLS, 2), axis=direction)
        hi = np.roll(lo, -1, axis=direction)
        return np.stack((lo, 2 * mid - (lo + hi) // 2, hi), axis=-1)

    return axis(axis(nodes, 0), 1)


def split_axis(data, axis):
    a, m, b = np.moveaxis(data, axis, 0)
    u, v = (a + m) // 2, (m + b) // 2
    w = (u + v) // 2
    return [np.stack(s, axis=axis) for s in ((a, u, w), (w, v, b))]


def verify_cells(balls, level=0):
    # [cells, y control, x control, six rows, five columns, real/imaginary]
    length = len(balls)
    central_control = balls[:, 1, 1]
    mid = (central_control[..., 0] + 1j * central_control[..., 1]) / SCALE
    proposed = np.linalg.pinv(mid)
    assert np.isfinite(proposed).all() and np.abs(proposed).max() < 10000
    left = np.rint(
        LEFT_SCALE * np.stack((proposed.real, proposed.imag), axis=-1)
    ).astype(np.int64)
    # The caps bound a complex dot product by 12*(10000*1024)*(2*SCALE)
    # and either squared left norm by 60*(10000*1024)**2, within int64.
    assert np.abs(left).max() <= 10000 * LEFT_SCALE and np.abs(balls).max() <= 2 * SCALE
    lreal, limag = left[..., 0], left[..., 1]
    norm2 = (left * left).sum(axis=(1, 2, 3))
    norms = np.fromiter(
        (isqrt(int(v)) + 1 for v in norm2), dtype=np.int64, count=length
    )
    unit = SCALE * LEFT_SCALE
    eye = np.eye(5, dtype=np.int64) * unit
    maxerr = np.zeros(length, dtype=np.int64)
    for j, k in product(range(3), repeat=2):
        r, i = balls[:, j, k, ..., 0], balls[:, j, k, ..., 1]
        re = lreal @ r - limag @ i - eye
        im = lreal @ i + limag @ r
        maxerr = np.maximum(maxerr, np.maximum(abs(re), abs(im)).max(axis=(1, 2)))
    # Frobenius norm <= sqrt(50) * largest real/imaginary entry of the defect.
    # Demand 1/4096 for every convex combination of the represented controls.
    ok = 8 * maxerr + norms * (SCALE // 4096) < unit
    counts = [int(ok.sum())]
    if ok.all():
        return counts
    assert level < 9
    failed = balls[~ok]
    for half in split_axis(failed, 1):
        for quarter in split_axis(half, 2):
            for d, n in enumerate(verify_cells(quarter, level + 1), 1):
                if d == len(counts):
                    counts.append(0)
                counts[d] += n
    return counts


def verify_character(pump, shell, k):
    basis = list(pump["rotation_bases"](shell))[k]
    pole = shell["soft"]["pole"]
    rows = [j for j in range(basis.ncols()) if basis[pole, j].is_zero()]
    assert len(rows) == 6
    nodes = 2 * CELLS
    control_net = np.empty((CELLS, CELLS, 3, 3, 6, 5, 2), dtype=np.int64)
    third_bounds = []
    residuals = []
    for col, q in enumerate(CHARGES):
        base, links = pump["clock_matrices"](shell, q)
        ops = [pump["adjoint"](basis) * m * basis for m in (base, *links)]
        v, third, residual = verify_line(pump, ops, q)
        third_bounds.append(third)
        residuals.append(residual)
        print("charge, residual bound", q, residual.str(6), flush=True)
        for r, j in enumerate(rows):
            grid = sample(v[:, :, j], nodes, q)
            for part in range(2):
                part_grid = [[z.imag if part else z.real for z in row] for row in grid]
                control_net[:, :, :, :, r, col, part] = controls(quantized(part_grid))
    h = arb.pi() / CELLS
    third = sum(s**2 for s in third_bounds).sqrt()
    # Tensor interpolation operator norm <= 5/4 on either coordinate.
    interpolation = (1 + arb(5) / 4) * third * h**3 / (9 * arb(3).sqrt())
    # Integer node error < 1 unit. Each control conversion has l1 norm <= 3
    # and rounding < 1 unit. Each dyadic subdivision adds at most four units.
    rounding = arb(60).sqrt() * (13 + 4 * 9) / SCALE
    assert interpolation + rounding < arb(".000177")
    counts = []
    for strip in control_net:
        for d, n in enumerate(verify_cells(strip)):
            if d == len(counts):
                counts.append(0)
            counts[d] += n
    assert sum(
        n * 4 ** (len(counts) - 1 - d) for d, n in enumerate(counts)
    ) == CELLS**2 * 4 ** (len(counts) - 1)
    # Each of the three softened edge Laplacians has norm two. The pinned
    # gap is > 1 at every parameter. Continuity keeps the selected Ritz energy
    # on the ground branch from (0,0,0), since the residual/norm is < 1/2.
    full_residual = max(residuals) + 6 * arb("1e-8") * arb("1.5")
    distance = full_residual / arb(".4")
    assert distance < arb(".5")
    off_line = full_residual / (1 - distance)
    floor = (arb(1) / 4096 - interpolation - rounding - arb(5).sqrt() * off_line) / arb(
        "1.5"
    )
    assert floor > arb(".00004")
    print(
        "character, interpolation, integer enclosure, normalized lower bound, leaves by depth",
        k,
        interpolation.str(8),
        rounding.str(8),
        floor.str(8),
        counts,
        flush=True,
    )


if __name__ == "__main__":
    if len(sys.argv) != 2 or not __debug__:
        raise SystemExit(
            "Run without -O, with the pinned replay directory as the single argument."
        )
    root = Path(sys.argv[1]).resolve()
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    pump["verify"](root)
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 160
    for character in range(3):
        verify_character(pump, shell, character)
