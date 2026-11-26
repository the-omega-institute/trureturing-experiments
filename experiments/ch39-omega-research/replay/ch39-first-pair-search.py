"""在前两枚实际零点的矩阵弦上探索全部不带钟的边对。浮点搜寻。"""
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import acb, ctx

CHARGES = (-3, -4, 2, 1, 6)
TIMES = np.linspace(0., 1., 37)
FRACTIONS = np.r_[1e-7, .0001, .001, np.linspace(.005, .995, 37), .999, .9999, 1-1e-7]


def prepare(shell, edges, heights, charge):
    first, second = [np.array(shell["laplacian"](charge, height).tolist(), complex)
                     for height in heights]
    columns = np.zeros((19, len(edges)), complex)
    for k, (i, j) in enumerate(edges):
        columns[i, k] = 1
        columns[j, k] = first[j, i]/abs(first[j, i])
    return first, second-first, columns


def resolvent_grid(first, direction, columns, time):
    levels, vectors = np.linalg.eigh(first + time*direction)
    projected = vectors.conj().T @ columns
    grids = []
    for count in (2, 3, 4):
        energies = levels[count-1] + (levels[count]-levels[count-1])*FRACTIONS
        factors = (levels[None, :] - energies[:, None])[:, None, :]
        resolvents = (projected.conj().T[None, :, :]/factors) @ projected
        grids.append((resolvents, energies))
    return grids


def refine(table, edge_indices, edges, time, energy, time_lo, time_hi):
    first, direction, columns = table
    indices = list(edge_indices)
    b = columns[:, indices]
    for _ in range(12):
        kernel = np.linalg.solve(first + time*direction - energy*np.eye(19), b)
        g = b.conj().T @ kernel
        dt = -(kernel.conj().T @ direction @ kernel)[0, 1]
        de = (kernel.conj().T @ kernel)[0, 1]
        try:
            step = np.linalg.solve([[dt.real, de.real], [dt.imag, de.imag]],
                                   [g[0, 1].real, g[0, 1].imag])
        except np.linalg.LinAlgError:
            return
        time -= step[0]
        energy -= step[1]
        if not time_lo-.02 < time < time_hi+.02:
            return
        if sum(abs(step)) < 1e-10:
            break
    base = first + time*direction
    kernel = np.linalg.solve(base - energy*np.eye(19), b)
    g = b.conj().T @ kernel
    if abs(g[0, 1]) > 1e-8 or not 0 < time < 1 or np.any(abs(g.diagonal()) < 1e-12):
        return
    increments = -1/g.diagonal().real
    weights = np.array([abs(base[edges[k]]) for k in indices]) + increments
    if min(weights) <= 0 or sum(np.linalg.eigvalsh(base) < energy) - sum(increments > 0) != 2:
        return
    frame = np.linalg.qr(kernel)[0]
    velocities = []
    for derivative in (direction, np.outer(b[:, 0], b[:, 0].conj()),
                       np.outer(b[:, 1], b[:, 1].conj())):
        compressed = frame.conj().T @ derivative @ frame
        velocities.append((compressed[0, 1].real, -compressed[0, 1].imag,
                           .5*(compressed[0, 0]-compressed[1, 1]).real))
    determinant = np.linalg.det(velocities)
    if abs(determinant) > 1e-9:
        return time, energy, weights.tolist(), increments, determinant


def search(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 80
    edges = [(i, j) for i, j, d, _ in shell["geo"]["edges"] if d not in (724, 804)]
    heights = [acb.zeta_zero(j).imag for j in (1, 2)]
    tables = {q: prepare(shell, edges, heights, q) for q in CHARGES}
    print("无钟边、无序边对、时间切片、区内能量样本", len(edges),
          len(edges)*(len(edges)-1)//2, len(TIMES), len(FRACTIONS), flush=True)
    for charge in CHARGES:
        table = tables[charge]
        previous = resolvent_grid(*table, TIMES[0])
        seeded = 0
        found = 0
        for k in range(1, len(TIMES)):
            t_lo, t_hi = TIMES[k-1], TIMES[k]
            current = resolvent_grid(*table, t_hi)
            for (low, e), (high, f) in zip(previous, current):
                corners = (low[:-1], low[1:], high[1:], high[:-1])
                reals = [z.real for z in corners]
                imags = [z.imag for z in corners]
                may_vanish = ((np.minimum.reduce(reals) < 0) & (np.maximum.reduce(reals) > 0)
                              & (np.minimum.reduce(imags) < 0) & (np.maximum.reduce(imags) > 0))
                for level, i, j in zip(*np.where(may_vanish)):
                    if i >= j:
                        continue
                    seeded += 1
                    e_middle = (e[level]+e[level+1]+f[level]+f[level+1])/4
                    root_data = refine(table, (i, j), edges, (t_lo+t_hi)/2, e_middle, t_lo, t_hi)
                    if root_data is None:
                        continue
                    time, energy, weights, increments, determinant = root_data
                    found += 1
                    print("候选", charge, edges[i], edges[j], time, energy, weights,
                          determinant, flush=True)
                    for q, (a, d, all_columns) in tables.items():
                        if q == charge:
                            continue
                        b = all_columns[:, [i, j]]
                        levels = np.linalg.eigvalsh(a + time*d + (b*increments) @ b.conj().T)
                        print("旁观荷、交点能量下态数、九以下态数", q,
                              sum(levels < energy), sum(levels < 9), flush=True)
            previous = current
        print("荷、牛顿初值数、通过筛选数", charge, seeded, found, flush=True)


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    search(root)
