"""Check blind phase boxes for every fifteen-site layout and every allowed spoke deficit."""

from itertools import combinations
import json
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx
import numpy as np


def omissions(shell):
    sites = shell["geo"]["sites"]
    rotation = [sites.index((b, c, a)) for a, b, c in sites]
    triples = {
        tuple(sorted((i, rotation[i], rotation[rotation[i]])))
        for i in range(19)
        if rotation[i] != i
    }

    def rotate(mask):
        return tuple(sorted(rotation[i] for i in mask))

    masks = list(combinations(range(19), 4))
    representatives = {min(mask, rotate(mask), rotate(rotate(mask))) for mask in masks}
    inherited = {
        mask for mask in representatives if any(set(o).issubset(mask) for o in triples)
    }
    assert (len(masks), len(representatives), len(inherited)) == (3876, 1296, 36)
    return representatives - inherited


def matrices(pump, shell, charges):
    frames = list(pump["rotation_bases"](shell))
    adj = pump["adjoint"]
    pole = shell["soft"]["pole"]
    star = acb_mat(19, 19)
    for site in shell["soft"]["star"]["orbit"]:
        star[pole, pole] -= 1
        star[site, site] -= 1
        star[pole, site] += 1
        star[site, pole] += 1
    blocks = []
    for frame in frames:
        drift = adj(frame) * star * frame
        for q in charges:
            base, links = pump["clock_matrices"](shell, q)
            h, a, b = [adj(frame) * m * frame for m in (base, *links)]
            blocks.append((frame, acb(0, q), h, a, b, drift))
    return blocks


def synthesis(point, delta, keep, blocks, pump, common, pivots=None):
    adj = pump["adjoint"]
    columns = [[], [], [], []]
    chosen = []
    for i, (frame, iq, base, first, second, drift) in enumerate(blocks):
        a = (iq * point[0]).exp() * first
        b = (iq * point[1]).exp() * second
        h = base + a + adj(a) + b + adj(b) + delta * drift
        if pivots is None:
            line = pump["first_line"](h)[0]
            pivot = max(range(len(line)), key=lambda j: float(abs(line[j])))
        else:
            pivot = pivots[i]
        chosen.append(pivot)
        slopes = [iq * (a - adj(a)), iq * (b - adj(b)), drift]
        vectors = common["eigenline_derivatives"](h, slopes, pump["first_line"], pivot)
        for out, v in zip(columns, vectors, strict=True):
            values = (frame * acb_mat([[z] for z in v])).entries()
            out.append([values[j] for j in keep])
    return [acb_mat(cols).transpose() for cols in columns], chosen


def transforms(matrix):
    # Arbitrary fixed coordinate changes are valid once interval determinants verify invertibility.
    mid = np.array([[complex(z) for z in row] for row in matrix.tolist()])
    u, s, vh = np.linalg.svd(mid)
    left = u.conj().T
    left[:-1] /= s[:-1, None]

    def exact(a):
        return acb_mat([[acb(z.real, z.imag) for z in row] for row in a])

    left, right = exact(left), exact(vh.conj().T)
    assert not left.det().contains(0) and not right.det().contains(0)
    return left, right


def schur(mats, left, right):
    matrix, *derivatives = [left * m * right for m in mats]
    block = acb_mat([[matrix[i, j] for j in range(14)] for i in range(14)])
    row = acb_mat([[matrix[14, i] for i in range(14)]])
    col = acb_mat([[matrix[i, 14]] for i in range(14)])
    assert not block.det().contains(0)
    x = block.solve(col)
    y = block.transpose().solve(row.transpose()).transpose()
    f = matrix[14, 14] - (row * x)[0, 0]
    covector = acb_mat([[-y[0, i] for i in range(14)] + [1]])
    vector = acb_mat([[-x[i, 0]] for i in range(14)] + [[1]])
    dx, dy, ds = [(covector * d * vector)[0, 0] for d in derivatives]
    return (
        arb_mat([[f.real], [f.imag]]),
        arb_mat([[dx.real, dy.real], [dx.imag, dy.imag]]),
        arb_mat([[ds.real], [ds.imag]]),
    )


def verify(root):
    here = Path(__file__).resolve().parent
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    common = run_path(str(here / "ch39-common-charge-amplitudes.py"))
    run_path(str(here / "ch39-five-charge-pole-readout.py"))["verify"](root)
    ctx.prec = 192
    expected = omissions(shell)
    entries = json.loads(Path(__file__).with_suffix(".json").read_text())
    assert len(entries) == len(expected) and {tuple(e[0]) for e in entries} == expected
    blocks = matrices(pump, shell, common["CHARGES"])
    delta_max = arb("1e-8")
    delta = arb(0).union(delta_max)
    eye = arb_mat([[1, 0], [0, 1]])
    worst = arb(0)
    for mask, a, b, r in entries:
        keep = [j for j in range(19) if j not in mask]
        center = list(map(arb, (a, b)))
        radius = arb(r)
        assert 0 < radius < arb("3e-8")
        mids, pivots = synthesis(center, arb(0), keep, blocks, pump, common)
        left, right = transforms(mids[0])
        value, jac, _ = schur(mids, left, right)
        pre = jac.mid().inv().mid()
        assert not pre.det().contains(0)
        box = [z + arb(0, radius.upper()) for z in center]
        whole, _ = synthesis(box, delta, keep, blocks, pump, common, pivots)
        _, jac, drift = schur(whole, left, right)
        move, rate = pre * value, pre * drift
        defect = eye - pre * jac
        ratio = max(
            (move[i, 0].abs_upper() + delta_max * rate[i, 0].abs_upper()) / radius
            + sum(z.abs_upper() for z in row)
            for i, row in enumerate(defect.tolist())
        )
        assert ratio < arb(".91"), (mask, ratio)
        worst = max(worst, ratio.upper())
    print(
        "four-site omissions, rotation classes, inherited, Schur boxes, worst self-map bound",
        3876,
        1296,
        36,
        len(entries),
        worst.str(12),
    )


if __name__ == "__main__":
    if len(sys.argv) != 2 or not __debug__:
        raise SystemExit("Run without -O and supply the pinned replay directory.")
    verify(Path(sys.argv[1]).resolve())
