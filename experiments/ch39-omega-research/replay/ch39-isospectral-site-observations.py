"""逐点返回矩在三个固定谱流形上的局部辨识阶数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


# 七个站轨按最小站号排列。在每轨的一至八阶行中取如下 26 行。
INVARIANT_ROWS = (0, 1, 2, 8, 9, 10, 11, 14, 24, 25, 28, 29, 33,
                  34, 35, 37, 38, 40, 41, 42, 43, 48, 49, 50, 53, 54)


def site_orbits(geometry):
    sites = geometry["sites"]
    rotation = [sites.index((y, z, x)) for x, y, z in sites]
    for i, j, distance, _ in geometry["edges"]:
        if distance in (724, 804):
            assert rotation[i] < rotation[j]
    return [i for i in range(len(sites)) if i == min(i, rotation[i], rotation[rotation[i]])], rotation


def matrix(weights, groups, phases):
    result = acb_mat(19, 19)
    for weight, group in zip(weights, groups):
        for i, j in group:
            result[i, i] += weight
            result[j, j] += weight
            result[i, j] = -weight * phases[i, j]
            result[j, i] = result[i, j].conjugate()
    return result


def eigen_derivatives(mid, groups, phases):
    values, vectors = mid.eig(right=True, algorithm="rump")
    order = sorted(range(19), key=lambda k: float(values[k].real))
    eigenvalues = [values[k].real for k in order]
    assert all(z.imag.contains(0) and z.real > 0 for z in values)
    gap = min((b-a).lower() for a, b in zip(eigenvalues, eigenvalues[1:]))
    assert gap > 0
    rows = []
    for k in order:
        norm = sum(abs(vectors[i, k])**2 for i in range(19))
        rows.append([
            sum(abs(vectors[i, k] - phases[i, j]*vectors[j, k])**2 for i, j in group)/norm
            for group in groups
        ])
    return eigenvalues, rows, gap


def newton_box(weights, variables, groups, phases, target):
    for pos, radius in enumerate(map(arb, ("2e-11", "1e-22"))):
        levels, jac, gap = eigen_derivatives(matrix(weights, groups, phases), groups, phases)
        inverse = arb_mat([[row[g] for g in variables] for row in jac]).inv()
        infnorm = max(sum(abs(inverse[i,j]) for j in range(19)).upper() for i in range(19))
        error = inverse * arb_mat([[a-b] for a, b in zip(levels, target)])
        correction = max(abs(x).upper() for x in error.entries())
        change = 6*19*radius
        assert gap > 2*change
        contraction = 19*infnorm*12*change/(gap-2*change)
        assert contraction < 1 and correction + contraction*radius < radius
        if pos == 0:
            weights = weights.copy()
            for i, g in enumerate(variables):
                old = weights[g]
                weights[g] = (old - error[i,0]).mid()
                assert abs(weights[g] - old) + arb("1e-22") < radius
        else:
            error_bound = 12*change/(gap-2*change)
            spectral = [[a + arb(0,error_bound) for a in row] for row in jac]
            boxes = [w + arb(0, radius) if g in variables else w for g,w in enumerate(weights)]
            assert all(w > 0 for w in boxes)
            return matrix(boxes, groups, phases), spectral


def site_moments(laplacian, groups, phases, reps, rotate):
    identity = acb_mat([[int(i == j) for j in range(19)] for i in range(19)])
    scaled = (laplacian - 40*identity) / 40
    powers = [identity]
    for _ in range(7):
        powers.append(powers[-1]*scaled)
    omega = acb(-1, arb(3).sqrt())/2
    invariant, nontrivial = [], []
    for r in reps:
        columns = []
        for group in groups:
            edges = []
            for i, j in group:
                # 单边导数等于 bb*，其中 b=e_i-conjugate(phase_ij)e_j。
                projected = [p[r,i]-p[r,j]*phases[i,j].conjugate() for p in powers]
                edges.append([
                    sum((projected[t]*projected[k-1-t].conjugate()).real for t in range(k))
                    for k in range(1,9)
                ])
            columns.append(edges)
        for k in range(8):
            invariant.append([sum(edge[k] for edge in group) for group in columns])
            if rotate[r] != r:
                nontrivial.append([g[0][k] + omega*g[1][k] + omega.conjugate()*g[2][k] for g in columns])
    return invariant, nontrivial


def nonsingular(rows):
    assert len(rows) == 45 and all(len(row) == 45 for row in rows)
    exact = acb_mat(rows)
    mid = acb_mat([[exact[i,j].mid() for j in range(45)] for i in range(45)])
    pre = mid.inv()
    defect = pre * exact - acb_mat([[int(i == j) for j in range(45)] for i in range(45)])
    bound = max(sum(abs(defect[i,j]) for j in range(45)).upper() for i in range(45))
    assert bound < 1
    return bound


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    towers = run_path(str(root / "ch39-three-even-positive-towers.py"))
    target_curve = run_path(str(root / "ch39-three-tower-phase-continuation.py"))
    ctx.prec = 320
    height = acb.zeta_zero(towers["NUMBER"]).imag
    groups = towers["orbits"](shell["geo"])
    reps, rotate = site_orbits(shell["geo"])
    assert len(reps) == 7 and sum(rotate[r] == r for r in reps) == 1
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i,i].real for i in range(19))/18
    for q in (2,4,6):
        nums, variables, perm = towers["DATA"][q]
        native = shell["laplacian"](q, height)
        phases = {(i,j): -native[i,j]/abs(native[i,j]) for group in groups for i,j in group}
        raw = native.eig(multiple=True)
        assert len(raw) == 19 and all(v.imag.contains(0) and v.real > 0 for v in raw)
        folded = sorted((z.real*abs(z.real-center) for z in raw), key=float)
        assert all(x < y for x,y in zip(folded, folded[1:]))
        target, _, _ = target_curve["targets"](folded, center, perm)
        assert all(x < y for x,y in zip(target, target[1:]))
        weighted, spectral = newton_box([arb(z)/10**15 for z in nums], variables, groups, phases, target)
        invariant, nontrivial = site_moments(weighted, groups, phases, reps, rotate)
        assert len(invariant) == 56 and len(nontrivial) == 48
        assert len(set(INVARIANT_ROWS)) == 26 and all(k % 8 < 7 for k in INVARIANT_ROWS)
        # 六个自由站轨均取一至七阶，另对前三个轨取第八阶。
        varying = [row for i, row in enumerate(nontrivial) if i < 24 or i % 8 < 7]
        trivial_bound = nonsingular(spectral + [invariant[i] for i in INVARIANT_ROWS])
        varying_bound = nonsingular(varying)
        print("荷", q, "两个 45 阶子式的预乘误差无穷范数上界",
              trivial_bound.str(12), varying_bound.str(12))
    free_reps = [r for r in reps if rotate[r] != r][:3]
    print("六个第八阶观测的零起算站号", [(r, rotate[r]) for r in free_reps])
    print("第七、八阶在精确对称解的固定谱切空间上实秩", 110, 116)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
