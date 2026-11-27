"""在三条极小弧构成的匹配图上，用对称态和反对称态的相位读取引力。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    real_gates = lib("fifth-seventh-shortest-site-gates")
    physics = lib("shortest-common-complex-gate")
    base, cosines, _, basis, station, weights, heat, field = real_gates["geometry"](
        root, fifth, columns, levels
    )
    pairs = ((6, 11), (7, 12), (8, 18))
    diagonal = sum(cosines[1:], cosines[0])
    swap = acb_mat([[sum(basis[a, i].conjugate()*station[b, j]
                         + basis[b, i].conjugate()*station[a, j] for a, b in pairs)
                     for j in range(2)] for i in range(19)])
    plus, minus = (diagonal+swap)/2, (diagonal-swap)/2
    fixed = base-diagonal

    def circle(q):
        denom = 1 + q**2
        return acb((1-q**2)/denom, 2*q/denom)

    def project(point):
        p, m = -circle(point[0]), circle(point[1])
        directions = [acb(0, 2)*w*block/(1+t**2) for t, w, block in
                      zip(point, (p, m), (plus, minus), strict=True)]
        return fixed+p*plus+m*minus, directions

    center = list(map(arb, ("-.36012685842176", ".54926774667943")))
    box, ratio, det = physics["certify"](project, center, arb(".00004"), levels, weights, heat, field)
    assert ratio < arb(".516") and det < -19

    # 匹配的邻接矩阵与图拉普拉斯各自的一参数酉演化均由两个半圆覆盖。
    low, heavy, top = levels
    gain = (heavy-low, top-low)
    normalized = lib("fifth-seventh-cutoff-separation")["moment_errors"](
        field, heat, weights, [".75553"]*19
    )
    assert top < arb(".75553")
    slopes = (sum(alpha for _, alpha, _ in field)/(2*arb.pi()),
              -sum(g for g, _, _ in field)/(8*arb.pi()))

    def gap(p, m):
        proj = fixed + p*plus + m*minus
        amps = [low+sum(k*abs(proj[i, j])**2 for j, k in enumerate(gain)) for i in range(19)]
        gaps = [abs(slope*(sum(w*a**n for w, a in zip(weights, amps, strict=True)) - target)).lower()
                for n, slope, target in zip((2, 4), slopes, normalized, strict=True)]
        if max(gaps) > arb(".025"):
            return max(gaps)
        vals, _ = physics["gravity_system"](proj, [], levels, weights, heat, field, True)
        return max(abs(v).lower() for v in vals.entries())

    for kind in ("邻接", "拉普拉斯"):
        least, covered = arb(100), 0
        for sign in (-1, 1):
            pending = [(arb(-1), arb(1))]
            while pending:
                left, right = pending.pop()
                z = sign*circle(left.union(right))
                common = z.conjugate() if kind == "邻接" else acb(1)
                margin = gap(common, z)
                if margin > arb(".025"):
                    least = min(least, margin)
                    covered += 1
                else:
                    assert right-left > arb(".00001"), (kind, left, right, margin)
                    mid = (left+right)/2
                    pending.extend(((left, mid), (mid, right)))
        assert least > (arb(".0258") if kind == "邻接" else arb("1.16"))
        print(kind, "完整圆的区段数和双矩较大残差的统一下界", covered, least.str(14), flush=True)
    return box


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root, *source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
