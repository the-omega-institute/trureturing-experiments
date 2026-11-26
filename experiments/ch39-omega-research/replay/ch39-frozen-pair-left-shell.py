"""在未参与两点共同拟合的左邻，逐层应用已固定的谱函数与夸克味空间。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify(root):
    training = run_path(str(root / "ch39-two-zero-simultaneous-quark-gravity.py"))
    frame, nodes = training["verify"](root)
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    gravity = run_path(str(root / "ch39-three-index-free-atoms.py"))
    spectrum = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    weak = run_path(str(root / "ch39-full-neutral-free-cutoff.py"))["weak_column"]
    left = 2_306_888_100_084
    moments, heat, fields = gravity["broken_fields"](root, (left,))
    ctx.prec = 320
    height = acb.zeta_zero(left).imag

    def response(level):
        weights = [1 / (level - t)**2 for t, _ in nodes]
        assert all(weight > 0 for weight in weights)
        return sum(weight * value for weight, (_, value) in zip(weights, nodes, strict=True)) / sum(weights)

    mass_squared = []
    top_squared = None
    for charge, projectors in zip((-4, 2), native["native"](root, left), strict=True):
        amplitudes = [response(level) for level in ordered(shell, charge, height)]
        four = sum((v * p for v, p in zip(amplitudes, projectors, strict=True)), acb_mat(4, 4))
        eigenvalues = sorted((v.real for v in four.eig(algorithm="rump")), key=lambda v: v.mid())
        assert all(x < y for x, y in zip(eigenvalues, eigenvalues[1:]))
        print("荷、任意三平面压缩后最大本征值的锐范围两端", charge,
              eigenvalues[-2].str(13), eigenvalues[-1].str(13), flush=True)
        if charge == -4:
            assert eigenvalues[-1] < arb(".252")
            top_squared = eigenvalues[-2].union(eigenvalues[-1])**2
        else:
            assert eigenvalues[-2] > arb(".1584")
        matrix = native["adj"](frame) * four * frame
        mass_squared.append(matrix * matrix)
    ckm, spectra, _ = flavor["mixing"](*mass_squared)
    for levels in spectra:
        print("左邻夸克有效汤川值", *[y.sqrt().str(12) for y in levels], flush=True)
    print("左邻 s12、s23、s13、有符号 J", *[x.str(12) for x in ckm], flush=True)
    assert ckm[2] > arb(".89")

    _, energies, pivots, _ = next(row for row in spectrum["shell_weights"](root, left) if row[0] == 6)
    cutoffs = [response(t) for t in energies]
    planck, mz = arb("1.22089e19"), arb("91.1876")
    assert all(mz / planck < x < arb(".948") for x in cutoffs)
    p0, p1 = (mz / planck)**2, (arb("1.176397e14") / planck)**2
    columns = [moments(fields[0], heat, x*x) for x in cutoffs]
    k, w = [sum(c.real * row[j] for c, row in zip(pivots[0], columns, strict=True)) for j in (0, 1)]
    diff = sum(c.real * weak(x*x, p0, p1) for c, x in zip(pivots[0], cutoffs, strict=True))
    print("左邻截止范围（GeV）", (min(x.lower() for x in cutoffs)*planck).str(12),
          (max(x.upper() for x in cutoffs)*planck).str(12), flush=True)
    print("左邻 K、W、D", k.str(12), w.str(12), diff.str(12), flush=True)
    assert k < arb("-.4357") and w > arb(".9909") and diff < arb("-.07425")
    # 对任意味空间，用 E1(s) <= 1+max(0,-log(s)) 包围轻顶质量替换的导数。
    scale = (arb("246.22") / planck)**2 / 2
    g, coefficient, old_z = fields[0][1]
    assert g == -12 and coefficient == 1 and old_z.overlaps(arb(".8616")**2 * scale)
    new_z = top_squared * scale
    floor_z = old_z.union(new_z).lower()
    assert 0 < floor_z < 1
    variation = sum(abs(c.real) for c in pivots[0])
    error = abs(old_z-new_z) * variation * (1-floor_z.log()) / (2*arb.pi())
    assert error < arb("1e-30") and k + error < arb("-.4357")
    print("任意复味空间下轻顶质量反馈绝对值的上界", error.upper().str(12), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
