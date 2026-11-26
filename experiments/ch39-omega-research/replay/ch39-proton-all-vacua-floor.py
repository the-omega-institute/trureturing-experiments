"""全部合格共同真空与 96 面任意非负概率的树级质子因子均方根下界。"""
import json
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


CLOCKS = ((1, 0), (0, 1), (1, -1))
SEPARATORS = ((335, 874, -1000), (861, -1000, 801), (206, 1000, 111))


def phase(t):
    return acb(0, t).exp()


def share_bound(records, target, threshold, modes, distance):
    assert records[0][0] == 815000 and records[-1][1] == 1280000
    assert all(a[1] == b[0] for a, b in zip(records, records[1:]))
    worst, leaves = arb(1), 0
    for left, right, directions, comparisons in records:
        assert left < right
        lo, hi = arb(left)/10**8, arb(right)/10**8
        gx, gy, hx, hy = [arb(v)/10**9 for v in directions]
        rivals = [(arb(x)/10**5, arb(y)/10**5, arb(w)/10**9) for x, y, w in comparisons]
        assert all(w > 0 for _, _, w in rivals)
        mass = sum(w for _, _, w in rivals)
        for q, (bare, soft) in enumerate(modes, 1):
            constant = arb(int(q in target)) + sum(w*sum((z*phase(q*t)).real
                       for z, t in zip(bare, (x, y, x-y), strict=True)) for x, y, w in rivals)
            combined = [(-mass-acb(0, q*g))*z-acb(0, q*h)*v for z, v, g, h in
                        zip(bare, soft, (gx, gy, gx-gy), (hx, hy, hx-hy), strict=True)]
            shift = phase(q*distance)
            constant += (combined[2]*shift.conjugate()).real
            harmonic = combined[0]+combined[1]*shift
            pending = [(lo, hi)]
            while pending:
                l, r = pending.pop()
                mid, radius = (l+r)/2, (r-l)/2
                z = harmonic*phase(-q*mid)
                lower = constant+z.real-radius*q*abs(z.imag)-(radius*q)**2*abs(harmonic)/2
                if lower > threshold:
                    leaves += 1
                    worst = min(worst, lower.lower())
                else:
                    assert r-l > arb("1e-15"), (left, right, q, lower)
                    pending.extend(((l, mid), (mid, r)))
    return worst, leaves


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    # 核对全角域归约和已有的 w_105 > 0.072 证书。
    load("all-native-spectral-masses")["verify"](root)
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    angles = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    geometry = load("mincut-real-triangle-tomography")
    ctx.prec = 224
    distance = -angles["reference_angles"](flavor)[0]
    tri = lib["terms"](shell)
    modes = [(lib["coeff"](tri, q, arb(0)),
              [sum((z**q for z, marked in tri[k] if marked), acb(0)) for k in CLOCKS])
             for q in range(1, 106)]
    records = json.loads(Path(__file__).with_suffix(".json").read_text())
    low, high, last = arb(".019"), arb(".015"), arb(".072")
    for name, target, threshold in (("104", {104}, high), ("low5", set(range(1, 6)), low)):
        bound, leaves = share_bound(records[name], target, threshold, modes, distance)
        print("绕数、覆盖记录、小段数、认证下界", name, len(records[name]), leaves, bound.str(15))

    faces = list(geometry["enumerate_faces"]().values())
    norms = []
    bounds = (".18546083", ".69944252", ".59620304")
    for clock, coefficients, expected, bound in zip(CLOCKS, SEPARATORS, (42, 42, 12), bounds, strict=True):
        labels = [pair for k, pair in faces if k == clock]
        assert len(labels) == expected
        phases = [geometry["unit"](pair) for pair in sorted(set(labels))]
        alpha, beta, gamma = [arb(v)/1000 for v in coefficients]
        minimum = min((alpha*(u**q).real+beta*(u**104).real+gamma*(u**105).real).lower()
                      for u in phases for q in range(1, 6))
        assert minimum > arb(bound)
        norm = arb(bound)/(5*(alpha/low)**2+(beta/high)**2+(gamma/last)**2).sqrt()
        norms.append(norm)
        print("钟类、逐面分离余量、类内傅里叶范数下界", clock, minimum.str(15), norm.str(15))
    combined = (2*sum(k**-2 for k in norms))**-1
    axis = (2*sum(k**-2 for k in norms[:2]))**-1
    assert combined.sqrt() > arb(".0016054") and axis.sqrt() > arb(".0016674")
    print("全部面概率与仅两轴面概率的均方根下界", combined.sqrt().str(16), axis.sqrt().str(16))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
