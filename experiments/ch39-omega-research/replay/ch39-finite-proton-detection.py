"""在两组指定的实际正零点中探测任一合格真空的树级质子因子。"""
from itertools import product
import json
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def compare_errors(records, errors, budget, modes, distance, harmonics):
    assert records[0][0] == 815000 and records[-1][1] == 1280000
    assert all(a[1] == b[0] for a, b in zip(records, records[1:]))
    upper = arb(0)
    for left, right, directions, rivals in records:
        assert left < right
        lo, hi = arb(left)/10**8, arb(right)/10**8
        # harmonics 的梯度与比较乘子单位为 10^-5，本表为 10^-8。
        constants, amplitudes = harmonics(modes, distance, [arb(v)/1000 for v in directions],
                                          [(x, y, arb(w)/1000) for x, y, w in rivals])
        for q, (e, k, h) in enumerate(zip(errors, constants, amplitudes, strict=True), 1):
            pending = [(lo, hi)]
            while pending:
                l, r = pending.pop()
                mid, radius = (l+r)/2, (r-l)/2
                z = h*acb(0, -q*mid).exp()
                high = e-k-z.real+radius*q*abs(z.imag)+(radius*q)**2*abs(h)/2
                if high < budget:
                    upper = max(upper, high.upper())
                else:
                    assert r-l > arb("1e-14"), (left, right, q, high)
                    pending.extend(((l, mid), (mid, r)))
    return upper


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    floor = load("proton-all-vacua-floor")
    floor["verify"](root)
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    angle = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    geometry = load("mincut-real-triangle-tomography")
    harmonics = load("common-vacuum-variance")["harmonics"]
    ctx.prec = 192
    clocks = floor["CLOCKS"]
    bounds = [arb(".18546083"), arb(".69944252"), arb(".59620304")]
    shares = [arb(".019"), arb(".015"), arb(".072")]
    scales = [[arb(s)/1000/b/w for s, w in zip(row, shares, strict=True)]
              for row, b in zip(floor["SEPARATORS"], bounds, strict=True)]
    tri = lib["terms"](shell)
    distance = -angle["reference_angles"](flavor)[0]
    modes = [(lib["coeff"](tri, q, arb(0)),
              [sum((z**q for z, marked in tri[k] if marked), acb(0)) for k in clocks])
             for q in range(1, 106)]
    labels = {k: sorted({pair for clock, pair in geometry["enumerate_faces"]().values() if clock == k})
              for k in clocks}
    data = json.loads(Path(__file__).with_suffix(".json").read_text())
    indices = [[row[0] for row in pack["readouts"]] for pack in data]
    assert all(numbers == sorted(set(numbers)) for numbers in indices)
    assert set(indices[0]) < set(indices[1]) and list(map(len, indices)) == [218, 245]
    assert min(indices[1]) == 33 and max(indices[1]) == 3965
    # 直接按编号生成实际高度；两个证书复用同一套相位。
    geo = shell["geo"]
    phases = {}
    for n in indices[1]:
        height = acb.zeta_zero(n).imag
        assert height > 14
        ph = geo["phases"](height)
        phases[n] = (ph[724], ph[804], ph[724]*ph[804].conjugate())
    assert shell["soft"]["point"]["study"]["offset"] >= 0
    delta = arb(23)/25*arb(2)**-20/14**2
    corners = list(product((0, 1), repeat=3))
    results = []
    for pack, numbers, limit, claim in zip(data, indices, (".87", ".69"), (".0001633", ".0003645"), strict=True):
        filters = [[arb(v)/10**6 for v in row[1:]] for row in pack["readouts"]]
        assert all(len(row)==3 for row in filters)
        costs = [sum(abs(sum(x*t for x, t in zip(row, vertex, strict=True))) for row in filters)
                 for vertex in corners]
        norm = max(c.upper() for c in costs)
        errors = [arb(v)/10**5 for v in pack["error_bounds"]]
        assert len(errors) == 105 and min(errors) >= 0
        weights = arb_mat(filters)
        for k, clock in enumerate(clocks):
            # 先合并零点上的三列实滤波器，最后分别乘每个面的精确相位。
            powers = [[phases[n][k]**q for n in numbers] for q in range(1, 106)]
            real = arb_mat([[z.real for z in row] for row in powers])*weights
            imag = arb_mat([[z.imag for z in row] for row in powers])*weights
            for pair in labels[clock]:
                u = geometry["unit"](pair)
                for q in range(1, 106):
                    uq = u**q
                    residual = []
                    for j, active in enumerate((q <= 5, q == 104, q == 105)):
                        target = scales[k][j] if active else 0
                        v = ((target-real[q-1, j])*uq.real+imag[q-1, j]*uq.imag).upper()
                        residual.append(max(arb(0), v))
                    assert sum(residual) < errors[q-1], (k, pair, q)
        upper = compare_errors(pack["comparisons"], errors, arb(limit), modes, distance, harmonics)
        bound = (1-upper)/norm-delta
        assert bound > arb(claim)
        print("零点枚数、编号范围、比较闭段、绕数总误差、滤波器范数、实际因子绝对值下界",
              len(numbers), (min(numbers),max(numbers)), len(pack["comparisons"]),
              upper.str(14), norm.str(14), bound.str(17))
        results.append((numbers, filters, upper, delta))
    return results


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
