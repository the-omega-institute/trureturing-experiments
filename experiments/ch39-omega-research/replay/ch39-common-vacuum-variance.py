"""两面字的任一合格共同全局真空都要求终点方差严格大于 831.05。"""
import json
import sys
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


ORDER = 105
CUTOFF = 370
FLOOR = "831.05"


def phase(t):
    return acb(0, t).exp()


def harmonics(modes, distance, directions, rivals):
    gx, gy, hx, hy = [arb(v)/10**5 for v in directions]
    cuts = [(arb(x)/10**5, arb(y)/10**5, arb(v)/10**5) for x, y, v in rivals]
    assert all(weight > 0 for _, _, weight in cuts)
    mass = sum(weight for _, _, weight in cuts)
    constants, amplitudes = [], []
    for q, (bare, soft) in enumerate(modes, 1):
        rival = sum(weight*sum((z*phase(q*t)).real for z, t in zip(bare, (x, y, x-y), strict=True))
                    for x, y, weight in cuts)
        combined = [(-mass-acb(0, q*g))*z-acb(0, q*h)*w for z, w, g, h in
                    zip(bare, soft, (gx, gy, gx-gy), (hx, hy, hx-hy), strict=True)]
        shift = phase(q*distance)
        constants.append(rival+(combined[2]*shift.conjugate()).real)
        amplitudes.append(combined[0]+combined[1]*shift)
    return constants, amplitudes


def segment(record, modes, distance, characters):
    left, right, directions, rivals, endpoints = record
    assert left < right and len(directions) == 4 and len(endpoints) == 2
    assert all(len(row) == ORDER and min(row) >= 0 for row in endpoints)
    lo, hi = arb(left)/10**8, arb(right)/10**8
    mid, radius = (hi+lo)/2, (hi-lo)/2
    nus = [[arb(v)/10**5 for v in row] for row in endpoints]
    constants, amplitudes = harmonics(modes, distance, directions, rivals)
    cap, curvature = arb(0), arb(0)
    for q, (k, h, n0, n1) in enumerate(zip(constants, amplitudes, *nus, strict=True), 1):
        center = k+(h*phase(-q*mid)).real
        cap += max(abs(center-n0).upper(), abs(center-n1).upper())+radius*q*abs(h)
        curvature += q*q*abs(h)
    tail = arb(CUTOFF**2)/2-cap
    assert tail > arb(FLOOR), (left, right, tail)

    # 光滑部分的二阶导数与整数 j 无关，非负乘子在整段上线性插值。
    count = 1
    while 32*(hi-lo)**2*curvature >= 8*count**2:
        count *= 2
    error = (hi-lo)**2*curvature/(8*count**2)
    columns = []
    for i in range(count+1):
        s = arb(i)/count
        t = lo+s*(hi-lo)
        columns.append([k+(h*phase(-q*t)).real-(1-s)*n0-s*n1 for q, (k, h, n0, n1)
                        in enumerate(zip(constants, amplitudes, *nus, strict=True), 1)])
    values = characters*arb_mat(columns).transpose()
    finite = min((arb(j*j)/2+values[j, i]-error).lower()
                 for j in range(CUTOFF) for i in range(count+1))
    assert finite > arb(FLOOR), (left, right, finite)
    return count, finite, tail.lower()


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    load("common-all-phase-global-order")["verify"](root)
    exact = load("mincut-real-triangle-tomography")
    flavor = load("chiral-neutral-flavor-phase")
    angles = load("mincut-reference-ckm-action")
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    ctx.prec = 224
    theta = exact["unit"]((108469, 109)).arg()
    distance = -angles["reference_angles"](flavor)[0]
    tri = lib["terms"](shell)
    modes = []
    for q in range(1, ORDER+1):
        bare = lib["coeff"](tri, q, arb(0))
        soft = [sum((z**q for z, marked in tri[k] if marked), acb(0))
                for k in ((1, 0), (0, 1), (1, -1))]
        modes.append((bare, soft))
    characters = arb_mat([[(j*q*theta).cos() for q in range(1, ORDER+1)] for j in range(CUTOFF)])
    records = json.loads(Path(__file__).with_suffix(".json").read_text())
    assert records[0][0] == 815000 and records[-1][1] == 1280000
    assert all(a[1] == b[0] for a, b in zip(records, records[1:]))
    results = [segment(rec, modes, distance, characters) for rec in records]
    print("有理闭段、插值小段、逐项整数差上限", len(results), sum(r[0] for r in results), CUTOFF-1)
    print("有限差值和全部尾部的严格下界",
          min(r[1] for r in results).str(18), min(r[2] for r in results).str(18))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
