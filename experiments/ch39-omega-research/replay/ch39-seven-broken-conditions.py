"""同一电弱背景的五权、共同门户和匹配能标，按对角有限动量核解七式。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def mixed_residual(system, weights, portal, moment, scalar_displacement):
    data = system["original"]
    _, a2, a3 = data["initial"]
    # 使用精确归一化，将色弱差的共同部分移出五权求和。
    columns = [col[4]-col[5] for col in data["cols"]]
    pivot = columns[2].mid()
    symmetric = a2-a3+pivot+sum(c*(d-pivot) for c, d in zip(weights, columns, strict=True))

    zw = system["weak_mass"]
    r0, r1 = zw/data["p0"], zw/data["p1"]
    high = moment(r1.lower()).union(moment(r1.upper()))
    error = (arb(43)/6)*sum(abs(c)*zw/x for c, x in zip(weights, data["sites"], strict=True))
    gauge = (high-moment(r0)+arb(0,error.upper()))/(4*arb.pi())
    portal_shift = scalar_displacement(system["masses"], system["signs"], portal*system["rho"],
                                       data, weights, panels=32)/2
    return symmetric+gauge+portal_shift


def verify(root):
    ctx.prec = 192
    build = run_path(str(root/"ch39-six-broken-conditions.py"))["broken_system"]
    mixed = run_path(str(root/"ch39-mixed-neutral.py"))
    scalar = run_path(str(root/"ch39-em-colour.py"))["scalar_displacement"]
    mixed["charge_algebra"]()
    lower, upper = map(arb, ("1.1763969e14", "1.17639745e14"))
    small, big = map(arb, (".42995401871", ".42995401873"))
    portal = small.union(big)
    faces = []
    for energy in (lower, upper):
        system = build(root, energy, portal)
        weights, _, _ = system["weights_at"](portal)
        residual = mixed_residual(system, weights, portal, mixed["mixed_log_moment"], scalar)
        print("能标面、混合残差的下界与上界", energy.str(17), residual.lower().str(14),
              residual.upper().str(14), flush=True)
        faces.append(residual)
    assert faces[0] > arb("5e-7") and faces[1] < arb("-5e-7")

    system = build(root, lower.union(upper), portal)
    candidate = system["weights_at"]
    _, below, _ = candidate(small)
    _, above, _ = candidate(big)
    weights, _, det = candidate(portal)
    print("门户面的径向导数", below.str(14), above.str(14), flush=True)
    print("全盒行列式及五权", det.str(14), *[c.str(14) for c in weights], flush=True)
    assert below > arb("1e-13") and above < arb("-1e-13")
    assert arb("-12.226") < det < arb("-12.224")
    initial = system["path"]["initial"]
    shifts = system["per_site"]
    endpoints = [initial[i]+sum(c*d[i] for c, d in zip(weights, shifts, strict=True))
                 for i in range(2)]
    assert all(sign*c > 0 for sign, c in zip((1,-1,1,-1,1), weights, strict=True))
    negative = -weights[1]-weights[3]
    print("负权总量、终点电磁及色逆耦合", negative.str(13),
          *[v.str(13) for v in endpoints], flush=True)
    assert arb("11.7315") < negative < arb("11.7317")
    assert arb("104.076") < endpoints[0] < arb("104.084")
    assert arb("39.027") < endpoints[1] < arb("39.033")
    path_bounds = run_path(str(root/"ch39-all-field-signed-crossing.py"))["running_floors"]
    floors = path_bounds(system["path"], weights, steps=144)
    print("电磁、色沿程下界", *[v.str(13) for v in floors[:2]], flush=True)
    assert floors[0] > arb("96.315") and floors[1] > arb("3.167")


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
