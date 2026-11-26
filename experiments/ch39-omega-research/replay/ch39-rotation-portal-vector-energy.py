"""十六层门户仅按三循环分号时，含轻弱矢量的电弱驻点与三规范交汇。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def verify(root):
    joint = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    rotation = run_path(str(root/"ch39-rotation-portal-top-energy.py"))
    signs, reconstructed = rotation["rotation_signs"](root)
    assert [j for j, s in enumerate(signs, 1) if s < 0] == [3, 5, 7, 8, 14, 16]
    lo, hi = arb("1.4513765e14"), arb("1.4513786e14")
    low_portal, high_portal = arb(".429954015"), arb(".429954021")
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    g1, g2, yt = map(arb, (".467453", ".63811", ".8616"))
    light = ((-12, yt**2), (6, g2**2/2), (3, (g2**2+3*g1**2/5)/2))

    for tag, energy in (("下端", lo), ("上端", hi), ("整盒", lo.union(hi))):
        data = joint["joint_system"](root, energy, order=16, last=46)
        mat = data["matrix"]
        assert not mat.det().contains(0)
        weights = mat.solve(data["target"]).entries()
        sites, heat = data["sites"], data["heat"]
        masses = [z for _, _, z in data["rows"][49:65]]
        assert all(z.overlaps(v) for z, v in zip(masses, reconstructed, strict=True))

        def potential(portal):
            slope, curvature = arb(0), arb(0)
            for x, c in zip(sites, weights, strict=True):
                # 先在每个节点合并各场斜率，以保留轻重层的抵消。
                first = sum(n*k*x*heat(2,k*rho/x)/4 for n,k in light)
                second = -sum(n*k**2*heat(1,k*rho/x)/4 for n,k in light)
                for sign, z in zip(signs, masses, strict=True):
                    shifted = z+sign*portal*rho
                    first += sign*portal*x*heat(2,shifted/x)
                    second -= portal**2*heat(1,shifted/x)
                slope += c*first
                curvature += c*second
            return slope/(8*arb.pi()**2), curvature/(8*arb.pi()**2)

        if tag == "整盒":
            assert potential(low_portal)[0] > 0 > potential(high_portal)[0]
            unified = data["initial"][0]+sum(c*col[3] for c, col in zip(weights, data["cols"], strict=True))
            assert arb("38.24") < unified < arb("38.36")
            assert all(s*c > 0 for s, c in zip((1,-1,1,-1,1), weights, strict=True))
            negative_mass = -weights[1]-weights[3]
            assert arb("11.56") < negative_mass < arb("11.57")
            print("整盒五权、负质量、统一逆耦合", *[w.str(10) for w in weights],
                  negative_mass.str(10), unified.str(10), flush=True)
            floors = joint["running_floors"](data, weights, steps=144)
            assert min(floors) > 1
            print("三支逆耦合全路径下界", *[v.str(12) for v in floors], flush=True)
        else:
            residual = potential(low_portal.union(high_portal))[1]-arb(".258")
            assert residual < 0 if tag == "下端" else residual > 0
            print(tag, "曲率残差", residual.str(14), flush=True)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
