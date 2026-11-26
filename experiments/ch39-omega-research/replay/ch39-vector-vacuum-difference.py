"""含轻顶和弱矢量的固定五截止联合解，从对称原点到电弱驻点的真空差。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def verify(root):
    joint = run_path(str(root / "ch39-all-field-signed-crossing.py"))
    rotation = run_path(str(root / "ch39-rotation-portal-top-energy.py"))
    signs, reconstructed = rotation["rotation_signs"](root)
    energy = arb("1.4513765e14").union(arb("1.4513786e14"))
    strength = arb(".429954015").union(arb(".429954021"))
    data = joint["joint_system"](root, energy, order=16, last=46)
    matrix = data["matrix"]
    assert not matrix.det().contains(0)
    weights = matrix.solve(data["target"]).entries()
    assert data["target"][0,0] == 1 and all(x[0] == 1 for x in data["cols"])
    abs_mass = sum(abs(c) for c in weights)
    inverse_x = sum(abs(c)/x for c, x in zip(weights, data["sites"], strict=True))
    masses = [z for _, _, z in data["rows"][49:65]]
    assert all(z.overlaps(other) for z, other in zip(masses, reconstructed, strict=True))

    r = (arb("246.22")/arb("1.22089e19"))**2/2
    segment = arb(0).union(r)
    pi = arb.pi()
    g1, g2, yt = map(arb, (".467453", ".63811", ".8616"))
    light = ((-12, yt**2), (6, g2**2/2), (3, (g2**2+3*g1**2/5)/2))
    # 归一化使轻场零端对数的系数精确独立于五权。
    logs = [-d*k**2/(128*pi**2) for d, k in light]
    light_lipschitz = inverse_x*sum(abs(d)*k**3 for d, k in light)/(32*pi**2)
    heavy_lipschitz = arb(0)
    for sign, z in zip(signs, masses, strict=True):
        portal = sign*strength
        floor = (z+portal*segment).lower()
        assert floor > 0
        heavy_lipschitz += abs_mass*abs(portal)**3/(8*pi**2*floor)

    main = 8*pi*r**2*(arb(".129")+sum(logs))
    err_light = 8*pi*r**3*light_lipschitz/6
    err_heavy = 8*pi*r**3*heavy_lipschitz/6
    error = err_light+err_heavy
    drop = main+arb(0, error.upper())
    epsilon = arb("1e-122")
    print("轻顶、W、Z 的对数项系数，以 8πρv² 为单位", *[v.str(18) for v in logs])
    print("W(0)−W(ρv) 的中心、轻场与重层绝对余量、总差与目标倍数",
          main.str(25), err_light.upper().str(12), err_heavy.upper().str(12),
          drop.str(22), (drop/epsilon).str(20))
    assert arb("1.39135134008521247798e-67") < drop < arb("1.39135134008521247801e-67")
    assert err_light < arb("4.927e-89") and err_heavy < arb("1.451e-91")
    assert error < arb("4.95e-89")
    assert epsilon-drop < 0


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
