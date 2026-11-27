"""按三循环给十六层门户分号，同时求规范交汇能标和电弱驻点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def rotation_signs(root, number=2_306_888_100_085, *, height=None):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    if height is None:
        height = acb.zeta_zero(number).imag
    H = shell["laplacian"](3, height)
    sites = shell["geo"]["sites"]
    rotation = [sites.index((j,k,i)) for i,j,k in sites]
    orbits = [{i, rotation[i], rotation[rotation[i]]} for i in range(19)]
    orbits = [o for i,o in enumerate(orbits) if i == min(o)]
    U = acb_mat([[arb(int(i in o))/arb(len(o)).sqrt() for o in orbits] for i in range(19)])
    fixed = [z.real for z in (U.transpose()*H*U).eig(multiple=True)]
    values = [z.real for z in H.eig(multiple=True)]
    delta = arb(23)/25*arb(2)**-20/(height**2 + shell["soft"]["point"]["study"]["offset"])
    center = 15-delta/3
    levels = sorted([(v*abs(v-center),v) for v in values], key=lambda pair:float(pair[0]))
    assert all(a[0] < b[0] for a,b in zip(levels,levels[1:]))
    assert all(sum(v.overlaps(k) for _,v in levels) == 1 for k in fixed)
    signs = [-1 if any(v.overlaps(k) for k in fixed) else 1 for _,v in levels]
    scale = arb("2.435e18")/(724**2*arb("1.22089e19"))
    assert len(fixed) == 7 and signs[3:].count(1) == 10
    return signs[3:], [(scale*m)**2 for m,_ in levels[3:]]


def verify():
    root = Path(__file__).resolve().parent
    joint = run_path(str(root / "ch39-all-field-signed-crossing.py"))
    signs, reconstructed = rotation_signs(root)
    y = arb(".8616")
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    lo, hi = arb("2.2709e14"), arb("2.2716e14")
    small, big = arb(".556760"), arb(".556772")
    for tag, energy in (("下端",lo), ("上端",hi), ("全盒",lo.union(hi))):
        d = joint["joint_system"](root, energy, order=16, last=47)
        mat, sites, heat = d["matrix"], d["sites"], d["heat"]
        assert mat.det() < 0
        weights = mat.solve(d["target"]).entries()
        masses = [z for _,_,z in d["rows"][49:65]]
        assert all(z.overlaps(v) for z,v in zip(masses,reconstructed,strict=True))
        i0 = sum(c*x for c,x in zip(weights,sites,strict=True))
        def diff(z):
            # x(E₂(z/x)−1) = x expm1(−z/x) − z E₁(z/x)。
            return sum(c*(x*(-z/x).expm1()-z*heat(1,z/x)) for c,x in zip(weights,sites,strict=True))
        def j(z):
            return sum(c*heat(1,z/x) for c,x in zip(weights,sites,strict=True))
        def slope(u):
            tail = u*sum(s*diff(z+s*u*rho) for s,z in zip(signs,masses,strict=True))
            return ((4*u-3*y*y)*i0 + tail-3*y*y*diff(y*y*rho))/(8*arb.pi()**2)
        assert slope(small)>0>slope(big)
        if tag != "全盒":
            u = small.union(big)
            curv = (3*y**4*j(y*y*rho)-u**2*sum(j(z+s*u*rho) for s,z in zip(signs,masses,strict=True)))/(8*arb.pi()**2)
            assert curv<arb(".258") if tag=="下端" else curv>arb(".258")
            print(tag, "径向曲率", curv.str(14), flush=True)
        else:
            common = d["initial"][0]+sum(c*col[3] for c,col in zip(weights,d["cols"],strict=True))
            print("整能标盒的五权及统一逆耦合", *[c.str(7) for c in weights],common.str(8), flush=True)
            assert common>0
            floors = joint["running_floors"](d, weights, steps=64)
            assert min(floors)>0
            print("三支沿程下界", *[v.str(9) for v in floors], flush=True)


if __name__ == "__main__":
    verify()
