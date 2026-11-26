"""固定能标与五个截止，在同一电弱背景拟合两项引力、未破缺群关系及驻点曲率。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def broken_system(root, energy, portal_box):
    joint = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    rotation = run_path(str(root/"ch39-rotation-portal-top-energy.py"))
    signs, reconstructed = rotation["rotation_signs"](root)
    pi = arb.pi()
    r = (arb("246.22")/arb("1.22089e19"))**2/2
    yt, g1, g2 = map(arb, (".8616", ".467453", ".63811"))
    top, zw, zz = yt**2, g2**2/2, (g2**2+3*g1**2/5)/2
    light = ((-12,arb(1),top), (6,-arb(1),zw), (3,-arb(1)/2,zz))
    data = joint["joint_system"](root, energy, order=16, last=47)
    heat, sites, rows = data["heat"], data["sites"], data["rows"]
    masses = [z for _, _, z in rows[49:65]]
    assert all(z.overlaps(v) for z, v in zip(masses, reconstructed, strict=True))

    def em(b):
        return 5*b[0]/3+b[1]

    def channels(photon, colour):
        return (photon, colour, 8*colour/3)

    # 三支数组存放 A_em、A_3、8A_3/3，便于复用单调减法核的路径认证。
    photon, colour = em(data["initial"]), data["initial"][2]
    fermion, scalar, vector, gauge = range(4)
    b_light = data["species"][0][2]
    top_em, top_colour = arb(16)/9, arb(2)/3
    species = [(fermion,arb(0),channels(em(b_light)-top_em,b_light[2]-top_colour)),
               (fermion,top*r,channels(top_em,top_colour)),
               (vector,zw*r,channels(-7,arb(0))),
               (gauge,arb(0),channels(0,arb(-11)))]
    # 重门户只平移十六个标量双重态，每层的净电磁 beta 系数为 1/3。
    nscalar = 0
    for kind, z, b in data["species"][3:]:
        if kind == scalar:
            assert z.overlaps(masses[nscalar]) and (em(b)-arb(1)/3).contains(0) and b[2] == 0
            z += signs[nscalar]*portal_box*r
            nscalar += 1
        species.append((kind,z,channels(em(b),b[2])))
    assert nscalar == 16
    rule = joint["angular_rule"](16,47)
    p0, p1 = data["p0"], data["p1"]
    log = (p1/p0).log()
    per_site = []
    for x in sites:
        de, dc = arb(0), arb(0)
        for kind, z, (be, bc, _) in species:
            integral = joint["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*pi)
            de -= be*integral
            dc -= bc*integral
        per_site.append((de,dc))

    # 最小曲率耦合，Feynman–'t Hooft 规范。每个破缺生成元合计 (3,-1/2)。
    g0, a0, _ = rows[0]
    g0 -= sum(g for g, _, _ in light)
    a0 -= sum(a for _, a, _ in light)
    assert g0 == -59 and (a0-arb(2)/3).contains(0)
    static = [(g0,a0,arb(0)), *rows[1:49], *rows[65:],
              *((g,a,k*r) for g,a,k in light)]
    for g, a, z in rows[49:65]:
        static.append((g-4,a-arb(2)/3,z))
    moments = run_path(str(root/"ch39-signed-cutoff-gravity.py"))["moments"]
    target = arb_mat([[1], [1], [arb("1e-122")], [8*colour/3-photon], [arb(".258")]])
    def fixed_portal(u):
        gravity = static+[(4,arb(2)/3,z+sg*u*r) for z, sg in zip(masses, signs, strict=True)]
        assert sum(g for g, _, _ in gravity) == sum(g for g, _, _ in rows)
        assert (sum(a for _, a, _ in gravity)-sum(a for _, a, _ in rows)).contains(0)
        cols, slopes = [], []
        for x, (de, dc) in zip(sites, per_site, strict=True):
            k, w = moments(gravity,heat,x)
            # 同一列内先合并斜率，保留轻重场的抵消。
            s = sum(g*h*x*heat(2,h*r/x)/4 for g, _, h in light)
            hess = -sum(g*h**2*heat(1,h*r/x)/4 for g, _, h in light)
            for sg, z in zip(signs,masses,strict=True):
                mass = z+sg*u*r
                s += sg*u*x*heat(2,mass/x)
                hess -= u**2*heat(1,mass/x)
            cols.append([1,k,w,de-8*dc/3,hess/(8*pi**2)])
            slopes.append(s/(8*pi**2))
        matrix = arb_mat(cols).transpose()
        det = matrix.det()
        assert not det.contains(0)
        c = matrix.solve(target).entries()
        derivative = sum(v*s for v, s in zip(c, slopes, strict=True))
        return c, derivative, det

    path = dict(data, species=species, initial=channels(photon,colour))
    return dict(original=data, path=path, weights_at=fixed_portal, per_site=per_site,
                masses=masses, signs=signs, rho=r, weak_mass=zw*r)


def verify(root):
    low, high = arb(".4299540174"), arb(".4299540176")
    portal_box = low.union(high)
    center = arb("1.45137757e14")
    energy = (center-arb("1e5")).union(center+arb("1e5"))
    system = broken_system(root, energy, portal_box)
    path, per_site, fixed_portal = (system[k] for k in ("path","per_site","weights_at"))
    photon, colour, _ = path["initial"]
    _, left, _ = fixed_portal(low)
    _, right, _ = fixed_portal(high)
    assert left > 0 > right
    c, _, det = fixed_portal(portal_box)
    assert -14 < det < -12
    assert all(sign*v > 0 for sign, v in zip((1,-1,1,-1,1),c,strict=True))
    negative = -c[1]-c[3]
    assert arb("11.5104") < negative < arb("11.5106")
    a_em = photon+sum(v*d[0] for v,d in zip(c,per_site,strict=True))
    a_3 = colour+sum(v*d[1] for v,d in zip(c,per_site,strict=True))
    assert arb("38.17") < a_3 < arb("38.18") and arb("101.79") < a_em < arb("101.80")
    print("两端驻点残差、五个截止权、负质量、能标处电磁与色逆耦合",
          left.str(15), right.str(15), *[v.str(14) for v in c], negative.str(12),
          a_em.str(12), a_3.str(12), flush=True)

    joint = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    floors = joint["running_floors"](path,c,steps=144)
    assert min(floors) > 0
    print("M_Z 至匹配能标之间电磁、色逆耦合下界", *[v.str(12) for v in floors[:2]])


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
