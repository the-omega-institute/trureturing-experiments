"""检验允许三个原壳轻能量出现交叉项时共同对称谱作用的最低总次数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def features(energy):
    x, y, z = [(e-5)/5 for e in energy]
    s = x+y+z
    t = x*x + y*y + z*z
    return [s, s*s, t, s**3, s*t, x*y*z]


def verify():
    replay = Path(__file__).resolve().parent
    selector = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = selector["light_spectra"](replay)
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    lowest = run_path(str(replay/"ch39-global-light-action-degree.py"))["lowest"]
    ctx.prec = 256
    adj = pump["adjoint"]
    frames = list(pump["rotation_bases"](shell))
    centers = {-4:2345666, 2:964890}
    reference = {q:features(light(q,n)) for q,n in centers.items()}

    def flat_light(q, factors):
        base, edges = pump["clock_matrices"](shell,q)
        link = edges[0]*factors[0] + edges[1]*factors[1]
        flat = base + link + adj(link)
        return [lowest(adj(b)*flat*b) for b in frames]

    def difference(q, energies):
        return [a-b for a,b in zip(features(energies), reference[q], strict=True)]

    quadratic = [difference(q,light(q,n))[:3] for q,n in
                 ((-4,29491),(2,5566),(2,17371),(2,24118))]
    r2 = selector["interior_radius"](quadratic)
    assert r2 > arb("8.38e-6")

    cubic = [difference(q,light(q,n)) for q,n in ((-4,13155),(-4,17712),(2,36))]
    for i,j in ((16,137),(18,138)):
        clocks = [acb(0,2*arb.pi()*k/160).exp() for k in (i,j)]
        cubic.append(difference(-4,flat_light(-4,clocks)))
    for q,offsets in ((-4,(1,-1)),(2,(0,0))):
        phases = shell["geo"]["phases"](acb.zeta_zero(centers[q]).imag)
        clocks = [phases[d]**q * acb(0,arb(a)/1000).exp()
                  for d,a in zip((724,804),offsets,strict=True)]
        cubic.append(difference(q,flat_light(q,clocks)))
    unit_vectors = []
    for v in cubic:
        length = sum(x*x for x in v).sqrt()
        assert length > 0
        unit_vectors.append([x/length for x in v])
    r3 = selector["interior_radius"](unit_vectors)
    assert r3 > arb("3.77e-6")
    assert sum(light(-4,centers[-4])) - sum(light(2,centers[2])) > arb(".09")
    print("二次有限位移与三次归一位移单形的原点面距下界",r2.str(14),r3.str(14))


if __name__ == "__main__":
    verify()
