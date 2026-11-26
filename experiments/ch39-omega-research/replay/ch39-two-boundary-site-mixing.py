"""用同一零点的边界站点和轻谱投影构造夸克混合，固定脉冲后作跨零点检验。"""
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


base = run_path(str(Path(__file__).with_name("ch39-moving-triality-link.py")))
edge = base["study"]["edge"]
shell = base["shell"]
omega = acb(-1, arb(3).sqrt()) / 2


def light_modes(number):
    height = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i, i].real for i in range(19)) / 18
    return edge["mode_columns"](edge["filters"]["spectra"](1, height, center))


def boundary_frames(modes):
    sites = shell["geo"]["sites"]
    start = shell["boundary_sites"]()[0]
    x, y, z = sites[start]
    orbit = (start, sites.index((y, z, x)), sites.index((z, x, y)))
    assert set(orbit) == set(shell["boundary_sites"]())
    anchor = [modes[start, k] for k in range(3)]
    assert all(abs(v)**2 > arb("3e-6") for v in anchor)
    # 用同一个基点固定三条谱线的相位。
    frames = [[modes[site, k] * anchor[k].conjugate() / abs(anchor[k])
               for k in range(3)] for site in orbit]
    for step, frame in enumerate(frames):
        assert all((frame[k] - abs(anchor[k]) * omega**(k*step)).contains(0)
                   for k in range(3))
    return frames


def pulse(i, j, frame, strength):
    amplitude = frame[i].conjugate() * frame[j]
    parameter = strength * abs(amplitude)
    denom = 1 + parameter**2
    cosine, sine = (1-parameter**2)/denom, 2*parameter/denom
    phase = amplitude / abs(amplitude)
    out = acb_mat([[int(a == b) for b in range(3)] for a in range(3)])
    out[i, i] = out[j, j] = cosine
    out[i, j] = sine * phase
    out[j, i] = -sine * phase.conjugate()
    return out


def calibration(frames):
    th12, th23, th13, delta = map(arb, (".22704", ".04275", ".00377", "1.139"))
    phase = acb(0, -delta).exp()
    s, c = th13.sin_cos()
    target = acb_mat([[c, s*phase], [-s*phase.conjugate(), c]])
    left, d, right = base["cartan_angles"](target)
    axis = omega**2
    angle = (d.sin() / abs(axis.imag)).asin()
    ss, cc = angle.sin_cos()
    middle = acb_mat([[cc, ss*axis], [-ss*axis.conjugate(), cc]])
    a, e, b = base["cartan_angles"](middle)
    assert abs(e-d) < arb("1e-40")
    # 矩阵按从左到右的乘法次序列出。位置 1 指定第二个站点。
    turns = [(1,2,0,th23), (0,2,0,left-a), (0,2,1,angle),
             (0,2,0,right-b), (0,1,0,th12)]
    controls = [(i,j,site,(theta/2).tan() / abs(frames[site][i]*frames[site][j]))
                for i,j,site,theta in turns]
    real = [acb(1)]*3
    target3 = acb_mat([[c,0,s*phase], [0,1,0], [-s*phase.conjugate(),0,c]])
    full = pulse(1,2,real,(th23/2).tan()) * target3 * pulse(0,1,real,(th12/2).tan())
    return controls, full


def read(frames, controls, orientation):
    out = acb_mat([[int(a == b) for b in range(3)] for a in range(3)])
    for i, j, site, coefficient in controls:
        out = out * pulse(i, j, frames[orientation*site], coefficient)
    return out


def observables(v):
    jarlskog = (v[0,1]*v[1,2]*v[0,2].conjugate()*v[1,1].conjugate()).imag
    return [abs(v[0,1]), abs(v[1,2]), abs(v[0,2]), jarlskog]


def boundary_orientation():
    sites = shell["geo"]["sites"]
    a = shell["boundary_sites"]()[0]
    x, y, z = sites[a]
    b, c = sites.index((y,z,x)), sites.index((z,x,y))
    bonds = {}
    for i, j, chord, link in shell["geo"]["edges"]:
        if i in (a,b,c) and j in (a,b,c):
            assert chord == 200
            bonds[i,j] = link
            bonds[j,i] = link.conjugate()
    flux = (bonds[a,b]*bonds[b,c]*bonds[c,a]).imag
    assert arb("-.052457") < flux < arb("-.052456")
    assert abs(flux**2 - arb(1500625)/545338513) < arb("1e-45")
    return flux


def endpoint_floor(modes, theta23, theta12):
    sites = shell["geo"]["sites"]
    pole = sites.index((17,17,17))
    # 在极点，两条非平凡轮换特征线分量恒为零。
    assert all(modes[pole,k].contains(0) for k in (1,2))
    sites_with_mixing = [s for s in range(19) if s != pole]
    assert all(abs(modes[s,k]) > 0 for s in sites_with_mixing for k in range(3))
    small23, small12 = (theta23/2).tan(), (theta12/2).tan()
    branches = list(product((small23, 1/small23), (small12, 1/small12)))
    bounds = []
    for i,j,k in permutations(range(3)):
        left = [abs(modes[s,j]*modes[s,k]) for s in sites_with_mixing]
        right = [abs(modes[s,i]*modes[s,j]) for s in sites_with_mixing]
        for first, last in product(left, right):
            for a, b in branches:
                ratio = a*last / (b*first)
                bounds.append(max(ratio.lower(), (1/ratio).lower()))
    return min(bounds)


def verify():
    ctx.prec = 448
    print("有向边界面和乐的虚部", boundary_orientation().str(12))
    spectra = {n: light_modes(n) for n in (11661,11934,12207,22027,28325,48067,58214)}
    controls, target = calibration(boundary_frames(spectra[28325]))
    for number in (28325,48067,58214):
        frames = boundary_frames(spectra[number])
        direct = read(frames, controls, 1)
        reverse = read(frames, controls, -1)
        assert max(abs(v) for v in (direct.conjugate()-reverse).entries()) < arb("1e-38")
        if number == 28325:
            assert max(abs(v) for v in (direct-target).entries()) < arb("1e-38")
        vals = observables(direct)
        print(number, "us, cb, ub, J", *(v.str(13) for v in vals))
        assert vals[-1] > 0 and observables(reverse)[-1] < 0

    wide23 = arb(".0415").union(arb(".0440"))
    wide12 = arb(".225").union(arb(".229"))
    for number, modes in spectra.items():
        exact = endpoint_floor(modes, arb(".04275"), arb(".22704"))
        wide = endpoint_floor(modes, wide23, wide12)
        print(number, "端点强度倍率下界，中心值与角度窗", exact.str(9), wide.str(9))
        assert exact > arb("1.0956") and wide > arb("1.0548")


if __name__ == "__main__":
    verify()
