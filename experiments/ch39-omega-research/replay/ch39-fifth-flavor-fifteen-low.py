"""第五壳共同第四步精确配平完整夸克读数的十五低层见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


NUMBER = 2_306_888_100_090
QUATERNION = (-11396407304, -13749932891, -73029510396, 997169871114)
DIRECTIONS = ((3, 0, "r"), (14, 9, "r"), (17, 6, "r"), (7, 5, "i"),
              (9, 7, "r"), (10, 9, "i"), (9, 8, "r"), (11, 1, "i"),
              (9, 8, "i"), (15, 7, "r"), (15, 14, "i"), (17, 9, "r"))


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    point, whole, deriv, widths, upstream, incoming, levels = lib("fourth-third-pivot-joint")["verify"](
        root, return_chart=True)
    quarks = lib("fifth-two-charge-full-quark")
    rotation = quarks["rotation"]
    turn = rotation(QUATERNION) * rotation(quarks["QUATERNION"])
    shape = lib("independent-unitary-quark-right-holdout")["coordinates"]
    make = lib("next-coherent-full-ckm")["block_from_coordinates"]
    projectors = lib("free-four-site-native-positive-spectrum")["native"](root, NUMBER)
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    shell = lib("three-ray-chiral-boundary")
    ctx.prec = 384
    left = (adj(whole**2), whole)
    mid = (adj(point**2), point)
    changes = ([adj(e*whole + whole*e) for e in deriv], deriv)
    sectors, mid_sources, sources, tangents = [], [], [], []
    for i, (charge, (src, _, w, _, goal), rows) in enumerate(zip((-4, 2), upstream, projectors, strict=True)):
        target = turn * make(goal.entries()) * adj(turn)
        sector = (adj(spectrum(shell, charge, NUMBER)[1]), w, arb_mat([shape(r) for r in rows]).transpose(),
                  arb_mat([[v] for v in shape(target)]))
        sectors.append(sector)
        mid_sources.append(mid[i]*src)
        sources.append(left[i]*src)
        tangents.append([e*src for e in changes[i]])
    ctx.prec = 384
    eye = lib("rank-two-generator-up-quark")["eye"]
    identity = eye(19)
    coupling = lib("next-common-ckm")
    local = coupling["local"]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"](
        root, eye, "ch39-fifth-flavor-fifteen-low-cayley.tsv")

    def system(coords, srcs, *, cross=False):
        z, dz = local(seed, coords, DIRECTIONS, identity, not cross)
        values, jac = [], []
        for i, ((dst, ws, cols, goal), s, op) in enumerate(zip(sectors, srcs, (adj(z*z), z), strict=True)):
            ov = dst*op*s
            output = arb_mat([[sum(w*abs(ov[k, j])**2 for j, w in enumerate(ws))] for k in range(19)])
            values.extend((cols*output - goal).entries())
            if cross:
                dv = [dst*op*ds for ds in tangents[i]]
            else:
                ops = [adj(e*z + z*e) for e in dz] if i == 0 else dz
                dv = [dst*e*s for e in ops]
            rows = []
            for do in dv:
                dw = arb_mat([[2*sum(w*(ov[k, j].conjugate()*do[k, j]).real for j, w in enumerate(ws))]
                              for k in range(19)])
                rows.append((cols*dw).entries())
            jac.extend(arb_mat(rows).transpose().tolist())
        return arb_mat([[v] for v in values]), arb_mat(jac)

    center = [arb(0)]*12
    for _ in range(3):
        f, j = system(center, mid_sources)
        step = j.mid().solve(f.mid()).entries()
        center = [(x-d).mid() for x, d in zip(center, step, strict=True)]
    f, j = system(center, mid_sources)
    pre = j.mid().inv()
    assert not pre.det().contains(0)
    _, cross = system(center, sources, cross=True)
    projected = pre*f + (pre*cross)*arb_mat([[arb(0, v.upper())] for v in widths])
    offset = max(abs(v).upper() for v in projected.entries())
    radius = 2*offset
    box = [v + arb(0, radius.upper()) for v in center]
    _, wide = system(box, sources)
    unit = arb_mat([[int(i==j) for j in range(12)] for i in range(12)])
    contraction = coupling["outer_sum"](unit-pre*wide)
    print("十二味方程的参数盒半宽与自映比上界", radius.str(9), (offset/radius+contraction).upper().str(9), flush=True)
    assert offset/radius + contraction < 1
    step, deriv_w = local(seed, box, DIRECTIONS, identity)
    step0, _ = local(seed, center, DIRECTIONS, identity, False)
    dst = adj(spectrum(shell, 6, NUMBER)[1])
    source = whole**3*incoming
    overlap = dst*step0**3*point**3*incoming

    def d_cube(unitary, tangent):
        return tangent*unitary**2 + unitary*tangent*unitary + unitary**2*tangent

    # 在两组真实 Cayley 坐标上分别用均值公式，保留矩阵乘法的抵消。
    fifth_changes = [dst*d_cube(step, e)*source for e in deriv_w]
    fourth_changes = [dst*step0**3*d_cube(whole, e)*incoming for e in deriv]
    for i in range(19):
        for j in range(2):
            err = sum(abs(e[i, j]).upper()*radius for e in fifth_changes)
            err += sum(abs(e[i, j]).upper()*r for e, r in zip(fourth_changes, widths, strict=True))
            delta = arb(0, err.upper())
            overlap[i, j] += acb(delta, delta)
    low, h, k = levels
    sites = [(low + (h-low)*abs(overlap[i, 0])**2 + (k-low)*abs(overlap[i, 1])**2)**2 for i in range(19)]
    for i, site in enumerate(sites, 1):
        print("层与截止平方", i, site.str(9), flush=True)
    assert [i+1 for i, s in enumerate(sites) if s > arb("1e-6")] == [1, 3, 4, 15]
    assert all(site < arb("1e-22") for i, site in enumerate(sites, 1) if i not in (1, 3, 4, 15))
    for i, cap in zip((15, 16, 17, 18), ("1.1e-33", "2.4e-34", "5.2e-34", "3.9e-30"), strict=True):
        assert sites[i] < arb(cap)
    weights = next(row[2][2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER) if row[0] == 6)
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [NUMBER])
    m = []
    for i, site in enumerate(sites, 1):
        if i in (1, 3, 4, 15):
            m.append(moments(field, heat, site))
        else:
            # 对非负质量有 E_2<=1 与 E_3<=1/2，避免在极低层对区间作除法。
            x = site.upper()
            k_cap = x*sum(abs(a) for _, a, _ in field)/(2*arb.pi())
            v_cap = x*x*sum(abs(g) for g, _, _ in field)/(8*arb.pi())
            m.append((arb(0, k_cap.upper()), arb(0, v_cap.upper())))
    thermal = [sum(c.real*v[j] for c, v in zip(weights, m, strict=True)) for j in range(2)]
    print("牛顿项与真空项", *[v.str(10) for v in thermal], flush=True)
    assert -25 < thermal[0] < -24 and arb("1.6") < thermal[1] < arb("1.7")
    return sites


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
