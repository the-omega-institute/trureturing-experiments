"""第五壳固定上游的共同味方程，可叠加第三枢轴的四项端点。"""
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


NUMBER = 2_306_888_100_090


def adj(matrix):
    return matrix.conjugate().transpose()


def solve(root, name, quaternion, directions, *, endpoints=None):
    size = len(directions)
    assert size == (16 if endpoints else 12)
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    point, whole, deriv, widths, upstream, incoming, levels = lib("fourth-third-pivot-joint")["verify"](
        root, return_chart=True)
    quarks = lib("fifth-two-charge-full-quark")
    rotation = quarks["rotation"]
    turn = rotation(quaternion) * rotation(quarks["QUATERNION"])
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
        root, eye, name+"-cayley.tsv")

    dst6 = adj(spectrum(shell, 6, NUMBER)[1])
    low, heavy, highest = levels
    def d_cube(unitary, tangent):
        return tangent*unitary**2 + unitary*tangent*unitary + unitary**2*tangent
    mid6 = point**3*incoming
    source6 = whole**3*incoming
    if endpoints:
        weights = next(row[2][2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER)
                       if row[0] == 6)
        moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [NUMBER])
        ctx.prec = 384
        exposure = endpoints(dst6*seed**3*mid6, (low, heavy, highest))
        ctx.prec = 384
        tangent6 = [d_cube(whole, e)*incoming for e in deriv]

    def system(coords, srcs, src6, *, cross=False):
        z, dz = local(seed, coords, directions, identity, not cross)
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
        if endpoints:
            ov = dst6*z**3*src6
            delta = ([dst6*z**3*e for e in tangent6] if cross else
                     [dst6*d_cube(z, e)*src6 for e in dz])
            amps = [low + (heavy-low)*abs(ov[i, 0])**2 + (highest-low)*abs(ov[i, 1])**2 for i in range(19)]
            totals = [-arb(1), -arb("1e-122")]
            gauge_values, gauge_slopes = exposure(amps)
            slopes = [[], []]
            for i, a in enumerate(amps):
                assert a > 0
                c = weights[i].real
                for j, v in enumerate(moments(field, heat, a*a)):
                    totals[j] += c*v
                ex = [(-m/(a*a)).exp() for _, _, m in field]
                slopes[0].append(c*a*sum(b*e for (_, b, _), e in zip(field, ex, strict=True))/arb.pi())
                slopes[1].append(-c*a**3*sum(g*e for (g, _, _), e in zip(field, ex, strict=True))/(2*arb.pi()))
            d_amp = arb_mat([[2*((heavy-low)*(ov[i, 0].conjugate()*e[i, 0]).real +
                                 (highest-low)*(ov[i, 1].conjugate()*e[i, 1]).real) for e in delta] for i in range(19)])
            values.extend([*totals, *gauge_values])
            jac.extend((arb_mat(slopes + gauge_slopes)*d_amp).tolist())
        return arb_mat([[v] for v in values]), arb_mat(jac)

    center = [arb(0)]*size
    for _ in range(3):
        f, j = system(center, mid_sources, mid6)
        step = j.mid().solve(f.mid()).entries()
        center = [(x-d).mid() for x, d in zip(center, step, strict=True)]
    f, j = system(center, mid_sources, mid6)
    pre = j.mid().inv()
    assert not pre.det().contains(0)
    _, cross = system(center, sources, source6, cross=True)
    projected = pre*f + (pre*cross)*arb_mat([[arb(0, v.upper())] for v in widths])
    offset = max(abs(v).upper() for v in projected.entries())
    radius = 2*offset
    box = [v + arb(0, radius.upper()) for v in center]
    _, wide = system(box, sources, source6)
    unit = arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
    contraction = coupling["outer_sum"](unit-pre*wide)
    print("联合方程维数、参数盒半宽与自映比上界", size, radius.str(9), (offset/radius+contraction).upper().str(9), flush=True)
    assert offset/radius + contraction < 1
    step, deriv_w = local(seed, box, directions, identity)
    step0, _ = local(seed, center, directions, identity, False)
    overlap = dst6*step0**3*mid6

    # 在两组真实 Cayley 坐标上分别用均值公式，保留矩阵乘法的抵消。
    fifth_changes = [dst6*d_cube(step, e)*source6 for e in deriv_w]
    fourth_changes = [dst6*step0**3*d_cube(whole, e)*incoming for e in deriv]
    for i in range(19):
        for j in range(2):
            err = sum(abs(e[i, j]).upper()*radius for e in fifth_changes)
            err += sum(abs(e[i, j]).upper()*r for e, r in zip(fourth_changes, widths, strict=True))
            delta = arb(0, err.upper())
            overlap[i, j] += acb(delta, delta)
    sites = [(low + (heavy-low)*abs(overlap[i, 0])**2 + (highest-low)*abs(overlap[i, 1])**2)**2 for i in range(19)]
    return sites
