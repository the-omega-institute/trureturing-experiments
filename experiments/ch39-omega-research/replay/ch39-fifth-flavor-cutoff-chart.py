"""第五壳固定上游的共同味方程，可叠加第三枢轴的四项端点。"""
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


NUMBER = 2_306_888_100_090


def adj(matrix):
    return matrix.conjugate().transpose()


def solve(root, name, quaternion, directions, *, endpoints=None, frame=None, rows=(), target_radius=0, target_sites=None, seed_rays=False, high_targets=()):
    size = len(directions)
    assert size == (16 if endpoints else 12) + 4*len(rows) + len(high_targets)
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    point, whole, deriv, widths, upstream, incoming, levels = lib("fourth-third-pivot-joint")["verify"](
        root, return_chart=True, tight_chart=bool(high_targets))
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
    blocks = [turn*make(sector[-1].entries())*adj(turn) for sector in upstream]
    if frame is not None:
        blocks = frame(root, lib, blocks)
    for i, (charge, (src, _, w, _, _), projectors_q) in enumerate(zip((-4, 2), upstream, projectors, strict=True)):
        target = blocks[i]
        sector = (adj(spectrum(shell, charge, NUMBER)[1]), w, arb_mat([shape(r) for r in projectors_q]).transpose(),
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
    if seed_rays:
        reference = dst6*seed**3*point**3*incoming
        input_goals = []
        for i,x in zip(rows,target_sites,strict=True):
            first,second = [acb(reference[i,j].real.mid(),reference[i,j].imag.mid()) for j in range(2)]
            scale = ((x.sqrt()-low)/((heavy-low)*abs(first)**2+(highest-low)*abs(second)**2)).sqrt()
            input_goals.append((first*scale,second*scale))
    else:
        amplitudes = ([((x.sqrt()-low)/(heavy-low)).sqrt() for x in target_sites]
                      if target_sites is not None else [arb(0)]*len(rows))
        input_goals = [(acb(v),acb(0)) for v in amplitudes]
    row_goals=[getattr(v,part) for part in ("real","imag") for pair in input_goals for v in pair]
    high_goals=[x.sqrt() for i,x in high_targets]
    goal_centers=[c.mid() for c in row_goals+high_goals]
    assert len(input_goals)==len(rows)
    def d_cube(unitary, tangent):
        return tangent*unitary**2 + unitary*tangent*unitary + unitary**2*tangent
    mid6 = point**3*incoming
    source6 = whole**3*incoming
    if endpoints or rows:
        tangent6 = [d_cube(whole, e)*incoming for e in deriv]
    if endpoints:
        weights = next(row[2][2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER)
                       if row[0] == 6)
        moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [NUMBER])
        ctx.prec = 384
        exposure = endpoints(dst6*seed**3*mid6, (low, heavy, highest))
        ctx.prec = 384

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
            columns = []
            for do in dv:
                dw = arb_mat([[2*sum(w*(ov[k, j].conjugate()*do[k, j]).real for j, w in enumerate(ws))]
                              for k in range(19)])
                columns.append((cols*dw).entries())
            jac.extend(arb_mat(columns).transpose().tolist())
        if endpoints or rows:
            ov = dst6*z**3*src6
            delta = ([dst6*z**3*e for e in tangent6] if cross else
                     [dst6*d_cube(z, e)*src6 for e in dz])
        if rows:
            keys=[(i,j,part) for part in ("real","imag") for i in rows for j in range(2)]
            for n,(i,j,part) in enumerate(keys):
                values.append(getattr(ov[i,j],part)-goal_centers[n])
                jac.append([getattr(e[i,j],part) for e in delta])
            for n,(i,_) in enumerate(high_targets,len(keys)):
                a=low+(heavy-low)*abs(ov[i,0])**2+(highest-low)*abs(ov[i,1])**2
                values.append(a-goal_centers[n])
                jac.append([2*((heavy-low)*(ov[i,0].conjugate()*e[i,0]).real+
                               (highest-low)*(ov[i,1].conjugate()*e[i,1]).real) for e in delta])
        if endpoints:
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
    for _ in range(14 if high_targets else 3):
        f, j = system(center, mid_sources, mid6)
        step = j.mid().solve(f.mid()).entries()
        center = [(x-d).mid() for x, d in zip(center, step, strict=True)]
        if high_targets and max(abs(v).upper() for v in step)<arb("1e-18"):
            break
    f, j = system(center, mid_sources, mid6)
    pre = j.mid().inv()
    assert not pre.det().contains(0)
    _, cross = system(center, sources, source6, cross=True)
    projected = pre*f + (pre*cross)*arb_mat([[arb(0, v.upper())] for v in widths])
    unit = arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
    if rows:
        if target_sites is not None:
            error=[arb(0)]*12+[v-c for v,c in zip(row_goals+high_goals,goal_centers,strict=True)]
            projected -= pre*arb_mat([[v] for v in error])
        if target_radius:
            arbitrary = [arb(0)]*12+[arb(0,arb(target_radius).upper())]*(4*len(rows))
            projected -= pre*arb_mat([[v] for v in arbitrary])
        assert not endpoints
        errs = [abs(v).upper() for v in projected.entries()]
        _, j_center = system(center, sources, source6)
        defect = (unit-pre*j_center).tolist()
        curvature = []
        for line in pre.tolist():
            end=12+4*len(rows)
            total = 18*sum(abs(v).upper() for v in line[12:end])
            total += 72*highest*sum(abs(v).upper() for v in line[end:])
            for iq, (_, ws, cols, _) in enumerate(sectors):
                diagonal = arb_mat([line[6*iq:6*iq+6]])*cols
                total += 8*(2-iq)**2*sum(map(abs,ws))*max(abs(v).upper() for v in diagonal.entries())
            curvature.append(total.upper())
        maximum = max(errs)
        radii = [2*max(v,maximum/10000) for v in errs]
        for _ in range(10):
            # 令 s 为反厄米坐标增量的 Frobenius 范数。酉 Cayley 预解式的算子范数至多一。
            # 第 p 次幂的二阶余项至多 2p²s²，酉共轭的迹余项至多 8p²s² 乘观测与源的范数。
            square = sum((1 if i==j else 2)*r*r for (i,j,_),r in zip(directions,radii,strict=True))
            need = [e+k*square+sum(abs(a).upper()*r for a,r in zip(row,radii,strict=True))
                    for e,k,row in zip(errs,curvature,defect,strict=True)]
            ratio = max((n/r).upper() for n,r in zip(need,radii,strict=True))
            print("逐坐标联合盒半宽与二阶自映比", max(radii).str(9),ratio.str(9),flush=True)
            if ratio < 1:
                break
            radii=[max(r,arb("1.03")*n) for r,n in zip(radii,need,strict=True)]
        assert ratio < 1
        box=[v+arb(0,r.upper()) for v,r in zip(center,radii,strict=True)]
    else:
        offset = max(abs(v).upper() for v in projected.entries())
        radius = 2*offset
        radii = [radius]*size
        box = [v+arb(0,radius.upper()) for v in center]
        _,wide = system(box,sources,source6)
        contraction = coupling["outer_sum"](unit-pre*wide)
        print("联合方程维数、参数盒半宽与自映比上界", size, radius.str(9),
              (offset/radius+contraction).upper().str(9), flush=True)
        assert offset/radius+contraction < 1
    step, deriv_w = local(seed, box, directions, identity)
    step0, _ = local(seed, center, directions, identity, False)
    overlap = dst6*step0**3*mid6

    # 在两组真实 Cayley 坐标上分别用均值公式，保留矩阵乘法的抵消。
    fifth_changes = [dst6*d_cube(step, e)*source6 for e in deriv_w]
    fourth_changes = [dst6*step0**3*d_cube(whole, e)*incoming for e in deriv]
    for i in range(19):
        for j in range(2):
            err = sum(abs(e[i,j]).upper()*r for e,r in zip(fifth_changes,radii,strict=True))
            err += sum(abs(e[i, j]).upper()*r for e, r in zip(fourth_changes, widths, strict=True))
            delta = arb(0, err.upper())
            overlap[i, j] += acb(delta, delta)
    sites = [(low + (heavy-low)*abs(overlap[i, 0])**2 + (highest-low)*abs(overlap[i, 1])**2)**2 for i in range(19)]
    return sites
