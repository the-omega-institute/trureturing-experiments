"""平方弦长不超过56时，每条边一次的同长共门构造与三壳六个引力端点。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


PARAMETERS = ("-2.875", ".30117250539529", "-.9375", ".875",
              "1.19824447190498", "1.0625", "-1.1875", "-.8125",
              ".75179797155578", "1.25", "-2.25", "2.38246246339496",
              "1.77884899248099", "-.9375", "2.0625", "-1.58897427908795")
ACTIVE = (1, 4, 8, 11, 12, 15)


def hessian_columns(station, layers, point, block):
    """共同门在整个角盒中的重列、梯度与对称二阶导数。"""
    from flint import acb, acb_mat
    angles = list(map(arb, PARAMETERS))
    for k, q in zip(ACTIVE, point, strict=True):
        angles[k] = q
    keys = [(i, j) for i in range(6) for j in range(i+1)]
    value = acb_mat(station)
    first = [acb_mat(19, 2) for _ in ACTIVE]
    second = [acb_mat(19, 2) for _ in keys]
    phase_signs = ((1, 0, 0, -1), (0, -1, 1, 0), (1, 1, 1, 1))
    for gate_id, edges in layers:
        g, dg = block(angles[4*gate_id:4*gate_id+4])
        local = {i: a%4 for i, a in enumerate(ACTIVE) if a//4 == gate_id}
        d2g = {}
        for i, a in local.items():
            for j, b in local.items():
                if i < j:
                    continue
                u, v = sorted((a, b))
                if u == v == 0:
                    d2g[i, j] = [-entry for entry in g]
                elif u == 0:
                    d2g[i, j] = [acb(0, k)*x for k, x in zip(phase_signs[v-1], dg[0], strict=True)]
                else:
                    d2g[i, j] = [-k*l*x for k, l, x in zip(phase_signs[u-1], phase_signs[v-1], g, strict=True)]

        def act(gate, col):
            u, v = col
            return (gate[0]*u + gate[1]*v, gate[2]*u + gate[3]*v)

        for a, b in edges:
            for col in range(2):
                old = value[a, col], value[b, col]
                grad = [(f[a, col], f[b, col]) for f in first]
                for m, (i, j) in zip(second, keys, strict=True):
                    terms = [act(g, (m[a, col], m[b, col]))]
                    if i in local:
                        terms.append(act(dg[local[i]], grad[j]))
                    if j in local:
                        terms.append(act(dg[local[j]], grad[i]))
                    if (i, j) in d2g:
                        terms.append(act(d2g[i, j], old))
                    m[a, col], m[b, col] = map(sum, zip(*terms, strict=True))
                for i, m in enumerate(first):
                    u, v = act(g, grad[i])
                    if i in local:
                        ra, rb = act(dg[local[i]], old)
                        u, v = u+ra, v+rb
                    m[a, col], m[b, col] = u, v
                value[a, col], value[b, col] = act(g, old)
    return value, first, second, keys


def hessian(station, layers, point, block, shells, levels):
    value, grad, second, keys = hessian_columns(station, layers, point, block)
    low, heavy, top = levels
    gains = (heavy-low, top-low)
    answer = []
    for left, weights, _, field in shells:
        p = left*value
        ds = [left*m for m in grad]
        d2 = [left*m for m in second]
        rows = [[arb(0) for _ in keys] for _ in range(2)]
        for i, w in enumerate(weights):
            a = low + sum(g*abs(p[i,k])**2 for k,g in enumerate(gains))
            left_end = max(low.lower(), a.lower())
            right_end = min(top.upper(), a.upper())
            assert 0 < left_end <= right_end
            domain = left_end.union(right_end)
            dx = [2*sum(g*(p[i,k].conjugate()*d[i,k]).real for k,g in enumerate(gains)) for d in ds]
            dk, dw, ddk, ddw = (arb(0) for _ in range(4))
            for g, alpha, mass in field:
                u_low, u_high = mass/right_end**2, mass/left_end**2
                exp_low, exp_high = (-u_low).exp(), (-u_high).exp()
                exp = exp_high.union(exp_low)
                # u exp(-u) 的唯一内部最大值为 1/e。
                scaled = (u_low*exp_low).union(u_high*exp_high)
                if not (u_low > 1 or u_high < 1):
                    scaled = scaled.union((-arb(1)).exp())
                dk += alpha*exp
                dw += g*exp
                ddk += alpha*(exp+2*scaled)
                ddw += g*(3*exp+2*scaled)
            slope = (domain*dk/arb.pi(), -domain**3*dw/(2*arb.pi()))
            curvature = (ddk/arb.pi(), -domain**2*ddw/(2*arb.pi()))
            for h, ((u,v), d) in enumerate(zip(keys, d2, strict=True)):
                mixed = 2*sum(g*(p[i,k].conjugate()*d[i,k]+ds[u][i,k].conjugate()*ds[v][i,k]).real
                              for k,g in enumerate(gains))
                for row, s, c in zip(rows, slope, curvature, strict=True):
                    row[h] += w*(s*mixed+c*dx[u]*dx[v])
        answer.extend(rows)
    return arb_mat(answer), keys


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    ordered = lib("two-shell-ordered-bond-gates")
    layers = ordered["LAYERS"] + ((3, ((11, 9), (12, 15), (18, 16))),
                                 (3, ((11, 13), (12, 10), (18, 17))))
    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    used = [((12, 32, 52, 56)[k], min(a, b), max(a, b)) for k, pairs in layers for a, b in pairs]
    lengths = {(sum((x-y)**2 for x, y in zip(sites[i], sites[j], strict=True)), i, j)
               for i, j in combinations(range(19), 2)}
    assert len(used) == len(set(used)) == 18 and set(used) == {e for e in lengths if e[0] <= 56}
    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    shells = []
    for n in range(2_306_888_100_092, 2_306_888_100_095):
        _, _, _, basis, _, weights, heat, field = geometry(root, fifth, columns, levels, number=n)
        shells.append((basis.conjugate().transpose(), weights, heat, field))

    def system(point, values):
        # 先乘两条重源列，直接在十九站上的两列传播门与六个导数。
        gates = ordered["transport"](fifth*columns, point, parameters=PARAMETERS, active=ACTIVE, layers=layers)
        fs, js = [], []
        for basis, weights, heat, field in shells:
            projection, *directions = [basis * g for g in gates]
            f, j = gravity(projection, directions, levels, weights, heat, field, values)
            fs.extend(f.tolist())
            js.extend(j.tolist())
        return arb_mat(fs), arb_mat(js)

    width = arb(".0002")
    center = [arb(PARAMETERS[j]) for j in ACTIVE]
    f, j = system(center, True)
    pre = j.mid().inv()
    assert not pre.det().contains(0)
    eye = arb_mat([[int(i == j) for j in range(6)] for i in range(6)])
    offset = max(z.abs_upper() for z in (pre*f).entries())
    print("角盒半径与预条件中心位移", width.str(8), offset.str(13), flush=True)
    assert offset < width
    full = [q+arb(0,width) for q in center]
    curv, keys = hessian(fifth*columns, layers, full, ordered["block"], shells, levels)
    projected = pre * curv
    center_defect = eye-pre*j
    ratios = []
    for i in range(6):
        derivative = sum(center_defect[i,k].abs_upper() for k in range(6))
        curvature = sum(projected[i,k].abs_upper()*(1 if a==b else 2)
                        for k,(a,b) in enumerate(keys))
        offset_i = (pre*f)[i,0].abs_upper()
        ratios.append(offset_i/width + derivative + curvature*width/2)
    print("二阶自映的逐行比例", *[q.str(12) for q in ratios], flush=True)
    assert max(ratios) < arb(".762")
    return [q + arb(0, width) for q in center]


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root, *source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
