"""以同一有符号固有时截止同时匹配引力和全部场的有限动量一圈交汇。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def angular_rule(order, last):
    rule = [arb.legendre_p_root(order, i, weight=True) for i in range(order)]
    nodes = []
    for j in range(2, last+1):
        left = arb(2)**-j
        for node, weight in rule:
            u = left*(3+node)/2
            v, a = u*(1-u), (1-2*u)**2
            omega = (6*v, 3*a, (8-3*a)/7, 3*(4-a)/11)
            nodes.append((v, [left*weight*f for f in omega]))
    ratio = arb(5)/14
    # 含端点小段、Cauchy 余项和 Gauss 正权带来的两个多项式误差。
    error = 6*arb(2)**-last+40*ratio**(2*order)/(1-ratio)
    return nodes, error


def angular_exposure(kind, z, x, p0, p1, logratio, rule):
    nodes, error = rule
    upper = logratio*(-z/x).exp()
    if upper < arb("1e-16"):
        return arb(0).union(upper)
    if z > 0:
        # 积分变量 p 上用 h/(z+ph) <= 1/(4z)。
        upper = upper.min((p1-p0)/(4*z)*(-z/x).exp())
        if upper < arb("1e-15"):
            return arb(0).union(upper)
    result = arb(0)
    for h, coefficients in nodes:
        result += coefficients[kind]*(((z+p0*h)/x).expint(1)-((z+p1*h)/x).expint(1))
    return result+arb(0, (error*logratio).upper())


def joint_system(root, energy_gev, *, scalar_shifts=None, order=14, last=41, number=2_306_888_100_085):
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, (rows,) = fields(root, (number,))
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 192
    fermion, scalar, vector, gauge = range(4)
    light = [3*sum(w[i] for w in model["weyl"]) for i in range(3)]
    species = [(fermion, arb(0), light), (scalar, arb(0), model["higgs"]),
               (gauge, arb(0), [0, -arb(22)/3, -11])]
    assert all((sum(b[i] for _, _, b in species)-model["beta"][i]).contains(0) for i in range(3))
    assert tuple(model["charges"]) == (1, 4, 2, 3, 6)
    # 合并场行保留原谱，另组实际背景的引力行。十六层标量可分别改变质量平方比。
    if scalar_shifts is None:
        scalar_shifts = [arb(0)]*16
    assert len(scalar_shifts) == 16
    gravity_rows = list(rows)
    for j, (q, weyl) in enumerate(zip(model["charges"], model["weyl"], strict=True)):
        for layer, (g, alpha, z) in enumerate(rows[1+16*j:17+16*j]):
            species.append((fermion, z, [2*v for v in weyl]))
            if q == 3:
                shift = arb(scalar_shifts[layer])
                assert z+shift > 0
                species.append((scalar, z+shift, model["higgs"]))
                if not shift.is_zero():
                    gravity_rows[1+16*j+layer] = (g-4, alpha-arb(4)/6, z)
                    gravity_rows.append((4, arb(4)/6, z+shift))
    species.extend((vector, z, [0, -7, -arb(21)/2]) for _, _, z in rows[-18:])
    assert len(species) == 117
    p0 = (model["mz"]/arb("1.22089e19"))**2
    p1 = (arb(energy_gev)/arb("1.22089e19"))**2
    logratio = (p1/p0).log()
    sites = list(map(arb, ("4e-15", "3.5e-11", "1.23e-6", ".5", "1")))
    fine = angular_rule(order, last)
    pi = arb.pi()
    cols = []
    for x in sites:
        k, w = single["moments"](gravity_rows, lib["expint"], x)
        shifts = [arb(0), arb(0), arb(0)]
        for kind, z, b in species:
            term = angular_exposure(kind, z, x, p0, p1, logratio, fine)/(4*pi)
            for i in range(3):
                shifts[i] -= b[i]*term
        cols.append([arb(1), k, w, *shifts])
    matrix = arb_mat([c[:3]+[c[3]-c[4], c[3]-c[5]] for c in cols]).transpose()
    target = arb_mat([[1], [1], [arb("1e-122")], [initial[1]-initial[0]], [initial[2]-initial[0]]])
    return dict(sites=sites, rows=rows, matrix=matrix, target=target, initial=initial, cols=cols,
                heat=lib["expint"], species=species, p0=p0, p1=p1, fine_error=fine[1]*logratio)


def verify(root, return_measure=False, *, scalar_shifts=None, check_path=None):
    data = joint_system(root, "1e16", scalar_shifts=scalar_shifts)
    sites, rows, matrix, target = (data[k] for k in ("sites", "rows", "matrix", "target"))
    cols, initial = data["cols"], data["initial"]
    det = matrix.det()
    assert -819 < det < -817
    weights = matrix.solve(target).entries()
    assert all(sign*w > 0 for sign, w in zip((1, 1, 1, -1, 1), weights, strict=True))
    assert sum(weights[1:]) > 0 and sum(weights[2:]) > 0 > sum(weights[3:])
    assert arb("1.4641238432") < -weights[3] < arb("1.4641238434")
    common = initial[0]+sum(c[3]*w for c, w in zip(cols, weights, strict=True))
    assert arb("39.12824") < common < arb("39.12826")
    if check_path is None:
        check_path = not return_measure
    if not check_path:
        if return_measure:
            return sites, weights, rows
        return
    print("精细角积分的绝对余量、行列式、五项权、统一逆耦合",
          data["fine_error"].str(12), det.str(12), *[w.str(13) for w in weights], common.str(13), flush=True)

    floors = running_floors(data, weights)
    assert min(floors) > 2
    print("从 M_Z 到交汇能标的三支逆耦合全程下界", *[v.str(12) for v in floors])
    if return_measure:
        return sites, weights, rows


def running_floors(data, weights, steps=28):
    species, sites, initial = (data[k] for k in ("species", "sites", "initial"))
    p0, p1 = data["p0"], data["p1"]
    pi, logratio = arb.pi(), (p1/p0).log()
    # 单个物种、单个截止的减法核随 p 单调。按各系数分离正负贡献，覆盖整条实能标路径。
    coarse = angular_rule(7, 17)
    previous_positive = [arb(0)]*3
    floors = list(initial)
    for j in range(1, steps+1):
        log = logratio*j/steps
        p = p1 if j == steps else p0*log.exp()
        positive, negative = [arb(0)]*3, [arb(0)]*3
        for weight, x in zip(weights, sites, strict=True):
            for kind, z, b in species:
                term = angular_exposure(kind, z, x, p0, p, log, coarse)/(4*pi)
                for i in range(3):
                    if b[i] == 0:
                        continue
                    coefficient = -weight*b[i]
                    if coefficient > 0:
                        positive[i] += coefficient*term
                    else:
                        assert coefficient < 0
                        negative[i] += coefficient*term
        floors = [min(f.lower(), (a+l+n).lower())
                  for f, a, l, n in zip(floors, initial, previous_positive, negative, strict=True)]
        previous_positive = positive
    return floors


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
