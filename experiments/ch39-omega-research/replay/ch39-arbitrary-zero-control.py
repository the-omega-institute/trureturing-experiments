"""全谱任意切换下的 CKM 与规范、端口联合检验。"""
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

root = Path(__file__).parent
local = run_path(str(root / "ch39-center-free-neighborhood-time.py"))
free = run_path(str(root / "ch39-center-free-three-pulses.py"))
study, box, ordered, flavor = (free[key] for key in ("study", "box", "ordered", "flavor"))
entries = free["entries"]
center = 11431


def nuclear(matrix):
    return (sum(abs(z)**2 for z in matrix.entries()) + 2*abs(matrix.det())).sqrt()


def operator_norm(matrix):
    mean = matrix.trace().real/2
    return abs(mean) + ((matrix[0, 0].real-mean)**2 + abs(matrix[0, 1])**2).sqrt()


def background():
    point, scale, eps = local["model"]()
    coefficients = [v + arb(0, arb(2)**-28) for v in box["centers"]]
    response = study["native"]["plane_response"]
    basis = flavor["basis"]
    center_height = point(center)[0]
    center_response = response(center_height)
    modes = free["three"]["spectral_modes"](center_response, coefficients)
    reference = flavor["radial"] + basis*ordered["evolution"](center_response, arb(1), coefficients)*basis.transpose()
    values = [reference[i, j] for i, j in entries]
    projectors = [acb_mat([[basis[i, a]*basis[j, b]/value for b in range(2)] for a in range(2)])
                  for (i, j), value in zip(entries, values)]
    return point, scale, eps, coefficients[0], center_response, modes, values, projectors


def base_tests(projectors):
    # 每项记录实线性泛函、规范乘子、端口乘子。泛函值是逐项乘积和的实部。
    tests = []
    def add(name, matrix, gauge, price):
        tests.append((name, matrix, arb(gauge), arb(price)))
    add("norm", None, ".000865", ".007314")
    for i, name, gauges, prices in (
        (0, "us", (".000942", ".000010"), (".00055", ".01417")),
        (1, "cb", (".000572", "-.000300"), (".000953", ".001410")),
        (3, "cs", (".000684", "-.000748"), (".003105", ".000411")),
    ):
        for sign, gauge, price in zip((1, -1), gauges, prices):
            add(name + str(sign), -acb(0, sign)*projectors[i], gauge, price)
    return tests


def integrated_bounds(data, tests, norm_bonus, center_price, cell_count, center_limit, budget_factor=1):
    point, scale, eps, eta, center_response, modes, _, _ = data
    cost, *_ = local["support"]["prev"]["setup"]()
    cap = 2*eps/scale
    response = study["native"]["plane_response"]
    dual_norms = [arb(1) if m is None else nuclear(m) for _, m, _, _ in tests]
    max_contrast = abs(eta)*2/(arb("1.3")-4*arb.pi()/32)
    floor = [arb(-1) if name != "norm" and name not in norm_bonus else
             arb(-150) if name in norm_bonus else arb("-0.1").upper() for name, *_ in tests]

    low, high = 9400, 14000
    for n, direction in ((low, 1), (high, -1)):
        _, x, c = point(n)
        _, derivative = cost(x)
        for (name, _, g, h), nn, cut in zip(tests, dual_norms, floor):
            assert direction*(10000*g + h*derivative/cap) > 0
            assert (nn+norm_bonus.get(name,0))*max_contrast - g*10000*x - h*c/cap < cut, n
    freq = modes[1][0] - modes[0][0]
    phases = [(acb(0, freq*(arb(k).union(k+1)/cell_count)).exp()) for k in range(cell_count)]
    envelopes = [[cut for _ in phases] for cut in floor]
    def linear(matrix, wave):
        return sum(a*b for a,b in zip(matrix.entries(), wave.entries()))
    def examine(x, c, block):
        diff = eta*(block-center_response)
        norm = operator_norm(diff)
        waves = []
        for a, qa in modes:
            for _, qb in modes:
                waves.append(acb(0, 1)*acb(0, a).exp()*qa*diff*qb)
        for j, ((name, matrix, g, h), nn) in enumerate(zip(tests, dual_norms)):
            penalty = 10000*g*x + h*c/cap - norm_bonus.get(name,0)*norm
            if nn*norm-penalty < floor[j]:
                continue
            if matrix is None:
                top = (norm - penalty).upper()
                if top > envelopes[j][0]:
                    envelopes[j] = [top for _ in phases]
                continue
            w0, w1, w2, w3 = [linear(matrix, wave) for wave in waves]
            constant = (w0 + w3).real - penalty
            harmonic = w1 + w2.conjugate()
            if constant + abs(harmonic) < floor[j]:
                continue
            # 包围整个时间格；背景本征值之差与当前编号无关。
            for k, phase in enumerate(phases):
                top = (constant + (harmonic*phase).real).upper()
                if top > envelopes[j][k]:
                    envelopes[j][k] = top
    # 中心附近先求一个全时段有效的下截值，远处被截掉的上界已包含于包络。
    for n in (11432, 11430, 11429):
        height, x, c = point(n)
        examine(x, c, response(height))
    floor = [min(row) for row in envelopes]
    for first in range(low+1, high, 128):
        for n, zero in enumerate(acb.zeta_zeros(first, min(128, high-first)), first):
            if n in (11429, 11430, 11431, 11432):
                continue
            x = study["center"]/(zero.imag**2+study["offset"])-1
            c, _ = cost(x)
            if all((nn+norm_bonus.get(name,0))*max_contrast - 10000*g*x - h*c/cap < cut
                   for (name, _, g, h), nn, cut in zip(tests, dual_norms, floor)):
                continue
            examine(x, c, response(zero.imag))

    bounds = {}
    for (name, _, g, h), row in zip(tests, envelopes):
        price = center_price.get(name, max(arb(0), -min(row)))
        row = [max(v, -price) for v in row]
        correction = center_limit*price
        bounds[name] = sum(row)/cell_count + budget_factor*(abs(g)*10000*eps + h) + correction
    return bounds


def verify(second_angle=False):
    data = background()
    values, projectors = data[-2:]
    beta = arb(932)/1000 if second_angle else arb(34)/100
    angle = 1 if second_angle else 0
    max_tilt = arb(".002") if second_angle else arb(".011")
    unit_angle, unit_phase = arb(".000055"), arb(".006")
    signs = (1, 1, -1, -1)
    tests = base_tests(projectors)
    phase_matrix = sum((-acb(0, s)*p for s, p in zip(signs, projectors)), acb_mat(2, 2))
    norm_bonus = {}
    center_price = {}
    if second_angle:
        cases = [(f"small{sign}", sign, "-2.418", "1.015", "1.715", "0") for sign in (-1, 1)]
        cases += [(f"large{sign}", sign, "-1.057", "1.161", "0", "1089") for sign in (-1, 1)]
    else:
        cases = [("score-1", -1, "-2.521", ".793", None, "0"),
                 ("score1", 1, "-2.344", ".691", None, "0")]
    for name, sign, gauge, price, occupation_price, bonus in cases:
        tilt = sign*max_tilt
        matrix = (beta*abs(values[angle])/unit_angle * (1-acb(0, tilt))*projectors[angle]
                  + (1-beta)/unit_phase*phase_matrix)
        tests.append((name, matrix, arb(gauge), arb(price)))
        norm_bonus[name] = arb(bonus)
        if occupation_price is not None:
            center_price[name] = arb(occupation_price)
    cell_count = 4096 if second_angle else 1024
    share = 185 if second_angle else 46
    bounds = integrated_bounds(data, tests, norm_bonus, center_price, cell_count, arb(share)/1000)
    dual_norms = [arb(1) if m is None else nuclear(m) for _, m, _, _ in tests]
    norm = bounds["norm"]
    assert 0 < norm < arb(".012")
    remainder = norm**2/2
    attenuation = [((1-arb((1,5,29)[i])**2/867)*(1-arb((1,5,29)[j])**2/867)).sqrt()
                   for i, j in entries]
    # 所有相对坐标 z=V/V0-1 落在 Re(1+z)>0 内，确保连续的主辐角。
    crude = [norm*a/abs(v) for a,v in zip(attenuation, values)]
    assert all(v < 1 for v in crude)
    ybounds = {i: bounds[name+"1"].max(bounds[name+"-1"]) + remainder*attenuation[i]/abs(values[i])
               for i,name in ((0,"us"), (1,"cb"), (3,"cs"))}
    assert ybounds[angle]/(2*(1-crude[angle])) < max_tilt
    ub_min, ub_max = arb(".003677"), arb(".003787")
    c13min, c13max = (1-ub_max**2).sqrt(), (1-ub_min**2).sqrt()
    min_us, max_us = arb(".224955")*c13min, arb(".225065")*c13max
    min_cb, max_cb = arb(".041775")*c13min, arb(".041885")*c13max
    ratios = [(min_us/abs(values[0]),max_us/abs(values[0])),
              (min_cb/abs(values[1]),max_cb/abs(values[1])),
              (ub_min/abs(values[2]),ub_max/abs(values[2])),
              (1-crude[3],1+crude[3])]
    theta = {i:(y/ratios[i][0]).asin() for i,y in ybounds.items()}
    quartet = values[0]*values[1]*values[2].conjugate()*values[3].conjugate()
    phase0 = quartet.arg()
    # q=quartet+正实数；两者相角差由该正实数与 |quartet| 之比包围。
    correction = (max_us*max_cb*ub_max/((1-ub_max**2)*abs(values[3])*ratios[3][0])).asin()
    lo, hi = arb("1.141"), arb("1.153")
    gap = abs(phase0-lo).max(abs(phase0-hi)) + correction
    theta[2] = theta[0]+theta[1]+theta[3]+gap
    assert theta[2] < arb(".1") and sum(crude[i].asin() for i in range(4)) < 1
    # |theta-Im z| <= |1-rho| |theta| + rho |theta|^3/6 .
    angle_error = sum(abs(1-a).max(abs(1-b))*theta[i] + b*theta[i]**3/6
                      for i,(a,b) in enumerate(ratios))
    angle_floor = min_cb if second_angle else min_us
    offset = beta*(abs(values[angle])-angle_floor)/unit_angle + (1-beta)*(phase0-lo)/unit_phase
    assert 0 < lo < hi < arb.pi()/2 and 0 < phase0 < arb.pi()/2
    results = []
    for j in range(len(tests)-len(cases), len(tests)):
        name = tests[j][0]
        norm_limit = arb(".003") if name.startswith("small") else norm
        bonus = norm_bonus[name]*arb(".003")
        upper = (offset + bounds[name] - bonus + dual_norms[j]*norm_limit**2/2
                 + (1-beta)*angle_error/unit_phase)
        results.append(upper)
    scores = results
    print("Testing center share", share, "/ 1000")
    print("Integrated norm upper", norm.upper().str(12), "phase error upper", angle_error.upper().str(12),
          "score upper bounds", *(v.upper().str(12) for v in scores))
    assert all(v < 0 for v in scores)
    print("Certified center share greater than", share, "/ 1000")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--second-angle", action="store_true", help="改用第二混合角与分段范数界")
    verify(second_angle=parser.parse_args().second_angle)
