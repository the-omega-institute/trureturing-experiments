"""由两处有限精度的轻子率排除任意固定方向的秩一裸交换匹配。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

root = Path(__file__).parent
contacts = run_path(str(root/"ch39-contact-current-direction.py"))
rates = run_path(str(root/"ch39-polarized-lepton-rates.py"))
light, poles, mat, scheme, q = (contacts[k] for k in ("light", "poles", "mat", "scheme", "q"))
window = scheme["uncertain"]


def measurement_error(graph, tolerance):
    signs = rates["sign_certificate"](graph, tolerance)
    rows = rates["rows"](rates["basis_names"])
    nominal = (arb_mat(rows)*rates["column"](graph)).entries()
    root_error = tolerance/(1+(1-tolerance).sqrt())
    assert all(s*a > 0 for s, a in zip(signs, nominal, strict=True))
    errors = arb_mat([[arb(0, (abs(a)*root_error).upper())] for a in nominal])
    return rates["as_matrix"](arb_mat(rows.inv())*errors)


def dual_basis(e, z):
    c = poles["colour"](-e, poles["unit"])/24
    r = poles["right"](-e, poles["unit"])/4
    s = 9*c+r
    return arb_mat([[6*r/s, 18*c/s, 0], [3/s, -1/s, 2/(z+s)], [0, 0, 1]])


def contact_constraints(graph, measurement, polar, polar0, nominal, e, z, length, targets):
    # 已知极化的预解式差只含观测误差和极化误差。
    eye = scheme["cal"]["eye"]
    varied = graph+measurement
    delta = (eye+varied*polar).inv()*(measurement-varied*(polar-polar0)*graph)*(eye+polar0*graph).inv()
    _, g, v, alpha, beta = contacts["axes"](e, z)
    inverse = dual_basis(e, z)
    g0, v0, kg0, p0 = nominal
    projected_g, projected_v = inverse*g0, inverse*v0
    # 超荷在新、旧基中均为第一根轴，对第二列贡献严格为零。
    observed = (inverse*delta*inverse.transpose() + projected_g*projected_g.transpose()*kg0
                + projected_v*projected_v.transpose()*p0)
    low, high = [1/(poles["offset"]+target/(15*length)) for target in targets]
    margin_low = observed[1, 1]-low*beta
    margin_high = high*beta-observed[1, 1]
    margin = min(margin_low.lower(), margin_high.lower())
    assert margin > 0
    t_y = arb(0, (abs(observed[0, 1])/margin).upper())
    t_v = arb(0, (abs(observed[2, 1])/margin).upper())
    y = mat["y"]
    first = g[0, 0]+t_y*y[0, 0]+t_v*v[0, 0]
    assert abs(first) > 0
    slope = g[1, 0]/g[0, 0]
    # v 的前两个坐标平行于 g，因此分子只含超荷投影。
    direction = slope+t_y*(y[1, 0]-slope*y[0, 0])/first
    schur_error = observed[0, 1]**2/margin
    hyper_shift = observed[0, 0]+arb(0, abs(schur_error).upper())
    return direction, alpha, hyper_shift


def verify():
    time = poles["study"]["nominal"][2]
    assert poles["verify_ratio_monotonicity"]() > 0
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11432)]
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        (scale, z0, rho, raw_d, raw_b, masses), _, _ = poles["retune"]["predictions"](logs,time)
        c1, c2 = poles["reference_readings"](logs)[:2]
        mass_window = window(q("1e-11"))
        length, z, _ = poles["charged_geometry"](c1*mass_window, c2*mass_window)
        e = scheme["anchor"]/length
        triplet = masses[0]*(scale/length).sqrt()*window(q("1e-7"))
        data = light["graphs"](logs,time)
        log_x = scheme["logs_x"][0]
        _, polarization, x = scheme["cal"]["shape"]["probe"](log_x,e,z,raw_b/rho,raw_d/rho,triplet)
        e0 = scheme["anchor"]/scale
        _, polar0, _ = scheme["cal"]["shape"]["probe"](log_x,e0,z0,raw_b/rho,raw_d/rho,masses[0])
        directions, alpha_values, hyper_values = [], [], []
        for point, multiplier, polar in ((0, q(1), arb_mat(3,3)), (1, x, polarization)):
            energy = multiplier*e0
            h0, g0, v0, alpha0, beta0 = contacts["axes"](energy, z0)
            nominal = (g0, v0, beta0/rho, energy/(raw_b+raw_d*energy+4*rho/h0))
            measurement = measurement_error(data[point], q("1e-11"))
            direction, alpha, shift = contact_constraints(
                data[point], measurement, polar, polar0 if point else polar,
                nominal, multiplier*e, z, length, targets)
            directions.append(direction)
            alpha_values.append(alpha)
            hyper_values.append(shift+alpha0/rho)
        first, second = directions
        print("case", case, "slopes", first.lower().str(9), first.upper().str(9),
              second.lower().str(9), second.upper().str(9), flush=True)
        assert q("-2.1790") < first < q("-2.1720")
        assert q("-2.1980") < second < q("-2.1920")
        assert first-second > q(".0130")
        delta_alpha = alpha_values[1]-alpha_values[0]
        assert delta_alpha > 0
        k = (hyper_values[1]-hyper_values[0])/delta_alpha
        arithmetic = 15*length*(1/k-poles["offset"])
        print("common hypercharge, target endpoints", arithmetic.lower().str(12),
              arithmetic.upper().str(12), flush=True)
        assert targets[1] < arithmetic < targets[0]


if __name__ == "__main__":
    verify()
