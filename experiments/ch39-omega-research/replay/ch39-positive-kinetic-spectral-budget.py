"""用正定局部动能约束两点族流所允许的全部有隙混合谱。"""
from dataclasses import dataclass
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

four = run_path(str(Path(__file__).with_name("ch39-four-current-probes.py")))
poles, mat, q, form = (four[key] for key in ("poles", "mat", "q", "form"))
family = arb_mat([[0], [0], [1]])
momenta = [q(".5").exp(), q(1).exp()]
tolerance = q("2e-8")


@dataclass
class Case:
    scale: arb
    family_shape: arb
    mass_ratio_squared: arb
    geometry: list
    hypercharge: arb
    reference_coupling: arb


def reference(logs, time):
    (scale, _, rho, eta, _, masses), _, _ = poles["retune"]["predictions"](logs, time)
    neutral = poles["reference_readings"](logs)[3]
    geometry, hyper, k0, _, _ = four["setup"](logs, time, [])
    return Case(scale, eta/rho, masses[0]**2*scale/neutral**2, geometry, hyper, k0)


def spectral_forms(case, target, error):
    k = 1/(poles["offset"] + target/(15*case.scale))
    metric = arb_mat([[k/24, 0, 0], [0, k/4, 0], [0, 0, k/case.family_shape]])
    metric += mat["hypercharge"]*(case.hypercharge-k/24)
    assert mat["positive"](metric)
    matrices = []
    for sign, row in zip((-1, 1), case.geometry[1:], strict=True):
        tree, polarization, _, _, _, reference_exchange, reference_row = row
        inverse = (tree*k + mat["hypercharge"]*case.hypercharge).inv()-polarization
        assert mat["positive"](inverse)
        graph = inverse.inv()
        desired = k*reference_exchange*(1+sign*error)
        column = tree*family
        dressed = (polarization.inv()-tree*k-mat["hypercharge"]*case.hypercharge).inv()
        difference = (reference_row*(column+tree*dressed*column*k))[0, 0]
        gap = k*((case.reference_coupling-k)*difference+sign*error*reference_exchange)
        assert gap.overlaps(desired-graph[2, 2])
        assert gap > 0
        column = graph*family
        matrices.append((graph+column*column.transpose()*(1/gap), inverse, desired))
    return metric, matrices


def at_mass(case, matrices, mass):
    factors = [(x-1)/(case.mass_ratio_squared*mass**2+x) for x in momenta]
    return [f*row[0] for f, row in zip(factors, matrices, strict=True)], factors


def lower_bound(metric, forms, budget, multipliers):
    # A positive dual matrix excludes every spectrum with cost <= budget.
    a, b = forms
    return any(mat["positive"](metric-(a+(a-b)*z)*budget) for z in multipliers)


def witness(metric, matrices, forms, factors, site_ratio):
    # Fix one rational site ratio and solve the other component exactly.
    delta = forms[0]-forms[1]
    site = arb_mat([[1], [site_ratio], [0]])
    a, b, c = delta[2, 2], form(family, delta, site), form(site, delta, site)
    disc = b*b-a*c
    assert disc > 0 and not a.contains(0)
    costs = []
    for sign in (-1, 1):
        direction = site+family*((-b+sign*disc.sqrt())/a)
        v = direction*(1/form(direction, forms[0], direction).sqrt())
        costs.append((form(v, metric, v), v))
    cost, v = min(costs, key=lambda pair: float(pair[0].mid()))
    for factor, (_, inverse, target) in zip(factors, matrices, strict=True):
        updated = inverse-v*v.transpose()*factor
        assert mat["positive"](updated)
        residual = updated.inv()[2, 2]/target-1
        assert residual.contains(0) and abs(residual) < q("1e-8")
    return cost


# Decimal witnesses, in the order minimum cost, dual multiplier, site ratio,
# then limiting mass/M3, its dual multiplier and site ratio.
adjacent = (
    ((".14693821", "9.208901", ".3984534", "4.364220", "7.156620", ".4080804"),
     (".11858080", "10.268555", ".3935635", "4.711691", "7.717737", ".4054333")),
    ((".14728960", "9.207293", ".3985078", "4.359872", "7.157374", ".4081172"),
     (".11886780", "10.266614", ".3936226", "4.706999", "7.718518", ".4054718")),
    ((".14882508", "9.200327", ".3987434", "4.341042", "7.160659", ".4082770"),
     (".12012200", "10.258203", ".3938786", "4.686682", "7.721921", ".4056390")),
)


def check_adjacent(case, target, error, values):
    cost, dual, site, mass, edge_dual, edge_site = map(q, values)
    metric, matrices = spectral_forms(case, target, error)
    forms, factors = at_mass(case, matrices, q(2))
    # If the second constraint were slack, a local minimum on the
    # first ellipsoid would be an unconstrained Rayleigh minimum.
    assert metric[2, 2]/forms[0][2, 2] < q(".01")
    for factor, (_, inverse, _) in zip(factors, matrices, strict=True):
        assert mat["positive"](metric-inverse.inv()*factor)
    radius = q("1e-8")
    assert lower_bound(metric, forms, cost-radius, [dual])
    assert witness(metric, matrices, forms, factors, site) < cost+radius
    step = q(".000001")
    forms, factors = at_mass(case, matrices, mass-step)
    assert witness(metric, matrices, forms, factors, edge_site) < 1
    forms, _ = at_mass(case, matrices, mass+step)
    assert lower_bound(metric, forms, q(1), [edge_dual])


def zero_targets(time):
    model = poles["study"]["comp"]["model"]
    first, last = 10946, 11430
    heights = [z.imag for z in acb.zeta_zeros(first, last-first+1)]
    offset = q(model["fixed"]["OFFSET"])
    numerator = model["target"](first, time)*(heights[0]**2+offset)
    targets = {n: numerator/(h**2+offset) for n, h in enumerate(heights, first)}
    assert targets[last].overlaps(model["target"](last, time))
    return targets


def verify():
    time = poles["study"]["nominal"][2]
    targets = zero_targets(time)
    # Search only for dual witnesses. Every acceptance is an Arb sign test.
    multipliers = [q(i)/20 for i in range(180, 221)]
    multipliers += [q(i)/4 for i in range(0, 401)]
    for case_index, logs in enumerate(poles["retune"]["log_sets"]):
        case = reference(logs, time)
        for label, target in targets.items():
            if label < 11430:
                cutoff = 2 if label < 11423 else q("4.33")
                metric, matrices = spectral_forms(case, target, tolerance)
                forms, _ = at_mass(case, matrices, cutoff)
                if label >= 11423:
                    at_two, _ = at_mass(case, matrices, q(2))
                    assert lower_bound(metric, at_two, q(".17"), multipliers)
                if not lower_bound(metric, forms, q(1), multipliers):
                    raise ArithmeticError(f"no exclusion for {case_index=} {label=} {cutoff=}")
            if 11423 <= label <= 11430:
                metric, matrices = spectral_forms(case, target, q(0))
                forms, factors = at_mass(case, matrices, q(2))
                assert witness(metric, matrices, forms, factors, q(".4")) < 1
        for error, values in zip((q(0), tolerance), adjacent[case_index], strict=True):
            check_adjacent(case, targets[11430], error, values)
            print("case", case_index, "relative error", error.str(4), "UV fraction",
                  values[0], "+/- 1e-8; gap/M3", values[3], "+/- 1e-6", flush=True)


if __name__ == "__main__":
    verify()
