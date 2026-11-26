"""Interval bounds for the portal derivative of the finite-vertex light eigenvalue."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
base, signed, mass = (finite[key] for key in ("base", "signed", "mass"))
q, loop, unit = (finite[key] for key in ("q", "loop", "unit"))
portal = mass["portal"]


def tensor_radius(groups, quartic):
    square = 288 * sum(n * abs(quartic["l", i]).upper()**2 for i, n in groups.items())
    for (kind, indices), entry in quartic.items():
        if kind == "a":
            i, j = indices
            pairs = groups[i]*groups[j] if i != j else groups[i]*(groups[i]-1)//2
            square += 96*pairs*(abs(entry).upper()+2*abs(quartic["b", indices]).upper())**2
    return square.sqrt().upper()


def coefficients(groups, quartic, gauge, yukawa, gradient=False):
    # First the class sums, then an arbitrary individual from each class.
    n = len(groups)
    totals = signed["flow_matrix"](groups, quartic, gauge, yukawa, True if gradient else None)
    result = arb_mat(2*n, 2*n)
    for i, count in groups.items():
        same = portal(quartic, i, i) if count > 1 else q(0)
        for j in groups:
            result[i, j] = totals[i, j]
            result[n+i, j] = 2*(same if i == j else portal(quartic, i, j))
        result[n+i, n+i] = totals[i, i] - 2*count*same
    return result


def add_at_threshold(k, groups, quartic, scale, direction, state):
    g = 1 if k < 3 else 2
    for key in quartic:
        if key not in direction:
            direction[key] = q(int(key == ("a", (0, 2))))
    new_mass = scale**2 / unit
    for total, individual, couplings, injection in (
        (*state[:2], quartic, new_mass), (*state[2:], direction, q(0))
    ):
        if len(total) < len(groups):
            total.append(q(0))
            individual.append(injection)
        for i, count in groups.items():
            old = count - int(i == g)
            if old:
                jump = portal(couplings, i, g) * new_mass / loop
                total[i] += old*jump
                individual[i] += jump
        total[g] += injection
        individual[g] = individual[g].union(injection)


def auxiliary_step(args, direction, state):
    groups, quartic, gauge, yukawa, dt = args
    probe = {j: acb(v, direction[j]+arb(0, ".018")) for j, v in quartic.items()}
    rate = {j: v.imag for j, v in base["rates"](groups, probe, gauge, yukawa).items()}
    assert all(abs(dt*v) < q(".018") for v in rate.values())
    tube = {j: v+q(0).union(dt)*rate[j] for j, v in direction.items()}
    for j, v in rate.items():
        direction[j] += dt*v
    generator = coefficients(groups, quartic, gauge, yukawa)
    variation = coefficients(groups, tube, gauge, yukawa, gradient=True)
    n = len(groups)
    block = arb_mat(4*n, 4*n)
    for i in range(2*n):
        for j in range(2*n):
            block[i, j] = block[2*n+i, 2*n+j] = generator[i, j]
            block[2*n+i, j] = variation[i, j]
    vector = arb_mat([[v] for part in state for v in part])
    # The same enclosure also covers a step stopped at any intermediate time.
    within = (signed["advance_matrix"](block, q(0).union(dt))*vector).entries()
    end = (signed["advance_matrix"](block, dt)*vector).entries()
    for i, part in enumerate(state):
        part[:] = end[i*n:(i+1)*n]
    def frobenius(offset):
        return 2*sum(num*abs(within[offset+i]).upper()**2 for i, num in groups.items()).sqrt()
    return tensor_radius(groups, tube), frobenius(n), frobenius(3*n)


def propagate(previous, coefficient, source, duration):
    growth = (coefficient*duration.upper()).exp()
    return (previous*growth + source*(growth-1)/coefficient).upper()


def bounds():
    events = []
    gain = q(81, 977)**5
    zeta = q("-.009").union(q("-.001"))
    total, B, A, margin, raw, time = base["check"](gain, mass_path=events, bulk_portal=zeta)
    Y = (12*total).sqrt()
    heavy_y = arb(12).sqrt()*raw.sqrt()*gain*(6*total*time).exp()
    quartic_rate = 12*A+3*q(".020")+4*Y**2
    quartic_source = (8*Y*B+96*Y**3)*heavy_y
    quartic_max = propagate(q(0), quartic_rate, quartic_source, time)
    assert quartic_max < min(q(".020"), 6*margin)
    low = q("-1e6").union(q("1e6"))/unit
    state = ([low], [low], [q(0)], [q(0)])
    direction = {}
    elapsed = quartic_error = direction_error = mass_error = slope_error = q(0)
    for op, *args in events:
        if op == "threshold":
            add_at_threshold(*args, direction, state)
            continue
        R, M, N = auxiliary_step(args, direction, state)
        groups, _, _, _, dt = args
        if len(groups) == 3:
            elapsed += dt.upper()
            quartic_error = propagate(q(0), quartic_rate, quartic_source, arb.min(time.upper(), elapsed))
            mass_rate = 2*A+quartic_error+2*Y**2
            mass_force = quartic_error+4*Y*heavy_y
            direction_error = propagate(direction_error, 12*A+6*quartic_error+4*Y**2,
                                        (6*quartic_error+8*Y*heavy_y)*R, dt)
            mass_error = propagate(mass_error, mass_rate, mass_force*M, dt)
            slope_error = propagate(slope_error, mass_rate, mass_force*N+direction_error*M
                                    +(R+direction_error)*mass_error, dt)
    masses, slopes = state[1], state[3]
    heavy_lower = min(v.lower() for v in masses[1:])
    band = q(125)**2/(2*unit)
    assert heavy_lower-mass_error > band
    tilt = mass_error/(heavy_lower-band)
    heavy_slope = max(abs(v).upper() for v in slopes[1:])
    light_slope = state[2][0]
    lower = light_slope.lower()*(1-tilt**2)-heavy_slope*tilt**2-slope_error
    upper = light_slope.upper()+heavy_slope*tilt**2+slope_error
    print("mass error, slope error, auxiliary heavy floor, angle bound", *(
        v.str(8) for v in (mass_error, slope_error, heavy_lower, tilt)))
    print("auxiliary light slope", light_slope.str(8), "absolute heavy slope below",
          heavy_slope.str(8), "finite slope endpoints in 1e24 GeV²",
          lower.lower().str(9), upper.upper().str(9))
    assert lower > q('41.248') and upper < q('47.990') and tilt < q('.045967')
    print("weak-band width between", (2*band/upper).lower().str(8), (2*band/lower).upper().str(8))
    return lower.lower(), upper.upper()


if __name__ == "__main__":
    bounds()
