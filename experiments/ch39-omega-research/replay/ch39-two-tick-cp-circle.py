"""Continuous CKM phase coverage for a fixed two-tick ground state and two shell arcs."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

two_ticks = run_path(str(Path(__file__).with_name("ch39-independent-two-tick-readout.py")))
flavor = two_ticks["contact"]["flavor"]
fit = flavor["arc"]["fit"]
shell = flavor["shell"]
identity = arb_mat([[1, 0], [0, 1]])


def modulus_equations(target):
    _, terms = shell["rephasing_equations"](target)
    p = fit["p"]
    constants = [sum((p[j]**2-p[i]**2)*abs(row[j])**2 for j in range(3) if i != j)
                 for i, row in enumerate(target[:2])]
    return constants, terms


def root_box(point_parameters, parameters, seed, domain):
    point = list(map(arb, seed))
    eq = fit["equations"]
    for _ in range(5):
        value, jacobian = eq(*point, point_parameters)
        step = jacobian.solve(value)
        point = [(pos-step[i,0]).mid() for i, pos in enumerate(point)]
    _, jacobian = eq(*point, point_parameters)
    inverse = jacobian.inv()
    residual = inverse*eq(*point, parameters)[0]
    err = max(abs(x).upper() for x in residual.entries())
    radius = 2*err + arb("1e-18")
    box = fit["around"](point, radius)
    remainder = identity-inverse*eq(*box, parameters)[1]
    lipschitz = max(sum(abs(remainder[i,j]).upper() for j in range(2)) for i in range(2))
    assert err + lipschitz*radius < radius
    assert all(a.contains_interior(b) for a,b in zip(domain,box))
    # Injectivity on a common rectangle identifies the same root across adjacent phase tiles.
    global_remainder = identity-inverse*eq(*domain, parameters)[1]
    assert all(sum(abs(global_remainder[i,j]) for j in range(2)) < 1 for i in range(2))
    return box


def coordinates(matrix):
    return [matrix[0,0].real, matrix[1,1].real, matrix[0,1].real, matrix[0,1].imag]


def coefficients_and_gap(target, root, inverse):
    basis = flavor["basis"]
    # Expand the row rephasing into phase differences so the large diagonal
    # terms have exactly zero phase width.
    angles = [*root, arb(0)]
    p = fit["p"]
    full = acb_mat([[target[i][j] * sum(
        acb(target[i][k]).conjugate()*p[k]/p[i] *
        (1 if j == k else fit["cis"](angles[j]-angles[k])) for k in range(3))
        for j in range(3)] for i in range(3)])
    w = basis.transpose()*full*basis
    trace = w.trace()
    gap = (abs(trace)/2).acos()
    assert 0 < gap < arb.pi()/2
    reduced = w*(trace.conjugate()/abs(trace))
    hermitian = (reduced-reduced.conjugate().transpose())*(
        gap/acb(0, 2*gap.sin()))
    hermitian += (trace.arg()+2*arb.pi())*flavor["identity"]
    entries = coordinates(hermitian)
    minor = (basis[0,0]*basis[2,1]-basis[0,1]*basis[2,0]).real
    skew = (full[0,2]-full[2,0])*trace.conjugate()/abs(trace)
    entries[3] = -gap*skew.real/(2*gap.sin()*minor)
    return (inverse*arb_mat([[v] for v in entries])).entries(), gap


def dyadic_accuracy(response, sines):
    # Round the four coefficients to the nearest multiple of 2^-18.
    norm = sum(abs(z)**2 for z in response.entries()).sqrt()
    error = ((norm+3)/2**19).upper()
    assert error < arb("6.31e-6")
    s,t,u = sines
    a,b,c = [(1-x*x).sqrt() for x in sines]
    umax = u.upper()+error
    cos_min = (1-umax**2).sqrt()
    sine_error = error*(1/cos_min+umax/cos_min**3)
    assert sine_error < arb("6.34e-6")
    # Use a rephasing-invariant complex number whose argument is the physical phase.
    bounds = [arb(1), umax, (s*b+a*t*u).upper()+error, (t*c).upper()+error]
    product_error = error*sum(prod(bounds[:i]+bounds[i+1:]) for i in range(4))
    cb = bounds[-1]
    ratio = umax**2/(1-umax**2)
    correction_error = error*((2*cb**2+2*cb)*ratio+cb**2*2*umax/(1-umax**2)**2)
    amplitude = (a*b*c*c*s*t*u).lower()
    phase_error = ((product_error+correction_error)/amplitude).asin()
    assert phase_error < arb(".0021")
    print("2^-18 coefficient mesh: sine error < 6.34e-6 and circular phase error < 0.0021.")


def verify(response):
    matrices = (response, flavor["arcs"][32], flavor["arcs"][152], flavor["identity"])
    inverse = arb_mat([coordinates(m) for m in matrices]).transpose().inv()
    centers = shell["centers"][:3]
    sines = [x+arb(0,shell["sine_radius"].upper()) for x in centers]
    divisions = 512
    envelopes = None
    domain = [arb("2.4 +/- .23"), arb("4.55 +/- .07")]
    for i in range(divisions):
        phase = (2*i+1-divisions)*arb.pi()/divisions
        box = phase + arb(0, (arb.pi()/divisions).upper())
        center_target = shell["standard_matrix"](*centers, phase)
        target = shell["standard_matrix"](*sines, box)
        root = root_box(modulus_equations(center_target), modulus_equations(target),
                        ("2.4","4.55"), domain)
        values, _ = coefficients_and_gap(target, root, inverse)
        envelopes = values if envelopes is None else [x.union(y) for x,y in zip(envelopes, values)]
    limits = ((".16","2.8"), (".5","3.7"), (".93","3.6"), ("2","5.2"))
    assert all(arb(lo) < v < arb(hi) for v,(lo,hi) in zip(envelopes,limits,strict=True))
    print("512 overlapping closed phase tiles; all four coefficients strictly positive.")
    dyadic_accuracy(response, sines)
    return inverse


def exclude_gain_window(inverse, gain_window):
    centers = shell["centers"]
    widths = [shell["sine_radius"]]*3 + [shell["phase_radius"]]
    domains = [fit["around"](list(map(arb, seed)), arb(".045")) for seed in fit["seeds"]]
    full_parameters = modulus_equations(shell["target_box"]())
    pending = [(0, 0, 0)]
    while pending:
        level, i, j = pending.pop()
        unit = arb.pi()/(2**level)
        box = fit["around"]([(2*i+1)*unit, (2*j+1)*unit], unit)
        if any(all(a.contains_interior(b) for a,b in zip(domain,box)) for domain in domains):
            continue
        values, _ = fit["equations"](*box, full_parameters)
        if all(0 in v for v in values.entries()):
            assert level < 17
            pending.extend((level+1, 2*i+a, 2*j+b) for a in (0,1) for b in (0,1))

    stack = [(list(centers), widths, 0)]
    accepted = 0
    while stack:
        center, half, level = stack.pop()
        target = shell["standard_matrix"](*(v+arb(0,r.upper()) for v,r in zip(center,half)))
        center_parameters = modulus_equations(shell["standard_matrix"](*center))
        parameters = modulus_equations(target)
        for seed, domain in zip(fit["seeds"], domains):
            xy = root_box(center_parameters, parameters, seed, domain)
            pars, gap = coefficients_and_gap(target, xy, inverse)
            winding = (gain_window/pars[0]-1)*gap/arb.pi()
            if winding.contains_integer():
                break
        else:
            accepted += 1
            continue
        assert level < 24
        axis = (2, 3, 0, 1)[level%4]
        narrow = half.copy()
        narrow[axis] /= 2
        for sign in (-1, 1):
            child = center.copy()
            child[axis] += sign*narrow[axis]
            stack.append((child, narrow, level+1))
    print("Physical angle/phase box:", accepted, "verified tiles; all logarithms excluded for", gain_window.str(5))


if __name__ == "__main__":
    exclude_gain_window(verify(two_ticks["verify"]()), arb("1 +/- .05"))
