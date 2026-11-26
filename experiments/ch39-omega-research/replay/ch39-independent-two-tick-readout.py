"""Certify two distinct time bins using the radial partition of the lifted shell."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx, fmpq, fmpq_mat, fmpz_mat

local = run_path(str(Path(__file__).with_name("ch39-local-product-contact.py")))
contact = local["contact"]


def radial_reflection():
    sites, pole = local["lifted_sites"]()
    distances = [local["geometry"]["squared"](site, pole) for site in sites]
    nearest = min(distance for distance in distances if distance)
    indicator = fmpz_mat([[int(distance == nearest)] for distance in distances])
    total = len(sites)
    marked = int((indicator.transpose()*indicator)[0, 0])
    assert (nearest, marked, total-marked) == (68, 9, 16)
    # The unnormalized bin indicators have norms three and four.
    first = fmpq(2*marked-total, total)
    second = fmpq(2*3*4, total)
    assert first**2 + second**2 == 1

    laplacian, _ = local["lifted_laplacian"]()
    near = [i for i, distance in enumerate(distances) if distance == nearest]
    rest = [i for i, distance in enumerate(distances) if distance != nearest]
    assert sum(all(laplacian[i,j] == 0 for i in near) for j in rest) == 6
    walk = total*laplacian**0 - laplacian
    assert all(walk[i,j] >= 0 for i in range(total) for j in range(total))
    two_steps = walk**2
    generator = fmpq_mat(total, total)
    for i in near:
        generator[i,i] = first
    for j in rest:
        generator[j,j] = -first
        neighbors = [i for i in near if two_steps[i,j] > 0]
        assert neighbors
        for i in neighbors:
            generator[i,j] = generator[j,i] = second*3/(4*len(neighbors))
    near_state = fmpq_mat([[fmpq(int(i in near), 3)] for i in range(total)])
    rest_state = fmpq_mat([[fmpq(int(i in rest), 4)] for i in range(total)])
    assert generator*near_state == first*near_state + second*rest_state
    assert generator*rest_state == second*near_state - first*rest_state

    powers = [indicator]
    for _ in range(4):
        powers.append(laplacian*powers[-1])
    ones = fmpz_mat([[1] for _ in sites])
    coefficients = (1225, -1568, 402, -36, 1)
    filtered = sum((a*col for a, col in zip(coefficients, powers)), 0*indicator)
    assert total*filtered == marked*coefficients[0]*ones
    cubic_and_constant_vector = fmpz_mat([
        [col[i, 0] for col in powers[:4]] + [1] for i in range(total)
    ])
    assert cubic_and_constant_vector.rank() == 5
    return arb(first), arb(second)


def verify():
    ctx.prec = 192
    early, late = radial_reflection()
    one_tick, two_ticks = [], []
    flavor = contact["flavor"]
    coefficients = [arb(k)/2048 for k in flavor["codes"]]
    reference = flavor["native"]["plane_response"](acb.zeta_zero(11431).imag)
    reference_entries = (reference[0,0].real, reference[1,1].real, reference[0,1].real, reference[0,1].imag)
    for n, generator, weight, observables in contact["candidate_sources"](include_tangent=True):
        da, db, dc, dd = [x-y for x, y in zip(observables[5:9], reference_entries, strict=True)]
        distance = coefficients[0]**2*(da*da + db*db + 2*(dc*dc + dd*dd))
        observables.insert(-1, distance)
        outputs = contact["pulses"]["port_outputs"](generator, 2, impulse=True)
        first, delayed = [forward for forward, _ in outputs]
        one_tick.append((n, abs(first)**2, weight, observables))
        score = abs(early*first + late*delayed)**2
        two_ticks.append((n, score, weight, observables))

    total_weight = sum(weight for _, _, weight, _ in two_ticks)
    modes = local["graph_spectrum"]()
    for length, records in ((1, one_tick), (2, two_ticks)):
        rows = [(n, p, w/total_weight, obs) for n, p, w, obs in records]
        n, _, _, _, probability, _, means, residual = contact["pulse_expectations"](rows, modes, length)
        values = contact["mixture"]["readouts"](means)
        if length == 1:
            assert n == 11143 and values[0] < arb(".203")
            print("single time slot: first sine", values[0].str(10))
            continue
        shell = contact["flavor"]["shell"]
        tolerances = [shell["sine_radius"]]*3 + [shell["phase_radius"]]
        assert all(abs(value-center) < tolerance
                   for value, center, tolerance in zip(values, shell["centers"], tolerances, strict=True))
        assert n == 11431 and probability > arb(".999806")
        assert residual < arb("4.24e-7") < arb(2)**-21
        print("two time slots: label", n, "probability", probability.str(11))
        print("flavor", *(value.str(11) for value in values), "gauge residual <", residual.upper().str(8))
        radius_squared = contact["mixture"]["phase_radius_squared"](means)
        assert arb(".99998413") < radius_squared < arb(".99998414")
        assert 0 < means[9] < arb("2.452e-8")
        a, b, real, imag = means[5:9]
        response = acb_mat([[a, acb(real, imag)], [acb(real, -imag), b]])
        coherent = flavor["mixing"](response, coefficients)
        direct = contact["mixture"]["direct_observables"](coherent)
        restored = contact["mixture"]["readouts"](direct)
        assert all(abs(v-c)<r for v,c,r in zip(restored, shell["centers"], tolerances, strict=True))
        diamond_bound = (2*means[9]/64).upper()
        assert diamond_bound < arb("7.67e-10")
        uncertainty = arb(0, diamond_bound)
        physical = [value + uncertainty for value in direct[:4]] + [direct[4] + 2*uncertainty]
        finite_values = contact["mixture"]["readouts"](physical)
        assert all(abs(v-c)<r for v,c,r in zip(finite_values, shell["centers"], tolerances, strict=True))
        print("phase circle R^2", radius_squared.str(14), "variance bound", means[9].str(10))
        print("coherent flavor", *(v.str(11) for v in restored), "64-step diamond bound", diamond_bound.str(8))
        return response


if __name__ == "__main__":
    verify()
