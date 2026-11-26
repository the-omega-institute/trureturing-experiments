"""Validated one-loop right-handed flavor directions in two RH-tagged matching windows."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx

ctx.prec = 256
ckm = run_path(str(Path(__file__).with_name('ch39-flavour-symmetry-gap.py')))['ckm']
mixing = acb_mat(ckm(*[arb(t).sin() for t in ('.22704', '.04275', '.00377')], arb('1.139')))
gauge0 = list(map(arb, ('.467453', '.63811', '1.0583')))
up = list(map(arb, ('.00000615', '.00311', '.8616')))
down = list(map(arb, ('.0000135', '.000268', '.01401')))
leptons = list(map(arb, ('.0000028227', '.00059465', '.010101')))
gauge_terms = [tuple(map(arb, numbers)) for numbers in (
    ('.85', '2.25', '8'), ('.25', '2.25', '8'), ('2.25', '2.25', '0'))]
beta = [arb(41)/10, -arb(19)/6, arb(-7)]
loop = 16*arb.pi()**2


def diagonal(values):
    return acb_mat([[v if i == j else 0 for j in range(3)] for i, v in enumerate(values)])


def adjoint(matrix):
    return matrix.conjugate().transpose()


cross = diagonal(up) * mixing * diagonal(down)
up_cross, down_cross = cross*adjoint(cross), adjoint(cross)*cross
initial = (diagonal([a*a for a in up]), diagonal([a*a for a in down]),
           *[a*a for a in leptons], arb(1))


def coefficients(state, time, order):
    # Power series for the two right Gram matrices, three charged-lepton
    # squares and the multiplier of Z(0)^* Z(0), relative to the current time.
    u, d, e, mu, tau, scale = ([x] for x in state)
    lepton = (e, mu, tau)
    au, ad, ae = [], [], []
    gauge = [1/(1/g**2 - 2*b*time) for g, b in zip(gauge0, beta, strict=True)]
    ratio = [2*b*g for b, g in zip(beta, gauge, strict=True)]
    for j in range(order):
        trace = 3*(u[j].trace().real + d[j].trace().real) + sum(l[j] for l in lepton)
        for line, weights in zip((au, ad, ae), gauge_terms, strict=True):
            line.append(trace - sum(a*b for a, b in zip(weights, gauge, strict=True)))
        for gram, scalar, constant in ((u, au, up_cross), (d, ad, down_cross)):
            value = -3*scale[j]*constant
            for i in range(j + 1):
                value += 2*scalar[i]*gram[j-i] + 3*(gram[i]*gram[j-i])
            gram.append(value / (j+1))
        for l in lepton:
            l.append(sum((2*ae[i] + 3*l[i])*l[j-i] for i in range(j+1)) / (j+1))
        scale.append(2*sum((au[i]+ad[i])*scale[j-i] for i in range(j+1)) / (j+1))
        gauge = [a*b for a, b in zip(gauge, ratio, strict=True)]
    return u, d, e, mu, tau, scale


def widen(value, radius):
    if isinstance(value, acb_mat):
        ball = acb(arb(0, radius), arb(0, radius))
        return acb_mat([[v + ball for v in row] for row in value.tolist()])
    return value + arb(0, radius)


def entry_bound(value):
    if isinstance(value, acb_mat):
        return max(v.abs_upper() for v in value.entries())
    return value.abs_upper()


def advance(state, time, step):
    radii = tuple(map(arb, ('.35', '.01', '.001', '.001', '.001', '.6')))
    tube = [widen(x, r) for x, r in zip(state, radii, strict=True)]
    remainder = coefficients(tube, time + arb(0).union(step), 23)
    for coeff, r in zip(remainder, radii, strict=True):
        assert entry_bound(coeff[1])*step.upper() < r
    result = []
    for series, envelope in zip(coefficients(state, time, 22), remainder, strict=True):
        value = series[-1] + step*envelope[-1]
        for term in reversed(series[:-1]):
            value = term + step*value
        result.append(value)
    return result


def mass_lines(gram):
    complements = [[j for j in range(3) if j != i] for i in range(3)]
    disks = [gram[i, i].real + arb(0, sum(gram[i, j].abs_upper() for j in cols).upper())
             for i, cols in enumerate(complements)]
    assert all(a.upper() < b.lower() for a, b in zip(disks, disks[1:]))
    result = []
    for i, cols in enumerate(complements):
        eigenvalue = disks[i]
        rhs = acb_mat([[-gram[j, i]] for j in cols])
        for _ in range(2):
            corner = acb_mat([[gram[j, k] - (eigenvalue if j == k else 0)
                               for k in cols] for j in cols])
            tail = corner.solve(rhs)
            eigenvalue = (gram[i, i] + sum(gram[i, j]*tail[k, 0]
                                          for k, j in enumerate(cols))).real
        vec = acb_mat([[1 if j == i else tail[cols.index(j), 0]] for j in range(3)])
        vec = vec / sum(abs(v)**2 for v in vec.entries()).sqrt()
        result.append((eigenvalue.sqrt(), vec))
    return result


def matching_data(state):
    up_lines, down_lines = mass_lines(state[0]), mass_lines(state[1])
    light, bottom = down_lines[0], down_lines[2]
    moment = (adjoint(light[1])*down_cross*light[1])[0, 0].real / down_cross.trace().real
    ckm_column = [state[-1].sqrt()*abs((adjoint(u)*cross*bottom[1])[0, 0])/(y*bottom[0])
                  for y, u in up_lines]
    electron, muon, tau = [v.sqrt() for v in state[2:5]]
    scale = (bottom[0] + light[0])/(muon + tau)
    down_cost = bottom[0] - scale*muon
    assert abs(down_lines[1][0] - scale*electron) < down_cost
    up_bound = ((up_lines[1][0]*ckm_column[1])**2 + (up_lines[2][0]*ckm_column[2])**2).sqrt()
    up_bound -= up_lines[0][0]*(1 - ckm_column[0]**2).sqrt()
    assert up_bound / up_lines[2][0] > arb('.9988')
    assert down_cost / bottom[0] > arb('.9443')
    return light[1], moment


def angle_sine(a, b):
    wedge = sum(abs(a[i, 0]*b[j, 0] - a[j, 0]*b[i, 0])**2 for i in range(3) for j in range(i))
    return (wedge / (sum(abs(v)**2 for v in a.entries())*sum(abs(v)**2 for v in b.entries()))).sqrt()


def within(values, lower, upper):
    assert min(x.lower() for x in values) > arb(lower)
    assert max(x.upper() for x in values) < arb(upper)


def main():
    state, start = initial, arb(0)
    first_line = acb_mat([[1], [0], [0]])
    windows = [
        ('7.41e10', '7.66e10', ('9.7511e-11', '9.7548e-11'), ('3.4231e-6', '3.4276e-6')),
        ('1.80e11', '1.861e11', ('9.8444e-11', '9.8480e-11'), ('3.5378e-6', '3.5422e-6')),
    ]
    directions = []
    for lo, hi, moment_bounds, angle_bounds in windows:
        end, right = [(arb(v)/1000).log()/loop for v in (lo, hi)]
        count = int(((end-start)*1200).ceil().unique_fmpz())
        step = (end-start)/count
        for i in range(count):
            state = advance(state, start+i*step, step)
        vectors, moments = [], []
        width = (right-end)/24
        for j in range(24):
            time = end+j*width
            vec, moment = matching_data(advance(state, time, arb(0).union(width)))
            vectors.append(vec)
            moments.append(moment)
            state = advance(state, time, width)
        start = right
        directions.append(vectors)
        within(moments, *moment_bounds)
        within([angle_sine(first_line, v) for v in vectors], *angle_bounds)
        print(lo, hi, 'invariant expectation', moment_bounds, 'angle sine', angle_bounds)
    gaps = [angle_sine(a, b) for a in directions[0] for b in directions[1]]
    within(gaps, '1.1026e-7', '1.1902e-7')
    print('Across both windows, sine of angular separation in (1.1026e-7, 1.1902e-7).')


if __name__ == '__main__':
    main()
