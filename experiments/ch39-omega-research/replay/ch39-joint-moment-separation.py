"""Enclose both polarization readouts along the rematched Planck center."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_series, ctx, fmpq_mat
from sympy import Poly, Rational, cancel, symbols

SOURCE = run_path(str(Path(__file__).with_name("ch39-softness-polarization.py")))
MODEL = SOURCE["RATIO"]["MODEL"]
a = SOURCE["a"]
Q = symbols("q")
WINDOW = F(1, 40000)


def scalar_remainders():
    """Exact inverse-power sums after removing the unsoftened weak triplet."""
    tables = [MODEL["mass_blocks"](r) for r in (0, 1, -1)]
    cache, sums = {}, [0, 0, 0, 0]
    for label, (base, metric) in tables[0].items():
        weight = sum(w*c for w, c in zip(MODEL["weight_indices"](label[1]), (5, 3, -8)))
        if not weight:
            continue
        pos, neg = (table[label][0] for table in tables[1:])
        key = tuple(m.tobytes() for m in (base, pos, neg, metric))
        if key not in cache:
            inv = fmpq_mat(metric.tolist()).inv()
            matrices = [inv*fmpq_mat(m.tolist()) for m in (base, (pos-neg)//2, (pos+neg)//2-base)]
            def characteristic(r):
                return (matrices[0]+r*matrices[1]+(r*r)*matrices[2]).charpoly()
            size = len(base)
            nullity = size-fmpq_mat(pos.tolist()).rank()
            coeff = [Poly(0, Q, domain="QQ") for _ in range(5)]
            basis = Poly(1, Q, domain="QQ")
            for r in range(size+1):
                poly = characteristic(r)
                assert poly == characteristic(-r)
                assert all(poly[j] == 0 for j in range(nullity))
                for j in range(5):
                    val = poly[nullity+j]
                    exact = Rational(int(val.p), int(val.q))
                    coeff[j] += basis*((exact-coeff[j].eval(r*r))/basis.eval(r*r))
                basis *= Poly(Q-r*r, Q, domain="QQ")
            c0, c1, c2, c3, c4 = (p.as_expr() for p in coeff)
            assert all(x >= 0 for x in Poly((-1)**(size-nullity)*c0, Q).all_coeffs()) and c0 != 0
            b, c, d, e = [x/c0 for x in (c1, c2, c3, c4)]
            cache[key] = [cancel(x) for x in (-b, b*b-2*c, -b**3+3*b*c-3*d,
                                                b**4-4*b*b*c+2*c*c+4*b*d-4*e)]
        for j, power in enumerate(cache[key]):
            sums[j] += weight*power
    functions = []
    for j, total in enumerate(sums, 1):
        # Remove the old weak and color terms; insert their physical masses below.
        num, den = cancel(total-6/Q**j+24/(1+Q)**j).as_numer_denom()
        functions.append([[a(v) for v in Poly(p, Q).all_coeffs()] for p in (num, den)])
    return functions


def polynomial(coefficients, x):
    result = coefficients[0]
    for coefficient in coefficients[1:]:
        result = result*x+coefficient
    return result


def check_scalar_reduction(functions, blocks):
    # Test all four powers directly by taking -p'(x)/p(x) after removing null modes.
    r = a("1/6")
    q = r*r
    sums = [arb(0) for _ in functions]
    for weight, (base, linear, quadratic), nullity in blocks:
        poly = (base+r*linear+q*quadratic).charpoly()
        reduced = [poly[nullity+j] for j in range(5)]
        deriv = arb_series([(j+1)*reduced[j+1] for j in range(4)], prec=4)
        logarithmic = -deriv/arb_series(reduced, prec=4)
        for i in range(4):
            sums[i] += weight*logarithmic[i]
    for j, (total, (num, den)) in enumerate(zip(sums, functions), 1):
        restored = polynomial(num, q)/polynomial(den, q)+6/q**j-24/(1+q)**j
        assert total.overlaps(restored)


def scaled_readouts(q, sigma, gamma, low2, functions, lam):
    # Multiply M_j by (lambda*sigma*q)^j before normalizing to reduce dependency.
    t = sigma*q
    unit = lam*t
    k = a("8/5")
    z, w = low2*(1+3/k**2), low2*(3+k**2)
    trace, det = (3+q)*z+w, (2+q)*z*w
    # Sums of powers of the two scaled inverse vector masses squared.
    first, product = unit*trace/det, unit**2/det
    vector_powers = (first, first**2-2*product, first**3-3*first*product,
                     first**4-4*first**2*product+2*product**2)
    ell = ((gamma/624)**2+a("1/4")).sqrt()
    real = ((ell+a("1/2"))/2).sqrt()
    d2, h2 = a("13/10")*q*z, a("7/2")*q*z
    plus, minus = ell+1+2*real, ell+1-2*real
    masses = (h2*plus, h2*ell, d2*minus, d2*ell)
    moments = []
    for j, ((num, den), pair) in enumerate(zip(functions, vector_powers), 1):
        scalar = (6 + t**j*polynomial(num, q)/polynomial(den, q)
                  -24*(t/(1+t))**j-6*unit**j)
        vector = 18*(pair-(unit/(3*z+w))**j)
        fs, fv = SOURCE["BOSON"]["spin_weights"](j)
        fermion = sum(sign*(unit/mass)**j for sign, mass in zip((1,-1,1,-1), masses))
        moments.append((fs*scalar+fv*vector)/3+fermion)
    assert moments[0][0] > 0
    return tuple(m/moments[0]**j for j, m in enumerate(moments[1:], 2))


class Center:
    def __init__(self, gamma, epsilon, products, powers, functions):
        self.gamma, self.functions = gamma, functions
        assert tuple(42*p for p in powers) == (-327, -107, 1107, -1)
        scalar = SOURCE["MATCH"]["SCATTER"]["scalar_thresholds"](products, arb(1), a("1/12"))
        low, v, _ = MODEL["matched"](scalar, epsilon, arb(1))
        geometry = (2*arb.pi()*(gamma**2+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        self.anchor = v/low*geometry
        self.low = low
        weights = (115, -333, 218)
        factors = SOURCE["RATIO"]["Q"]
        exponents = {p: sum(w*c for w, c in zip(weights, charges))/504
                     for p, charges in products[1].items()}
        exponents[factors+2] -= F(1, 2)
        powers84 = tuple(84*exponents[p] for p in (factors, factors+1, factors+2, 2*factors+1))
        assert powers84 == (339, 109, -1164, 2)
        self.target42 = (SOURCE["PLANCK"]["TARGET"]/self.anchor)**42

    @staticmethod
    def mass_partials(q, sigma):
        # Logarithmic derivatives of the positive 42nd-power mass, at fixed eta.
        return (-327/q-107*sigma/(1+sigma*q)+1107/(q+2)-2/(2*q+1),
                108/sigma-107*q/(1+sigma*q))

    def root(self, delta, eta=0, refine=False):
        sigma = (1+a(delta))/256
        def residual(r):
            q = r*r
            lhs = (1+a(eta))*sigma**108*((q+2)/3)**1107
            rhs = self.target42*q**327*((1+sigma*q)/2)**107*((2*q+1)/3)
            return lhs-rhs
        lo, hi = SOURCE["PLANCK"]["enclose_root"](residual, 0, F(".16"), F(".17"), False)
        if refine:
            for _ in range(46):
                mid = (lo+hi)/2
                value = residual(a(mid))
                assert value > 0 or value < 0
                if value > 0:
                    lo = mid
                else:
                    hi = mid
        return a(lo).union(a(hi))**2

    def readouts(self, delta, q, relative_lambda=arb(1)):
        sigma = (1+delta)/256
        factors = (q, (1+sigma*q)/2, (q+2)/3, (2*q+1)/3, sigma, relative_lambda)
        factors = [arb_series([x], prec=2) if isinstance(x, arb) else x for x in factors]
        powers = (339, 109, -1164, 2, -111, -2)
        power42 = arb(1)
        logarithmic_slope = arb(0)
        for factor, p in zip(factors, powers):
            assert factor[0] > 0
            power42 *= factor[0]**p
            logarithmic_slope += p*factor[1]/factor[0]
        low2 = self.low**2*power42.root(42)*arb_series([1, logarithmic_slope/42], prec=2)
        return scaled_readouts(q, sigma, self.gamma, low2, self.functions, relative_lambda/12)

    def slopes(self, left, right, qleft, qright):
        delta = arb_series([a(left).union(a(right)), arb(1)], prec=2)
        q = qleft.union(qright)
        sigma = (1+delta)/256
        partial_q, partial_sigma = self.mass_partials(q, sigma[0])
        assert partial_sigma > 0 and partial_q < 0
        dq = -partial_sigma/(256*partial_q)
        return self.readouts(delta, arb_series([q, dq], prec=2))


def box_jacobian(center, delta, eta, q):
    sigma = (1+delta)/256
    partial_q, partial_sigma = center.mass_partials(q, sigma)
    partial_delta = partial_sigma/256
    assert partial_q < 0 and partial_delta > 0
    columns = []
    for d, e, partial in ((1, 0, partial_delta), (0, 1, 1/(1+eta))):
        qjet = arb_series([q, -partial/partial_q], prec=2)
        values = center.readouts(arb_series([delta, d], prec=2), qjet,
                                 arb_series([1+eta, e], prec=2))
        columns.append([x[1] for x in values[:2]])
    jac = arb_mat(columns).transpose()
    return arb_mat([[x[0]] for x in values]), jac


def shared_pair(center, targets):
    # Round only the search centers. The final inclusion test verifies them.
    p = [F(0), F(0)]
    for _ in range(3):
        q = center.root(*p, refine=True)
        values, jac = box_jacobian(center, *map(a, p), q)
        residual = arb_mat([[values[i, 0]-targets[i]] for i in range(2)])
        step = jac.solve(residual)
        for i in range(2):
            p[i] -= F(round(float(step[i, 0])*10**19), 10**19)
    radius = F(1, 10**18)
    assert abs(p[0])+radius < F(1, 2*10**6) and abs(p[1])+radius < F(118, 10**6)
    corners = [[v+side*radius for v in p] for side in (-1, 1)]
    qbox = center.root(*corners[0], refine=True).union(center.root(*corners[1], refine=True))
    delta, eta = [a(lo).union(a(hi)) for lo, hi in zip(*corners)]
    box_values, box_jac = box_jacobian(center, delta, eta, qbox)
    q = center.root(*p, refine=True)
    values, jac = box_jacobian(center, *map(a, p), q)
    # Freeze the midpoint of the preconditioner as an exact dyadic matrix.
    inv = jac.inv()
    B = arb_mat([[inv[i, j].mid() for j in range(2)] for i in range(2)])
    assert abs(B.det()) > 0
    identity = arb_mat([[1, 0], [0, 1]])
    # The center root increases in both coordinates throughout this rectangle.
    bounds = [WINDOW, F(120, 10**6)]
    lower = [-v for v in bounds]
    full_q = center.root(*lower).union(center.root(*bounds))
    rectangle = [a(lo).union(a(hi)) for lo, hi in zip(lower, bounds)]
    _, full_jac = box_jacobian(center, *rectangle, full_q)
    global_diff = identity-B*full_jac
    assert all(sum(abs(global_diff[i, j]) for j in range(2)) < a("1/4") for i in range(2))
    diff = identity-B*box_jac
    contraction = max(sum(abs(diff[i, j]) for j in range(2)).upper() for i in range(2))
    residual = B*arb_mat([[values[i, 0]-targets[i]] for i in range(2)])
    assert contraction < a(".01")
    assert all(abs(x)+contraction*a(radius) < a(radius) for x in residual.entries())
    print("matched relative coefficients, ppm", [float(v*10**6) for v in p],
          "fourth moment", box_values[2, 0].str(18))
    return box_values[2, 0]


def verify():
    ctx.prec = 180
    ratio = SOURCE["RATIO"]
    products = ratio["product_polynomials"]()
    powers = SOURCE["PLANCK"]["mass_exponents"](products)
    functions = scalar_remainders()
    domain = a(".16").union(a(".17"))**2
    sigmas = (1+a(-WINDOW).union(a(WINDOW)))/256
    coefficient, target = a("9.276"), a("5.36647005")
    previous = None
    fourths = []
    blocks = list(SOURCE["BOSON"]["charged_blocks"]())
    check_scalar_reduction(functions, blocks)
    for n, gamma, epsilon in ratio["high_zero_readings"]():
        center = Center(gamma, epsilon, products, powers, functions)
        # These signs on the large domain justify bounding q by its endpoint roots.
        partial_q, partial_sigma = center.mass_partials(domain, sigmas)
        assert partial_q < 0 and partial_sigma > 0
        edges = [WINDOW*F(i, 32) for i in range(-32, 33)]
        roots = [center.root(edge) for edge in edges]
        zero = arb_series([arb(0)], prec=2)
        def point(i):
            return center.readouts(zero+a(edges[i]), zero+roots[i])
        H, J, _ = (v[0] for v in point(32))
        # Cross-check the polynomial reduction against the original matrix traces.
        reference = SOURCE["center_function"](gamma, epsilon, products, powers, blocks)
        h_ref, j_ref = reference(0)
        assert H.overlaps(h_ref) and J.overlaps(j_ref)
        assert point(0)[0][0] < target < point(-1)[0][0]
        for left, right, qleft, qright in zip(edges, edges[1:], roots, roots[1:]):
            hj, jj, _ = center.slopes(left, right, qleft, qright)
            assert a(".03711") < hj[1] < a(".03713")
            assert abs(jj[1]-coefficient*hj[1]) < a(".000030")
        drift = a(WINDOW)*a(".000030")
        projected = J-coefficient*(H-target)+arb(0, drift)
        if previous is not None:
            assert projected-previous > a("5.3e-8")
        previous = projected
        fourths.append(shared_pair(center, [target, a("32.30328974")]))
        print(n, "projected readout", projected.str(15),
              "H center", H.str(15), "J center", J.str(15))

    gaps = [v-u for u, v in zip(fourths, fourths[1:])]
    assert all(gap > a("4.31e-12") for gap in gaps)
    print("fourth moment gaps", [gap.str(10) for gap in gaps])


if __name__ == "__main__":
    verify()
