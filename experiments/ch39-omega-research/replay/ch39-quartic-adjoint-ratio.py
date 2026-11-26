"""Charged scalar pseudo-determinants as exact polynomials in (u/v)^2."""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq_mat
from sympy import Poly, Rational, Symbol, factor_list

ROOT = Path(__file__).parent
MODEL = run_path(str(ROOT/"ch39-quartic-lock-thresholds.py"))
Q = Symbol("q")
a = MODEL["a"]


def product_polynomials():
    # The entries of J are affine in r. Three integer evaluations fix J.T J.
    blocks = [MODEL["mass_blocks"](r) for r in (0, 1, -1)]
    exponents = defaultdict(lambda: [Fraction(0)]*3)
    constants = [defaultdict(Fraction) for _ in range(3)]
    cache = {}
    for label, (base, metric) in blocks[0].items():
        _, weight = label
        charges = MODEL["weight_indices"](weight)
        if not any(charges):
            continue
        positive, negative = [block[label][0] for block in blocks[1:]]
        linear = (positive-negative)//2
        quadratic = (positive+negative)//2-base
        assert (base+linear+quadratic == positive).all()
        assert (base-linear+quadratic == negative).all()
        key = (base.shape, *(m.tobytes() for m in (base, linear, quadratic, metric)))
        if key not in cache:
            inverse = fmpq_mat(metric.tolist()).inv()

            def characteristic(r):
                return (inverse*fmpq_mat((base+r*linear+r*r*quadratic).tolist())).charpoly()

            polynomial = characteristic(1)
            nullity = next(i for i in range(len(polynomial)) if polynomial[i])
            rank = len(base)-nullity
            # The coefficient is even in r by the sign change on the adjoint.
            # Its degree in q is at most rank. Interpolate at the integer squares.
            result, basis = Poly(0, Q, domain="QQ"), Poly(1, Q, domain="QQ")
            for r in range(rank+1):
                value = characteristic(r)[nullity] * (-1)**rank
                exact = Rational(int(value.p), int(value.q))
                result += basis*((exact-result.eval(r*r))/basis.eval(r*r))
                basis *= Poly(Q-r*r, Q, domain="QQ")
            constant, factors = factor_list(result.as_expr())
            cache[key] = Fraction(constant), factors
        constant, factors = cache[key]
        for i, charge in enumerate(charges):
            constants[i][constant] += charge
        for factor, power in factors:
            for i, charge in enumerate(charges):
                exponents[factor][i] += charge*power
    return constants, dict(exponents)


def thresholds(products, ratio):
    constants, exponents = products
    q = ratio**2
    logs = []
    for i in range(3):
        total = sum(a(weight)*a(c).log() for c, weight in constants[i].items())
        for factor, indices in exponents.items():
            value = arb(0)
            for coefficient in Poly(factor, Q).all_coeffs():
                value = value*q+a(coefficient)
            total += a(indices[i])*value.log()
        logs.append(-total/12)
    return logs

def high_zero_readings():
    from flint import acb
    from sympy import fraction
    graph = run_path(str(ROOT/"ch39-weighted-arc-zero-phase.py"))
    numerator, denominator = fraction(graph["response"](3))
    evaluated = graph["evaluated"]
    for n in range(11485, 11491):
        zero = acb.zeta_zero(n)
        assert zero.real == a("1/2")
        height = zero.imag
        phase = 1-1/(2*(height**2+a("1/4")))
        yield n, height, evaluated(numerator, phase)/evaluated(denominator, phase)


def verify():
    from sympy import cancel, diff
    ctx.prec = 160
    products = product_polynomials()
    unit = thresholds(products, arb(1))
    reference = MODEL["scalar_thresholds"](MODEL["spectrum"]()[1])
    assert all((x-y).contains(0) for x, y in zip(unit, reference))
    _, powers = products
    assert powers == {Q: [20,8,11], Q+1: [0,0,3], Q+2: [10,42,28], 2*Q+1: [5,3,2]}
    weights = (115, -333, 218)
    exponents = {p: sum(w*b for w,b in zip(weights, index))/504 for p,index in powers.items()}
    exponents[Q+2] -= Fraction(1, 2)
    derivative = sum(c*diff(p,Q)/p for p,c in exponents.items())
    cubic = 1428*Q**3+562*Q**2-1435*Q-678
    assert cancel(derivative + cubic/(84*Q*(Q+1)*(Q+2)*(2*Q+1))) == 0
    assert Poly(cubic,Q).count_roots(0, 2) == 1
    assert cubic.subs(Q, Rational(10148,10000)**2) < 0 < cubic.subs(Q, Rational(10149,10000)**2)

    readings = [epsilon for _, _, epsilon in high_zero_readings()]
    assert all(a("4.6086e-5") < e < a("4.6087e-5") for e in readings)
    for r in (a("11/25"), a(2)):
        scalar = thresholds(products, r)
        for eps in [a(Fraction(6, 19**4)), *readings]:
            low, v, log_x = MODEL["matched"](scalar, eps, arb(1), r)
            x = log_x.exp()
            assert a(".63") < low < a(".78") and low**2/arb.pi() < a("1/4")
            assert x**2 < 4*arb.pi() and v*eps*x > 1000
        print("u/v", r.str(6), "final high-zero sample", low.str(9), v.str(9), x.str(9))
    print("Exact threshold factors", powers)
    print("Logarithmic powers for g_low", exponents)


if __name__ == "__main__":
    verify()
