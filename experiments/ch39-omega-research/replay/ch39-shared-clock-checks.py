"""复算共同算术钟的第七项、三段谱权和十九点中性壳的刚度响应。"""
import json
from fractions import Fraction
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx, fmpq, fmpq_poly, fmpz_poly


ROOT = Path(__file__).resolve().parent
SIXTH_SHIFT = "1.360900010e-8"
SCALE = 2_903_600_000
HEIGHTS = (38 * 10**9, 45 * 10**9, 70 * 10**9, 10**11, 10**12)


def rational(decimal):
    number = Fraction(decimal)
    return fmpq(number.numerator, number.denominator)


def held_out(extrapolation, reference, moments, rules):
    lower = extrapolation["readout"](rules[1], 7)
    upper = extrapolation["readout"](rules[0], 7)
    assert arb(".93341795164287258400188") < lower < arb(".93341795164287258400189")
    assert arb(".93341795164501234976009") < upper < arb(".93341795164501234976010")
    assert arb("1.20207618035599e-7") < lower-reference < upper-reference < arb("1.20209757801358e-7")
    for value, code in ((lower, 32071996606), (upper, 32071996606), (reference, 32071992476)):
        assert (2**35 * value).floor().unique_fmpz() == code
    codes_36 = [(2**36 * value).floor().unique_fmpz() for value in (lower, upper)]
    assert codes_36 == [64143993212, 64143993213]
    hankel = arb_mat([[moments[i+j] for j in range(4)] for i in range(4)])
    assert hankel.det() > arb("8.95e-34")
    print("第七项径向两端", lower.str(23), upper.str(23))
    print("35 位预报与参考码", 32071996606, 32071992476, "36 位径向端点码", *codes_36)
    print("四阶矩行列式", hankel.det().str(14))


def height_rule(height, bands, moments, denominator):
    cutoff = 1/(1+(arb(height)/SCALE)**2)
    nodes, weights = bands["cutoff_rule"](cutoff, moments, denominator)
    for upper in (False, True):
        bands["verify_indicator"](nodes, denominator, upper)
    return cutoff, nodes, weights


def band_penalty(radius):
    y = fmpq_poly([0, 1])
    rates = [j*j-1 for j in range(2, 7)]
    factors = [1+rate*y for rate in rates]
    denominator = prod(factors)
    a, b = map(rational, (".001212594", ".004902081"))
    shape = y*(y-a)**2*(y-b)**2*(1-y+rational("4e-9")*y)
    low_scale, high_scale = SCALE*(1-radius), SCALE*(1+radius)

    def cutoff(height, scale):
        return scale**2/(scale**2+height**2)

    h1, h2, h3, h4, h5 = HEIGHTS
    t = cutoff(h5, low_scale)
    numerator = shape*(denominator(t)/shape(t))
    gap = numerator-denominator*(1+18*y)/(1+18*t)
    intervals = (
        (t, cutoff(h4, high_scale)),
        (cutoff(h3, low_scale), cutoff(h2, high_scale)),
        (cutoff(h1, low_scale), fmpq(1)),
    )
    for left, right in intervals:
        assert t <= left < right
        powers = gap(fmpq_poly([left, right-left]))
        coefficients = [sum(powers[i]*comb(k, i)/comb(6, i) for i in range(k+1)) for k in range(7)]
        assert all(c >= 0 for c in coefficients)

    constant = numerator(1)/denominator(1)
    time_coefficients = []
    for i, rate in enumerate(rates):
        pole = -fmpq(1, rate)
        remaining = prod(factors[:i]+factors[i+1:])
        time_coefficients.append(numerator(pole)/((1-pole)*remaining(pole)))
    time_coefficients.insert(0, numerator(0)-constant-sum(time_coefficients))
    basis = [(1-y)*denominator, *[(1-y)*prod(factors[:i]+factors[i+1:]) for i in range(5)]]
    assert numerator == constant*denominator+sum(c*p for c, p in zip(time_coefficients, basis))
    return t, constant, time_coefficients


def robust_bands():
    data = json.loads((ROOT / "ch39-minimum-tick-finite-inputs.json").read_text())
    assert data["depth"] == 288
    cell = fmpq(1, 2**data["depth"])
    target_floors = [cell*level for level in data["levels"][1:7]]
    target_floors[5] += rational(SIXTH_SHIFT)
    for radius, cases in (
        (fmpq(0), (
            ("精确六项", "0", "0", ".99931149998", ".99931139549"),
            ("固定尺度六项误差", "1.67e-15", "1.67e-15", ".99915077450", ".99915077450"),
            ("仅第六项误差", "0", "2.10e-13", ".99915043084", ".99915043084"),
        )),
        (fmpq(1, 25), (
            ("尺度与六项同时变动", Fraction(1, 2**50), Fraction(1, 2**50), ".99916144733", ".99916133005"),
        )),
    ):
        t, constant, coefficients = band_penalty(radius)
        pinned_upper = constant+sum(c*v+max(c, 0)*cell for c, v in zip(coefficients, target_floors))
        for label, first_five, sixth, arithmetic, flavor in cases:
            upper = pinned_upper+sum(abs(c) for c in coefficients[:5])*rational(first_five)
            upper += abs(coefficients[5])*rational(sixth)
            assert 1-upper > rational(arithmetic)
            assert 1-(1+18*t)*upper > rational(flavor)
            print(label, "三段算术和味概率下界", arithmetic, flavor)


def stiffness(neutral, denominator, rules, band_masses):
    polynomial = neutral["neutral_polynomial"]()
    exact = fmpz_poly([int(c.unique_fmpz()) for c in polynomial])
    eigenvalues = []
    for factor, multiplicity in exact.factor()[1]:
        for root, repeat in factor.complex_roots():
            assert root.imag.contains(0) and -root.real > 8
            eigenvalues.extend([-root.real] * (multiplicity*repeat))
    assert len(eigenvalues) == 18

    shell = neutral["shell"]
    soft = shell["soft"]
    pole, orbit = soft["pole"], soft["star"]["orbit"]
    spokes = [i+j-pole for i, j, _, _ in shell["geo"]["edges"] if pole in (i, j)]
    softened = [other for other in spokes if other in orbit]
    assert len(softened) == len(set(softened)) == 3
    assert soft["point"]["study"]["offset"] > 0 and acb.zeta_zero(8359).imag > 8462
    delta = arb(23)/25*arb(2)**-20/arb(8462)**2
    assert min(eigenvalues) - 4*delta > 8
    error = 18*4*delta/(6*arb(3).sqrt()*8)

    def unit_sensitivity(y):
        return sum(lam*y*(1-y)/(lam-(lam-1)*y)**2 for lam in eigenvalues)

    values, baselines = [], []
    for nodes, probabilities in rules:
        endpoint, *interior = nodes
        brackets = []
        for lam in eigenvalues:
            r = lam/(lam-1)
            # 取自由节点 y=0，给差商符号所用括号的上界。
            bracket = denominator.derivative()(r)/denominator(r) + 1/(r-1)
            bracket -= 1/(r-endpoint) + 2*sum(1/(r-v) for v in interior)
            assert bracket < -1
            brackets.append(bracket)
        values.append(sum(p*unit_sensitivity(y) for y, p in zip(nodes, probabilities)))
        baselines.append(sum(p*(1+sum(y/(lam-(lam-1)*y) for lam in eigenvalues))
                             for y, p in zip(nodes, probabilities)))
        print("固定端点、差商括号上界", endpoint.str(2), max(z.upper() for z in brackets).str(14))

    lower, upper = values[1]-error, values[0]+error
    assert arb(".00208104145403") < lower < upper < arb(".00208104148066")
    print("固定算术权的最大一阶绝对变化", lower.str(23), upper.str(23))

    assert (arb(SCALE)/HEIGHTS[0])**2 < min(eigenvalues)
    for index, (mass, h, limit) in enumerate(zip(band_masses, (HEIGHTS[1], HEIGHTS[3], HEIGHTS[4]),
                                                 (".606546", ".196015", ".00205319"))):
        x = (arb(SCALE)/h)**2
        value = unit_sensitivity(x/(1+x))
        if index < 2:
            share = mass*(value-error)/upper
            assert share > arb(limit)
            print("刚度响应份额下界", limit)
        else:
            share = mass*(value+error)/lower
            assert share < arb(limit)
            print("刚度响应份额上界", limit)

    baseline_low = arb("1.0020816770805")
    baseline_high = arb("1.0020816771057")
    assert baseline_low < baselines[0] < baselines[1]+18*delta < baseline_high
    eps = arb("1e-7")
    drop = eps/(1+eps)*lower-(baseline_high-baseline_low)
    assert drop > arb("1.82e-10")
    print("同向加硬后任意重选算术权的平均响应降幅下界", drop.str(18))


def verify():
    extrapolation = run_path(str(ROOT / "ch39-six-time-extrapolation.py"))
    bands = run_path(str(ROOT / "ch39-six-time-height-bands.py"))
    radial = bands["radial"]
    ctx.prec, ctx.cap = 512, 3
    neutral = radial["neutral"]
    # times() 独立求值归一 ξ 对数导数，并核对归档的 288 位小数格。
    readings = neutral["model"]["shared"]["clock"]["clock"]["times"]()[:8]
    reference = readings.pop()
    readings[6] += arb(SIXTH_SHIFT)
    denominator, moments = radial["compact_moments"](readings)
    rules = [extrapolation["quadrature"](end, moments, denominator) for end in (0, 1)]
    held_out(extrapolation, reference, moments, rules)
    robust_bands()

    cutoff, nodes, weights = height_rule(200_000_000, bands, moments, denominator)
    low_mass = sum(w for y, w in zip(nodes, weights) if y > cutoff)
    assert low_mass > arb("2.44066534e-11")
    print("高度严格低于 2e8 的算术概率下界", "2.44066534e-11")

    tails = []
    for height in HEIGHTS:
        cutoff, nodes, weights = height_rule(height, bands, moments, denominator)
        below = sum(w for y, w in zip(nodes, weights) if y < cutoff)
        tails.append(below.union(below+weights[0]))
    band_masses = (tails[0]-tails[1], tails[2]-tails[3], tails[4])
    stiffness(neutral, denominator, rules, band_masses)


if __name__ == "__main__":
    verify()
