"""Count compatible labels using monotonicity and the exact corner shifts."""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import fraction


def verify():
    ctx.prec = 180
    source = Path(__file__).parent
    graph = run_path(str(source/"ch39-weighted-arc-zero-phase.py"))
    a, initial, fit = run_path(str(source/"ch39-shell-zero-mass-hierarchy.py"))["parameters"]()
    num, den = fraction(graph["response"](3))
    evaluate = graph["evaluated"]
    target = a("1.2209e19")
    pi = arb.pi()

    def epsilon(height):
        u = 1-1/(2*(height**2+a("1/4")))
        return evaluate(num,u)/evaluate(den,u)

    @cache
    def deficit(j):
        h = acb.zeta_zero(j).imag
        low, v, _, _, _ = fit(epsilon(h), a(1))
        mass = v/low*((2*pi/(3*arb(5).sqrt()/2))*(h**2+a("1/4")+30+a("2/15"))).sqrt()
        return (mass/target).log()

    coefficients = (-220, 648, -428)
    def corner(e, direction, all_three):
        total = arb(0)
        for i, (inverse, c) in enumerate(zip(initial, coefficients)):
            if all_three or i == 1:
                sign = 1 if c > 0 else -1
                total += c*inverse*((1-direction*sign*e)**(-2)-1)
        return total*pi/21

    # Bounds of the possible index interval for one fluctuating input and for all.
    windows = [("1e-6", (6954, 7051), (6913, 7094)),
               ("1e-7", (6998, 7007), (6994, 7011)),
               ("1e-8", (7003, 7003), (7002, 7004)),
               ("1e-9", (7004, 7003), (7004, 7003))]
    for radius, weak, all_couplings in windows:
        for all_three, (first, last) in ((False, weak), (True, all_couplings)):
            lower = corner(a(radius), -1, all_three)
            upper = corner(a(radius), 1, all_three)
            assert deficit(first-1)+upper < 0 < deficit(first)+upper
            assert deficit(last)+lower < 0 < deficit(last+1)+lower
            print(radius, "three" if all_three else "weak", max(0, last-first+1))

    # Enclose masses and perturbative conditions throughout the widest input box.
    width = a("1e-6")
    inv_box = [(x/(1+width)**2).union(x/(1-width)**2) for x in initial]
    ebox = epsilon(acb.zeta_zero(6913).imag).union(epsilon(acb.zeta_zero(7094).imag))
    low, v, log_x, g5, g10 = fit(ebox, a(1), inv_box)
    assert g5**2 < pi and g10**2 < pi and log_x.exp()**2 < 4*pi
    assert ebox < 1 and v*log_x.exp()*ebox > 1000


if __name__ == "__main__":
    verify()
