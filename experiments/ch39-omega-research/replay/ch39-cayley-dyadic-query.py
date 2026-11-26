"""Certify the Cayley oracle cost of the dyadic phase on the six candidate zeros."""
from itertools import combinations

from flint import acb, arb, ctx


def verify():
    ctx.prec = 168
    heights = [z.imag for z in acb.zeta_zeros(11485, 6)]
    half = arb(1)/2
    quarter = half**2
    log_two = arb(2).log()
    tolerance = arb(1)/8192
    lower = []
    for i, j in combinations(range(6), 2):
        g, h = heights[i], heights[j]
        target_gap = 2*abs(((h-g)*log_two/2).sin())
        input_gap = (h-g)/((g*g+quarter)*(h*h+quarter)).sqrt()
        lower.append(((target_gap-2*tolerance)/input_gap, i, j))
    winner, i, j = max(lower, key=lambda v: float(v[0].mid()))
    assert (i, j) == (4, 5)
    assert all(v < winner for v, i, j in lower if (i,j) != (4,5))
    assert winner.ceil().unique_fmpz() == 85562026

    midpoint = (heights[0]+heights[-1])/2
    degree = int((log_two*(midpoint**2+quarter)+half).floor().unique_fmpz())
    assert degree == 85970718
    def phase(g):
        return g*log_two + 2*degree*(half/g).atan()
    critical = (degree/log_two-quarter).sqrt()
    assert heights[0] < critical < heights[-1]
    low = phase(critical)
    left, right = phase(heights[0]), phase(heights[-1])
    assert left > right
    center = (low+left)/2
    slope = arb(178262298)/2**28
    fixed = acb(1, slope) / acb(1, -slope)
    rounding = abs(acb(0, center).exp()-fixed)
    uniform_error = (left-low)/2 + rounding
    assert uniform_error < arb('0.00010592') < tolerance
    for g in heights:
        input_phase = acb(g, half)/acb(g, -half)
        assert abs(fixed*input_phase**(-degree)-acb(0, g*log_two).exp()) < uniform_error

    k = arb(4)/9
    h = arb(2).sqrt()
    a, b = 1+2*k*k, 2*h*k
    def envelope(v):
        return (3*h-4*v)*(a+b*v).sqrt()
    peak = (3*h*b-8*a)/(12*b)
    assert -1 < peak < 1 and envelope(1) < envelope(-1)
    radius = (envelope(peak)-envelope(1))/2
    assert arb('2.56221809') < radius < arb('2.56221810')
    print('Cayley calls at tolerance 1/8192:', int(winner.ceil().unique_fmpz()), degree)
    print('Uniform error upper bound', uniform_error.upper().str(14))
    print('Full-spectrum uniform envelope error floor', radius.str(14))


if __name__ == '__main__':
    verify()
