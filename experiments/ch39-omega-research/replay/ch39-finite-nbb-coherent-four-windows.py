"""Bound one real two-clock filter against the four-atom RH calibration."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, acb_series, arb, arb_mat, arb_series, ctx

record = run_path(str(Path(__file__).with_name("ch39-finite-rh-pivot-mass.py")))
two = arb(2)


def readings(count):
    def xi(j):
        s = arb_series([arb(1)/2+j, 1])
        return s*(s-1)/2 * ((-s/2)*arb.pi().log()).exp() * (s/2).gamma()*s.zeta()
    anchor = xi(0)
    mass = anchor[2]/anchor[0]
    levels = [arb(1)]
    for j in range(1, count):
        value = xi(j)
        levels.append(value[1]/(2*j*value[0]*mass))
    return levels


def slow_clock():
    with ctx.workprec(2048):
        n, e = (arb(3).log()/two.log()).mid().man_exp()
        ratio = (Fraction(int(n))*Fraction(2)**int(e)).limit_denominator(2**860)
        p, q = ratio.numerator, ratio.denominator
        beat = abs(p*two.log()-q*arb(3).log())
        assert two**-865 < beat and beat < two**-864
        return p, q, beat


def bound_filter(heights, weights):
    delta, guard = two**-180, two**-160
    B, N = two**850, 2**1800
    assert all(h-g > 1 for g,h in zip(heights,heights[1:]))
    p, q, beat = slow_clock()
    assert beat*(B+100)<1 and max(2*p*N+8, 2*q*N)<2**2663
    # On the circle the even step is one gain per signed band, zero elsewhere.
    # At distance guard*delta from every edge, convolving with normalized Fejer
    # costs at most (oscillation of the step) * pi/(N*beat*guard*delta).
    gain_limit = 128
    spill = gain_limit*arb.pi()/(N*beat*delta*guard)+2*N*two**-2400
    assert spill < two**-580
    gains = []
    for g, w in zip(heights, weights):
        derivative = abs(acb_series([acb(arb(1)/2, g), 1]).zeta()[1])
        region = g+arb(0, 2*delta)
        second = abs(2*acb_series([acb(arb(1)/2, region), 1]).zeta()[2])
        assert derivative > arb(".75") and second < 8
        pole = lambda x: 1-arb(8)/17*(arb(16).log()*x).cos()
        gain = (w*(g*g+arb(1)/4)).sqrt()/(pole(g)*derivative)
        assert gain > 1 and gain < gain_limit-1 and region < 100
        gains.append(gain)
        # The squared density divided by w*u^2/(2*pi) when the gain
        # equals its ideal value. Bound the Taylor quotient also at u=0.
        zeta_ratio = 1+arb(0, 8*delta/derivative)
        density_ratio = (pole(region)/pole(g)*zeta_ratio)**2 * (g*g+arb(1)/4)/(region**2+arb(1)/4)
        assert abs(density_ratio-1) < two**-170
    # Each band pair has ideal mass w_i*S. Core lower bounds have their own
    # gain. Up to each outer edge the other arcs contribute at most spill
    # (Fejer is positive and each other arc is separated by the guard).
    S = 2*delta**3/(3*arb.pi())
    inner = (1-guard)**3*(1-two**-170)*(1-spill)**2*(1-100*two**-200)**2
    outer = (1+guard)**3*(1+two**-170)*(1+spill)**2
    assert inner>1-two**-157 and outer<1+two**-157
    # Raw spectral mass is <=2. The high-frequency tail estimate follows
    # from the real fractional-part dilation difference, averaged at 0..2/B.
    tail = 80*(1+B.log())/B
    s, T = 2/B, B/4
    difference2 = 2*(2*s*s*T+2*s*(1+T.log())+1/T)+4*s*s
    assert difference2<tail
    outside = (arb(25)/17)**2*(2*spill**2+(gain_limit+1)**2*tail)/S
    assert outside<two**-270
    # For bounded even Lipschitz tests (absolute value and Lipschitz
    # constant at most 1), normalize the four masses plus the exterior.
    mass_error = two**-157 + outside
    expectation_error = (2*mass_error+(1+guard)*delta*outer)/(1-mass_error)
    assert expectation_error < two**-155
    return expectation_error, gains


def main():
    ctx.prec, ctx.cap = 480, 3
    levels = readings(5)
    heights = [acb.zeta_zero(j).imag for j in (1, 2, 13, 15)]
    matrix = arb_mat([[g*g/(g*g+j*j) for g in heights] for j in range(4)])
    weights = matrix.solve(arb_mat([[f] for f in levels[:4]])).entries()
    assert min(weights) > arb(".1665") and sum(weights).contains(1)
    _, gains = bound_filter(heights, weights)
    # f_j(t)=t^2/(t^2+j^2), j>=1, has |f'_j|<=1.
    model = [arb(1)]+[levels[j]+arb(0,two**-155) for j in (1,2,3)]
    for j in range(1,4):
        for i in range(j):
            h = (i*model[i]+j*model[j])/((i+j)*(model[i]*model[j]).sqrt())
            assert (h*record["den"]).floor().unique_fmpz()==record["bins"][j][i]
    fifth = sum(w*g*g/(g*g+16) for w,g in zip(weights,heights))
    gap = levels[4]-fifth+arb(0,two**-155)
    assert not gap.contains(0)
    print("Four coherent paired windows. Filter order 2^1800; each prime exponent < 2^2663.")
    print("Absolute error < 2^-155 for the bounded even 1-Lipschitz readings.")
    print("All six Gram bins reproduced. Full RH minus filtered fifth squared overlap",gap.str(16))
    print("Positive passband gains", [a.str(12) for a in gains])


if __name__ == "__main__":
    main()
