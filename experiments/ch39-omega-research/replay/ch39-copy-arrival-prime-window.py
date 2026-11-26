"""Recover the four-tap prime filter from copy-clock arrivals; compare raw curves."""
from fractions import Fraction
from math import gcd
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
arithmetic = run_path(str(root/'ch39-dyadic-prime-window-pivot.py'))


def convolve(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def check_arrivals():
    active = Fraction(5,6)
    for retained in (Fraction(1,2), Fraction(1,3), Fraction(3,4)):
        times = [(1/retained-1)/f for f in (4,8)]
        assert times[0] == 2*times[1]
        for f, t in zip((4,8), times, strict=True):
            assert (1-active) + active/(1+f*t) == (1-active)+active*retained
    taps = convolve(convolve([-1,1], [-2,1]), [432,192])
    taps = [v//gcd(*taps) for v in taps]
    assert taps == [18,-19,-3,4]
    assert sum(taps) == sum(v*2**j for j,v in enumerate(taps)) == 0
    for x in (Fraction(-4), Fraction(3,4), Fraction(10,3), Fraction(-9,4)):
        assert sum(v*x**j for j,v in enumerate(taps)) == (x-1)*(x-2)*(9+4*x)


def check_raw_curves():
    ctx.prec, ctx.cap = 240, 3
    heights = [z.imag for z in acb.zeta_zeros(1,8)]
    k = arb(4)/9
    h = arb(2).sqrt()
    log_two = arb(2).log()
    bare_total = 1+arb.const_euler()/2-(4*arb.pi()).log()/2
    bare_tail = bare_total-sum(1/(g*g+arb(1)/4) for g in heights[:-1])
    assert bare_tail > 0
    def window(g):
        v = (log_two*g).cos()
        return (3*h-4*v)*(1+2*k*k+2*h*k*v).sqrt()
    for power in (1,2):
        def common(g):
            return ((arb.pi()*heights[0]).cosh()/(arb.pi()*g).cosh())**(arb(power)/2)
        masses = [window(g)*common(g)/(g*g+arb(1)/4) for g in heights[:-1]]
        tail = (3*h+4)*(1+h*k)*common(heights[-1])*bare_tail
        def integrate(f):
            return (sum(r*f(g) for r,g in zip(masses, heights[:-1], strict=True))
                    +arb(0).union(tail*f(heights[-1])))
        a = integrate(lambda g: arb(1))
        b = integrate(lambda g: 1/(g*g+1))
        c = integrate(lambda g: 4/((g*g+1)*(g*g+4)))
        q = b/a
        pivot = arithmetic['schur_readout'](q,c/b)
        if power == 1:
            assert arb('.99998674') < masses[0]/a < arb('.99998675')
            assert 50846 < pivot < 50847
        else:
            assert arb('.99999999973') < masses[0]/a < arb('.99999999974')
            assert 11372954 < pivot < 11372955
        assert arb('.99501968') < 1-q < arb('.99501973')
        print('Kernel power', power, 'first-zero share', (masses[0]/a).str(14),
              'T1', (1-q).str(14), 'pivot', pivot.str(13))
    assert arb('.99839735') < arithmetic['arithmetic_times']()[0] < arb('.99839736')


if __name__ == '__main__':
    check_arrivals()
    check_raw_curves()
