"""有限计数链的详细平衡、实际零点初始几何、复相位与平直变分。"""
from fractions import Fraction
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def rational_chain():
    n, x = 9, Fraction(1,5)
    raw = [(n*n*x)**k/factorial(k)**2 for k in range(n+1)]
    total = sum(raw)
    pi = [w/total for w in raw]
    next_pi = [Fraction(0)] * (n+1)
    for k in range(n+1):
        up = x/2 if k<n else 0
        down = Fraction(k*k,2*n*n)
        next_pi[k] += pi[k]*(1-up-down)
        if k < n:
            next_pi[k+1] += pi[k]*up
            assert pi[k]*up == pi[k+1]*Fraction((k+1)**2,2*n*n)
        if k:
            next_pi[k-1] += pi[k]*down
    assert next_pi == pi
    assert sum(pi[k]*Fraction(k*k,n*n) for k in range(n+1)) == x*(1-pi[n])


def moments_phase(n, area, theta):
    x = area**2
    weight = arb(1)
    total, centered, square = arb(0), arb(0), arb(0)
    numerator = acb(0)
    for k in range(n+1):
        delta = arb(k)/n-area
        total += weight
        centered += weight*delta
        square += weight*delta**2
        numerator += weight*acb(0,k*theta).exp()
        ratio = n*n*x/(k+1)**2
        weight *= ratio
    shift = centered/total
    noise = n*(square/total-shift**2)
    return area+shift, noise/x, noise, numerator/total


def verify():
    ctx.prec = 640
    rational_chain()
    q = 2*arb(3).sqrt()-3
    b = 4*q/(arb(5)/3).sqrt()
    expected = -54*b/(arb(5).sqrt()*8*arb.pi())
    for n in (64, 512, 5432):
        shape_force, real_part, imaginary_part = arb(0), arb(0), arb(0)
        for count, area, dx, angle1, angle2 in [
            (30, 1-q, b*b/6, arb(18)/5, -arb(108)/25),
            (20, q, b*b/4, -18*arb(3).sqrt()/5, 162*arb(3).sqrt()/25)
        ]:
            mean, slope, noise, factor = moments_phase(n,area,arb(1)/2)
            shape_force += count*mean*angle1/(8*arb.pi())
            real_part += count*(mean*angle2+2*slope*dx*angle1)/(8*arb.pi())
            imaginary_part += count*noise*(angle1/(8*arb.pi()))**2
            z = acb(0,arb(1)/4).exp()
            ratio = factor/((2*n*area*(z-1)).exp()/acb(0,arb(1)/8).exp())
            assert abs(ratio-1) < arb(".1")/n
            assert abs(mean-area+arb(1)/(4*n)) < arb(".1")/n**2
            assert abs(noise-area/2) < arb(".1")/n**2
        print("槽位、平直形状力、实 Hessian 与 Regge 之比、虚 Hessian",
              n,shape_force.str(12),(real_part/expected).str(12),imaginary_part.str(12))

    # 配到由实际 11366 号零点触发的 5432 槽几何，比较固定非零角亏的相位。
    model = run_path(str(Path(__file__).with_name("ch39-clock-hinge-ground.py")))
    n = 5432
    _, four, six, _ = model["matched_geometry"](n,2911)
    total = acb(0)
    for count, k, deficit in [(30,2911,four),(20,2521,six)]:
        area, theta = arb(k)/n, deficit/(8*arb.pi())
        _, _, _, factor = moments_phase(n,area,theta)
        assert factor.real > 0
        total += count*factor.log()/acb(0,n)
    print("原整数驻点上每胞归一复作用的实部与虚部",total.real.str(12),total.imag.str(12))


if __name__ == "__main__":
    verify()
