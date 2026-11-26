"""在同一组等权路径数下，从原生壳边重放重子与味的两条闭路。"""
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
field = run_path(str(root / "ch39-five-charge-rational-populations.py"))
spectrum = run_path(str(root / "ch39-port-coupling-competition.py"))
sites, pole, radius, dot = (field[k] for k in ("sites", "pole", "radius", "dot"))
a, b = 71672243, 294016886
baryon_loops = (
    ((0, 10, 14), 28), ((0, 14, 13), -12), ((6, 8, 7), -12),
    ((6, 13, 7), -2), ((6, 14, 7), -2), ((6, 9, 13), 1),
    ((6, 9, 14), 8), ((6, 11, 15), 7), ((6, 15, 14), 7),
    ((9, 13, 14), -5), ((9, 14, 18), 9),
)
flavor_loops = (((6, 12, 15), 1), ((6, 12, 17), 1), ((0, 10, 16), 5))


def primitive(x, y):
    divisor = gcd(x, y)
    return x // divisor, y // divisor


def multiply(z, w):
    # Keep the oriented ray in Z[omega], where omega^2 = -1-omega.
    x, y = z
    u, v = w
    return primitive(x*u-y*v, x*v+y*u-y*v)


def native_edge(i, j):
    x, y = sites[i], sites[j]
    distance = 2*(radius-dot(x, y))
    assert 0 < distance <= 804 and distance not in (724, 804)
    A = radius*(radius+dot(pole, x)+dot(x, y)+dot(y, pole))
    B = -17*sum(pole[k]*(x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3])
                for k in range(3))
    sign = -1 if distance == 200 else 1
    return primitive(sign*(A+B), sign*2*B)


def word(loops):
    result = (1, 0)
    for (i, j, k), count in loops:
        # Bring every triangle to the same basepoint (the pole).
        if count < 0:
            j, k = k, j
        walk = [i, j, k, i]
        start = sites.index(pole)
        if i != start:
            walk = [start, *walk, start]
        triangle = (1, 0)
        for source, target in zip(walk, walk[1:]):
            triangle = multiply(triangle, native_edge(source, target))
        for _ in range(abs(count)):
            result = multiply(result, triangle)
    return result


def verify():
    ctx.prec = 192
    field["verify"]()
    unit = word(baryon_loops)
    assert unit == (-35, 13)
    power = (1, 0)
    for _ in range(5):
        power = multiply(power, unit)
    x, y = power
    norm = x*x-x*y+y*y
    length = isqrt(norm)
    assert length**2 == norm
    assert Fraction(2*x-y, 2*length) == -Fraction(a, b)
    assert b % 2 == 0 and a+b == 365689129
    # Current entries in Z[omega], scaled by (a+b)*length.
    current = [(0, 0)]*len(sites)
    rotate = [sites.index((v[1], v[2], v[0])) for v in sites]
    assert rotate[sites.index(pole)] == sites.index(pole)
    orbit = (0, rotate[0], rotate[rotate[0]])
    assert len(set(orbit)) == 3
    current[orbit[0]] = (a*length, 0)
    current[orbit[1]] = ((b//2)*x, (b//2)*y)
    current[orbit[2]] = ((b//2)*(x-y), -(b//2)*y)
    for i in range(len(sites)):
        rows = [current[j] for j in (i, rotate[i], rotate[rotate[i]])]
        assert tuple(map(sum, zip(*rows))) == (0, 0)
    print("Baryon oriented ray, real part at charge five, equal reverse count",
          unit, -Fraction(a, b), b//2, flush=True)

    x, y = word(flavor_loops)
    phase = acb(2*x-y, arb(3).sqrt()*y)
    phase /= abs(phase)
    assert arb(".89752504451") < phase.arg() < arb(".89752504452")
    print("Flavor phase", phase.arg().str(18), flush=True)
    vector = acb_mat([[(arb(a)/(a+b)).sqrt()], [phase*(arb(b)/(a+b)).sqrt()]])
    height = acb.zeta_zero(11431).imag
    gap = spectrum["full_spectrum"](height, vector, arb(2)**-10,
                                    curvature=arb(1)/10000, allowance=arb("7.8e-9"))
    assert gap > arb("3e-16")


if __name__ == "__main__":
    verify()
