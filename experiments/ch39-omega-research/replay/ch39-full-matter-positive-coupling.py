"""任意完全激活的五、十、十五维壳物质不容许正的共同逆耦合。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_poly, ctx

fifteen = run_path(str(Path(__file__).with_name("ch39-symmetric-fifteen-threshold.py")))
matter = fifteen["matter"]
soft = matter["soft"]


def fifteen_cost():
    h = acb.zeta_zero(11431).imag
    base = soft["vector_pencil"](h)[0]
    blocks = {q:matter["nonspokes"](q, h) for q in (1, 4, 6)}
    tables = {q:matter["coefficients"](base, b) for q, b in blocks.items()}
    maximum = fifteen["derivative_bounds"](tables)
    # A test vector with no support on the pole: its star Rayleigh quotient is one.
    v = acb_mat([[int(i == 7)] for i in range(19)])
    for _ in range(10):
        v = blocks[1]*v/24
    adj = soft["star"]["point"]["adjoint"]
    assert v[soft["pole"], 0] == 0
    assert (adj(v)*blocks[1]*v)[0, 0].real > 21*(adj(v)*v)[0, 0].real
    radius = (matter["geo"]["theta"](804)/matter["geo"]["theta"](12)
              *arb(".999").union(arb("1.001")))
    assert radius**2 > arb("72.69") and radius**2 < 74

    def polynomials(x):
        return {q:arb_poly([sum(value*x**j for (l, j), value in t.items() if l == i)
                            for i in range(19)]) for q, t in tables.items()}

    def at_mass(x, upper, *, increasing):
        data = polynomials(x)
        a, b, c = (data[q] for q in (1, 4, 6))
        if increasing:
            deriv = 11*b.derivative()*a*c-4*a.derivative()*b*c-7*c.derivative()*a*b
            # The degree-53 term is exactly zero from 11-4-7. On [0,upper]
            # substitute e=upper*t/(1+t) and clear the positive denominator.
            assert deriv[53].contains(0)
            for k in range(53):
                assert sum(deriv[i]*upper**i*comb(52-i, k-i) for i in range(k+1)) > 0
        return (11*b(upper).log()-4*a(upper).log()-7*c(upper).log())/(5*arb.pi())

    def cost(mean):
        return 7/(3*arb.pi())*arb(19)/2*(arb("72.69")/mean).log()

    # If e <= 1, use mass monotonicity and the global massless bound.
    assert cost(5+arb(270)/19) > 16*maximum
    high_mass = at_mass(5, arb(1), increasing=True)
    assert cost(74-19+arb(270)/19) > 16*high_mass

    # If e >= 1, full activation gives x+1+21e <= r^2 and x+19 <= r^2.
    unit = at_mass(0, arb(1), increasing=False)
    for lo, hi in ((0, 5), (5, 15), (15, 30), (30, 55)):
        bound = unit if lo == 0 else at_mass(lo, arb(73-lo)/21, increasing=True)
        mean = hi+(36+234*arb(73-hi)/21)/19
        assert cost(mean) > 16*bound
    return radius


def verify():
    ctx.prec = 224
    radius = fifteen_cost()
    # S = (-1,1,1) dot the three inverses. S+16C equals the common inverse at a crossing.
    coeff = (-1+16, 1-48, 1+32)
    assert sum(coeff) == 1 and sum(v*w for v, w in zip(coeff, (5,3,2))) == 0
    mz = arb("91.1876")
    em, w = arb("127.930"), arb(".23129")
    bare = (3*(1-w)*em/5, w*em, 1/arb(".1180"))
    beta = (arb(41)/10, -arb(19)/6, -7)
    rate = -sum(x*y for x, y in zip(coeff, beta))/(2*arb.pi())
    assert rate > 0
    light = sum(x*y for x, y in zip(coeff, bare))+rate*(arb("1e17")/mz).log()
    neutral = 7*5/(4*arb.pi())*(18*radius.log()-arb(19).log()/2)
    print("Light and vector upper, common inverse at a possible equality <=",
          light.upper().str(13), neutral.upper().str(13), (light+neutral).upper().str(13))
    assert light+neutral < -8


if __name__ == "__main__":
    verify()
