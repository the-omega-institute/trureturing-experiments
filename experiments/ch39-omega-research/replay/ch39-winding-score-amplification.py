"""Uniform coherent amplification of the bare Cayley factor on the certified matching band."""
from fractions import Fraction
from math import factorial, prod

from flint import acb, arb, ctx


def ball(x):
    q = Fraction(x)
    return arb(q.numerator)/q.denominator


def interpolation(m, order, scale):
    nu = Fraction(1, m)
    nodes = tuple(range(-order, order+1))
    weights = [prod((nu-i)/(j-i) for i in nodes if i!=j) for j in nodes]
    anchor = acb(22000,1)/acb(22000,-1)
    # Combining the two identity coefficients before taking the l1 norm is essential.
    l1 = abs(1-anchor*ball(weights[order])) + sum(ball(abs(w)) for j,w in zip(nodes,weights) if j)
    assert scale*l1 < 1
    theta_anchor = 2*(arb(1)/22000).atan()
    offsets = [abs(2*(arb(1)/(2*g)).atan()-theta_anchor) for g in (10570,11644)]
    angle = m * max(offsets)
    error = ball(abs(prod(nu-i for i in nodes))) * angle**len(nodes)/factorial(len(nodes))
    relative = error * (arb(11644)**2+arb(1)/4).sqrt()
    return nodes, weights, anchor, relative


def verify():
    ctx.prec = 160
    heights = [z.imag for z in acb.zeta_zeros(10855,1238)]
    for m,order,scale,limit,minimum in [(16384,1,6500,".0000663",".065"),
                                        (131072,3,9500,".0000343",".29")]:
        nodes, weights, anchor, relative = interpolation(m, order, scale)
        assert relative < arb(limit)
        intervals = []
        for g in heights:
            d = arb(2).sqrt()*acb(0,g*arb(2).log()).exp()
            raw_modulus = abs((d-1)*(d-2)*(9+4*d))
            bare = 1/(g*g+arb(1)/4)
            ideal = scale**4 * bare**2 * raw_modulus**2 / 5256
            intervals.append((ideal*(1-relative)**4, ideal*(1+relative)**4))
        target = 10875-10855
        low, high = intervals[target]
        closest = max(top for j,(_,top) in enumerate(intervals) if j!=target)
        assert low > closest and low > arb(minimum)
        g = heights[target]
        w = acb(g,arb(1)/2)/acb(g,-arb(1)/2)
        phase = (w/anchor)**m
        amplitude = scale*(1-anchor*sum(ball(c)*phase**j for c,j in zip(weights,nodes)))
        d = arb(2).sqrt()*acb(0,g*arb(2).log()).exp()
        success = abs(amplitude)**4*abs((d-1)**2*(9*d+8))**2/5256
        assert low < success < high
        print("Queries",2*m*order,"winner probability",success.str(13),
              "relative amplitude error <",relative.upper().str(9),
              "certified probability gap >", (low-closest).lower().str(9))


if __name__ == "__main__":
    verify()
