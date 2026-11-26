"""同一非辐边刚度下，三类物质的全部平滑一圈阈值无法产生正的共同值。"""
from functools import cache
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx, fmpz_mat

root = Path(__file__).parent
smooth = run_path(str(root/"ch39-smooth-two-stiffness-crossing.py"))
fifteen = run_path(str(root/"ch39-symmetric-fifteen-threshold.py"))
matter, soft = smooth["matter"], smooth["soft"]


def convolve(first, second):
    result = [arb_poly() for _ in range(len(first)+len(second)-1)]
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i+j] += a*b
    return result


def derivative_numerator(tables, K):
    polys = [[arb_poly([t.get((i, j), 0) if (i, j) != (0, 0) else 0 for j in range(20-i)])
              for i in range(19)] for t in tables]
    weights = (arb("2.4")*(1-K), arb(".9")+K*arb("6.6"), arb(".2")-K*arb("4.2"))
    total = [arb_poly() for _ in range(55)]
    for q, weight in enumerate(weights):
        product = convolve(convolve([p.derivative() for p in polys[q]], polys[(q+1)%3]),
                           polys[(q+2)%3])
        for i, value in enumerate(product):
            total[i] += weight*value
    return total


def positive_region(terms, x0, e0, *, width=None):
    shifted = [p(arb_poly([x0, 1])) for p in terms]
    for j in range(57):
        poly = arb_poly([shifted[i][j] for i in range(min(54, 56-j)+1)])(arb_poly([e0, 1]))
        coeffs = poly.coeffs()
        if width is not None:
            degree = poly.degree()
            coeffs = [sum(coeffs[i]*width**i*comb(degree-i, k-i) for i in range(k+1))
                      for k in range(degree+1)]
        assert all(c > 0 for c in coeffs)


def compact_region(base, blocks, radius):
    grid = 4096
    @cache
    def eigen(q, i):
        return smooth["masses"](base+(arb(i)/grid)*blocks[q])

    @cache
    def exposure(q, i, j, top):
        factor = arb("1.001") if top else arb(".999")
        return sum(smooth["sheet"](z+arb(j)/grid, (radius*factor)**2, False) for z in eigen(q,i))
    # Positive weights take the least exposure; negative weights the greatest.
    weights = ((1, arb("-16.8")), (4, arb("53.7")), (6, arb("-33.4")))
    todo = [(grid, 64*grid, 0, 5*grid)]
    visited = 0
    while todo:
        a, b, c, d = todo.pop()
        visited += 1
        lower = sum(w*(exposure(q,b,d,False) if w > 0 else exposure(q,a,c,True))
                    for q, w in weights)
        if lower > arb("1e-10"):
            continue
        if (b-a)*(grid//2+c) > (d-c)*(grid//2+a):
            mid = (a+b)//2
            assert a < mid < b
            todo.extend(((a,mid,c,d), (mid,b,c,d)))
        else:
            mid = (c+d)//2
            assert c < mid < d
            todo.extend(((a,b,c,mid), (a,b,mid,d)))
    # Each unshifted unit-stiffness mass in these three sectors exceeds one.
    for q in blocks:
        assert min(v.lower() for v in eigen(q,grid)) > 1
    ref2 = (radius*arb("1.001")*arb("91.1876")/arb("1e13"))**2
    assert sum(abs(w) for _,w in weights)*19*ref2/8 < arb("1e-10")
    print("Compact-region angular-integral boxes", visited)


def verify():
    ctx.prec = 192
    height = acb.zeta_zero(11431).imag
    base, neutral, _ = soft["vector_pencil"](height)
    blocks = {q:matter["nonspokes"](q,height) for q in (1,4,6)}
    tables = {q:matter["coefficients"](base,b) for q,b in blocks.items()}
    maximum = fifteen["derivative_bounds"](tables)
    K = arb(49)/4
    positive_region(derivative_numerator(tables.values(), K), 5, 0, width=1)
    larger = derivative_numerator(tables.values(), arb(8))
    positive_region(larger, 5, 1)
    positive_region(larger, 0, 64)
    geo = soft["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    rbox = radius*arb(".999").union(arb("1.001"))
    mean = 5+arb(270)/19
    ref2 = (rbox*arb("91.1876")/arb("1e13"))**2
    # log((y+z)/(y+z0)) is convex and decreasing in y: apply Jensen before subtracting.
    cost = 7*19/(6*arb.pi())*(arb(4)/5*(1+rbox**2/(4*mean)).log()-ref2/(5*mean))
    assert cost > K*maximum
    compact_region(base, blocks, radius)

    matrix = base+neutral
    integral = fmpz_mat([[int(v.real.unique_fmpz()) for v in row] for row in matrix.tolist()])
    unit = [z.real for z,k in integral.charpoly().complex_roots() if z != 0 for _ in range(k)]
    star = [arb(1)]*17+[arb(19)]
    em, w = arb("127.930"), arb(".23129")
    initial = (3*(1-w)*em/5, w*em, 1/arb(".1180"))
    L = (arb("1e17")/arb("91.1876")).log()
    for projection, spectrum, bound in ((K, star, arb("-1.3068")), (arb(8), unit, arb("-6.5068"))):
        coeff = (-1+projection, 1-3*projection, 1+2*projection)
        slope = -sum(v*b for v,b in zip(coeff, smooth["beta"]))/(2*arb.pi())
        assert slope > 0
        light = sum(a*x for a,x in zip(coeff, initial))+slope*L
        heavy = 35/(8*arb.pi())*sum((1+rbox**2/(4*x)).log() for x in spectrum)
        assert light+heavy < bound
        print("Separating projection coefficient and maximum common inverse",
              projection.str(6), (light+heavy).upper().str(14))


if __name__ == "__main__":
    verify()
