"""相对逆高度平方先验的熵极小谱权及全部三窗候选上的代价界。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

root = Path(__file__).parent
ranges = run_path(str(root / "ch39-eleven-neutral-cp-range.py"))
joint = ranges["joint"]


def column(values):
    return arb_mat([[v] for v in values])


def projection_weights(rows, times):
    basis = arb_mat(list(zip(*(r[0] for r in rows[:8]), strict=True)))
    extra = arb_mat(list(zip(*(r[0] for r in rows[8:]), strict=True)))
    start = column(basis.solve(column([*times[:6], 0, 0])).entries()+[arb(0)]*3)
    front = -basis.solve(extra)
    tangent = arb_mat([[front[i,j] if i < 8 else int(i == 8+j)
                        for j in range(3)] for i in range(11)])
    xs = [r[1] for r in rows]

    def derivatives(t):
        p = (start+tangent*t).entries()
        assert all(v > 0 for v in p)
        g = tangent.transpose()*column([(v/x).log() for v, x in zip(p, xs, strict=True)])
        h = arb_mat([[sum(tangent[k,i]*tangent[k,j]/p[k] for k in range(11))
                      for j in range(3)] for i in range(3)])
        return g, h

    center = column([arb(s) for s in (".06458424719409188", ".16361116312614066", ".04697193110629481")])
    radius = arb("1e-16")
    box = center+column([arb(0,radius)]*3)
    gradient, hessian = derivatives(center)
    inverse = hessian.inv()
    _, full_hessian = derivatives(box)
    identity = arb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    remainder = identity-inverse*full_hessian
    norm = max(sum(abs(remainder[i,j]) for j in range(3)) for i in range(3))
    step = max(abs(v) for v in (inverse*gradient).entries())
    assert step+norm*radius < radius
    return (start+tangent*box).entries()


def finite_projection(rows, times, probabilities, reference):
    p = projection_weights(rows, times)
    xs = [r[1] for r in rows]
    loss = ranges["cp_loss"](rows, p, reference)
    flavor = [v*r[2] for v, r in zip(p, rows, strict=True)]
    charm = sum(v*(r[3][7]-probabilities[7]) for v, r in zip(flavor, rows, strict=True))/sum(flavor)
    divergence = sum(v*(v*sum(xs)/x).log() for v, x in zip(p, xs, strict=True))
    assert arb(".068381186952") < loss < arb(".068381186954")
    assert arb("2.21363e-12") < charm < arb("2.21365e-12")
    assert arb("34.7133097987") < divergence < arb("34.7133097988")
    print("十一枚熵极小解的 CP 亏额、第二概率超额与相对熵",
          loss.str(15), charm.str(12), divergence.str(15))
    return divergence, sum(xs)/2_903_600_000**2


def inverse_square_tail(h):
    h = arb(h)
    return ((h/(2*arb.pi())).log()+1)/(2*arb.pi()*h)+arb(".28")*(2*h.log()+arb(".5"))/h**2


def spectrum_bounds(times, finite_divergence, finite_mass):
    short = run_path(str(root / "ch39-short-arc-born-weights.py"))
    ctx.prec = 176
    prior_mass = sum(1/acb.zeta_zero(n).imag**2 for n in short["indices"])
    assert prior_mass > arb("2.95816e-8")
    heights = [30_000_000_000, 5*10**12, 15*10**12]
    x = (arb("2899267716.58")/heights[0])**2
    # Π (j²x)/(1+j²x) 随 x 增加，期望由五个已知时间经部分分式给出。
    ids = range(1,6)
    factors = [arb(j*j)*x for j in ids]
    at_cut = prod(v/(1+v) for v in factors)
    expected = arb(1)
    for j in ids:
        pole = -1/arb(j*j)
        coefficient = prod(k*k*pole for k in ids)/prod(1+k*k*pole for k in ids if k != j)
        expected += coefficient*times[j]
    mass = 1-expected/at_cut
    assert mass > arb(".988029")
    # 另两条尾质量下界由 ch39-neutral-five-high-ordinates.py 在同一尺度范围给出。
    bounds = [mass, arb(".37426"), arb(".23893")]
    levels = [arb(i) for i in (0,9,13,14)]
    partition = arb(1)
    mean = arb(0)
    for h, bound, low, high in zip(heights, bounds, levels[:-1], levels[1:], strict=True):
        partition += (high.exp()-low.exp())*inverse_square_tail(h)/prior_mass
        mean += (high-low)*bound
    lower = mean-partition.log()
    upper = finite_divergence+(inverse_square_tail(8000)/finite_mass).log()
    assert acb.zeta_zero(8359).imag > 8000
    assert lower > arb("6.60014") and upper < arb("44.07337")
    print("三百亿高度以上的算术概率下界与全谱熵代价界",
          mass.str(12), lower.str(13), upper.str(13))


def verify():
    rows, times, _, _, probabilities, reference = joint["columns_at"](2_903_600_000, ranges["indices"])
    ctx.prec = 640
    value, mass = finite_projection(rows, times, probabilities, reference)
    spectrum_bounds(times, value, mass)


if __name__ == "__main__":
    verify()
