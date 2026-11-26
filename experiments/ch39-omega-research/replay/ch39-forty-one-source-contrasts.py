"""五至四十一项精确时间的读出倍率与中性壳响应最少因子数。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb_poly, arb, arb_mat, arb_poly, ctx


root = Path(__file__).parent
neutral = run_path(str(root / "ch39-neutral-clock-response.py"))
model = neutral["model"]
shared = model["shared"]


def quadrature(count, times):
    odd = count % 2
    degree = count//2
    moments = [v/200**j for j,v in enumerate(shared["clock"]["clock"]["cleared_moments"](times,count))]
    gram = arb_mat([[moments[odd+i+j] for j in range(degree)] for i in range(degree)])
    fit = gram.solve(arb_mat([[-moments[odd+degree+i]] for i in range(degree)])).entries()
    zeros = acb_poly([*fit,1]).roots(tol=arb("1e-600"))
    assert len(zeros)==degree and all(z.imag.contains(0) for z in zeros)
    roots = sorted(z.real for z in zeros)
    assert roots[0] > 0 and all(a<b for a,b in zip(roots,roots[1:]))
    nodes = ([arb(0)] if odd else []) + roots
    vandermonde = arb_mat([[x**i for x in nodes] for i in range(len(nodes))])
    tilted = vandermonde.solve(arb_mat([[v] for v in moments[:len(nodes)]])).entries()
    weights = [v*prod(1+j*j*x for j in range(1,count+1)) for x,v in zip(nodes,tilted,strict=True)]
    tail = prod(j*j for j in range(1,count+1))*(moments[count]-sum(v*x**count for v,x in zip(tilted,nodes,strict=True)))
    assert tail>0 and all(w>0 for w in weights) and (tail+sum(weights)).contains(1)
    for j in range(1,count+1):
        assert sum(w/(1+j*j*x) for x,w in zip(nodes,weights,strict=True)).overlaps(times[j])
    return nodes,weights,tail,([arb(0)] if odd else [])+roots+roots


def residues(count, pivots, poles, amplitudes):
    return [a*prod(1-z/(j*j) for j in range(1,count+1))/prod(1+z*x for x in pivots)
            for a,z in zip(amplitudes,poles,strict=True)]


def three_pole_sign(values):
    assert values[0]>0 and values[1]>0>values[2] and sum(values[:3])>0


def verify():
    ctx.prec, ctx.cap = 6144, 3
    times = shared["clock"]["clock"]["times"]()
    invariants,c0,poles,slopes = shared["flavor_data"]()
    assert 0<poles[0]<poles[1]<poles[2]<1
    amplitudes = [s/z for s,z in zip(slopes,poles,strict=True)]
    three_pole_sign(amplitudes)
    limit = sum(amplitudes)
    d1,d2,d3 = [p.derivative()(c0) for p in invariants]
    assert d1<0<d2 and d3>0
    crossing = (-d2+(d2*d2-4*d1*d3).sqrt())/(2*d3)
    # D'^2-DD'' 除以 x² 后，每一项的 x 系数在全物理相位上严格为正。
    es = [arb_poly([1]),*invariants]
    for k in range(2,7):
        coefficient = sum(a.derivative()*b.derivative()-a*b.derivative().derivative()
                          for i,a in enumerate(es) for j,b in enumerate(es) if i+j==k)
        assert coefficient(arb(0,1))>0
    poly = neutral["neutral_polynomial"]()
    for count in range(5,42):
        nodes,weights,tail,pivots = quadrature(count,times)
        assert max(nodes)<crossing
        scores = [model["gradient"](x,poles,slopes) for x in nodes]
        ratio = -sum(v*w for v,w in zip(scores,weights,strict=True))/(tail*limit)
        three_pole_sign(residues(count,pivots,poles,amplitudes))
        error = sum(pivots)+sum(1/arb(j*j) for j in range(1,count+1))
        surplus = crossing/(1+poles[2]*crossing)-error
        assert surplus>0
        factors = 1
        while 19**factors<ratio:
            factors += 1
        assert 19**(factors-1)<ratio<19**factors

        def endpoint_score(k):
            return tail*19**k*limit+sum(w*neutral["efficiency"](x,poly)**k*v
                                        for x,w,v in zip(nodes,weights,scores,strict=True))
        if count==15:
            assert factors==12
            # 对负 g，用谱均值 15 的 Jensen 界和一次 Bernoulli 界。
            upper = tail*19**12*limit+sum(w*(1+216*x/(15+x))*v
                                         for x,w,v in zip(nodes,weights,scores,strict=True))
            assert upper < arb("-1.27e-12")
            scaled = [a*(1+216/(1-15*z)) for a,z in zip(amplitudes,poles,strict=True)]
            extra = 216*model["gradient"](arb(-15),poles,slopes)
            corrected = residues(count,pivots,[*poles,arb(1)/15],[*scaled,extra])
            three_pole_sign(corrected)
            assert corrected[3]>0 and (19**12-217)*surplus>3240
            factors += 1
        assert endpoint_score(factors)>0
        print(count,"项的效率倍率下确界、中性响应最少因子数",ratio.str(13),factors)


if __name__ == "__main__":
    verify()
