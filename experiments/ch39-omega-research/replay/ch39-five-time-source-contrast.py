"""五个算术时间与味真空共用零点时，两种读出的相对效率锐下界。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx


root = Path(__file__).parent
model = run_path(str(root / "ch39-global-fifth-time.py"))
shared = model["shared"]


def moments_five(times):
    return [v/200**i for i,v in enumerate(shared["clock"]["clock"]["cleared_moments"](times,5))]


def compact_rule(moments):
    m = moments
    p,q = arb_mat([[m[1],m[2]],[m[2],m[3]]]).solve(arb_mat([[-m[3]],[-m[4]]])).entries()
    nodes = [arb(0)]+[(-q+sign*(q*q-4*p).sqrt())/2 for sign in (-1,1)]
    tilted = arb_mat([[x**j for x in nodes] for j in range(3)]).solve(arb_mat([[v] for v in m[:3]])).entries()
    weights = [w*prod(1+j*j*x for j in range(1,6)) for w,x in zip(tilted,nodes,strict=True)]
    tail = 14400*(m[5]-sum(w*x**5 for w,x in zip(tilted,nodes,strict=True)))
    return [*nodes,None],[*weights,tail]


def verify():
    ctx.prec, ctx.cap = 768, 3
    times = shared["clock"]["clock"]["times"]()
    invariants, c0, lambdas, slopes = shared["flavor_data"]()
    nodes, weights = compact_rule(moments_five(times))
    assert 0 < nodes[1] < nodes[2] and all(w>0 for w in weights)
    assert sum(weights).contains(1)
    for j in range(1,6):
        assert model["integrate"]((nodes,weights),lambda x:model["time_kernel"](x,j)).overlaps(times[j])
    d1,d2,d3 = [p.derivative()(c0) for p in invariants]
    assert d1 < 0 < d2 and d3 > 0
    crossing = (-d2+(d2*d2-4*d1*d3).sqrt())/(2*d3)
    assert crossing > nodes[2]
    g = [model["gradient"](x,lambdas,slopes) for x in nodes]
    ratio = -sum(w*v for w,v in zip(weights[:3],g[:3]))/(weights[3]*g[3])
    assert arb("18.49684338265") < ratio < arb("18.49684338267")
    tilted = [*weights[:3],weights[-1]*ratio]
    model["certify_maximum"]((nodes,[w/sum(tilted) for w in tilted]),invariants,c0)

    matrix, target = [], []
    for x,v in zip(nodes,g,strict=True):
        matrix.append([1,*[model["time_kernel"](x,j) for j in range(1,6)]])
        target.append([v*(ratio if x is None else 1)])
    for x in nodes[1:3]:
        matrix.append([0,*[-j*j/(1+j*j*x)**2 for j in range(1,6)]])
        target.append([sum(s/(1+z*x)**2 for s,z in zip(slopes,lambdas,strict=True))])
    coefficients = arb_mat(matrix).solve(arb_mat(target)).entries()
    var = arb_poly([0,1])
    factors = [1+v*var for v in [1,4,9,16,25,*lambdas]]

    def without(i):
        return prod(factors[:i]+factors[i+1:])

    majorant = coefficients[0]*prod(factors)+sum(coefficients[i+1]*without(i) for i in range(5))
    gradient = var*sum(s*without(5+i) for i,s in enumerate(slopes))
    contact = var*(var-nodes[1])**2*(var-nodes[2])**2
    quotient, remainder = divmod(majorant-gradient,contact)
    assert all(v.contains(0) for v in remainder) and all(v>0 for v in quotient)
    high = (majorant-ratio*gradient)(arb_poly([crossing,1]))
    assert high[8].contains(0) and all(high[i]>0 for i in range(8))
    print("五项精确的最小效率倍率下确界", ratio.str(16))


if __name__ == "__main__":
    verify()
