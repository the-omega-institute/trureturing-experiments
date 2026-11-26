"""十九点中性壳的预解式迹为第五个算术时间提供两种读出的效率差。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx, fmpz_mat


root = Path(__file__).parent
contrast = run_path(str(root / "ch39-five-time-source-contrast.py"))
model = contrast["model"]
shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))


def neutral_polynomial():
    laplacian = fmpz_mat(19,19)
    for i,j,_,_ in shell["geo"]["edges"]:
        laplacian[i,i] += 1
        laplacian[j,j] += 1
        laplacian[i,j] -= 1
        laplacian[j,i] -= 1
    coefficients = (-laplacian).charpoly().coeffs()
    assert len(coefficients)==20 and coefficients[0]==0 and all(c>0 for c in coefficients[1:])
    return arb_poly(coefficients[1:])


def efficiency(x, polynomial):
    return 1+x*polynomial.derivative()(x)/polynomial(x)


def fixed_ends(moments, low, high):
    tilted = [-moments[i+2]+(low+high)*moments[i+1]-low*high*moments[i] for i in range(4)]
    matrix = arb_mat([[tilted[0],tilted[1]],[tilted[1],tilted[2]]])
    p,q = matrix.solve(arb_mat([[-tilted[2]],[-tilted[3]]])).entries()
    nodes = [low,*[(-q+sign*(q*q-4*p).sqrt())/2 for sign in (-1,1)],high]
    vandermonde = arb_mat([[x**j for x in nodes] for j in range(4)])
    tilted_weights = vandermonde.solve(arb_mat([[v] for v in moments[:4]])).entries()
    weights = [w*prod(1+j*j*x for j in range(1,6)) for w,x in zip(tilted_weights,nodes,strict=True)]
    assert all(a<b for a,b in zip(nodes,nodes[1:])) and all(w>0 for w in weights)
    return nodes, weights


def verify_maximum(nodes, weights, invariants, c0):
    # weights 已经是味侧的归一概率。
    t = arb(0,1)-c0
    average = arb_poly([0])
    penalty = arb(0)
    e1,e2,e3 = invariants
    for x,w in zip(nodes,weights,strict=True):
        det = 1+x*e1+x**2*e2+x**3*e3
        ratio = (det*(1/det(c0)))(arb_poly([c0,1]))
        average += w*ratio
        if x>1000:
            divided = arb_poly(ratio.coeffs()[1:])(t)
            maximum = max(arb(1),(1+t*divided).upper())
            assert divided>0
            penalty += w*divided.lower()**2/(2*maximum**2)
    assert average[1].contains(0)
    assert arb_poly(average.coeffs()[2:])(t)<penalty


def verify():
    ctx.prec, ctx.cap = 768, 3
    shared = model["shared"]
    times = shared["clock"]["clock"]["times"]()
    invariants,c0,lambdas,slopes = shared["flavor_data"]()
    moments = contrast["moments_five"](times)
    poly = neutral_polynomial()
    low = arb("1e-6")
    rules = [fixed_ends(moments,low,arb(v)) for v in ("1e9","1e13")]

    def weighted_gradient(x):
        return efficiency(x,poly)*model["gradient"](x,lambdas,slopes)

    left, right = [model["integrate"](rule,weighted_gradient) for rule in rules]
    assert left < 0 < right
    fraction = right/(right-left)
    nodes = [low,*rules[0][0][1:],*rules[1][0][1:]]
    weights = [fraction*rules[0][1][0]+(1-fraction)*rules[1][1][0],
               *[fraction*w for w in rules[0][1][1:]], *[(1-fraction)*w for w in rules[1][1][1:]]]
    columns = [[1,*[model["time_kernel"](x,j) for j in range(1,6)],weighted_gradient(x)] for x in nodes]
    matrix = arb_mat(list(zip(*columns,strict=True)))
    actual = matrix.solve(arb_mat([[v] for v in [*times[:6],0]])).entries()
    assert all(v>0 and v.overlaps(w) for v,w in zip(actual,weights,strict=True))
    tilted = [v*efficiency(x,poly) for x,v in zip(nodes,actual,strict=True)]
    verify_maximum(nodes,[w/sum(tilted) for w in tilted],invariants,c0)
    print("两组有限求积的味导数与第一组份额",left.str(13),right.str(13),fraction.str(12))
    print("七节点模型的平均中性响应",sum(tilted).str(12))


if __name__ == "__main__":
    verify()
