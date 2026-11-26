"""前四项时间与参考味驻点给第五项的全域锐界。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx


shared = run_path(str(Path(__file__).with_name("ch39-shared-flavor-clock-weights.py")))


def time_kernel(x, order):
    # None 表示有限坐标 u=1/(1+x) 的左端点 u=0。
    return arb(0) if x is None else 1/(1+order**2*x)


def gradient(x, lambdas, slopes):
    if x is None:
        return sum(s/z for s,z in zip(slopes,lambdas,strict=True))
    return sum(s*x/(1+z*x) for s,z in zip(slopes,lambdas,strict=True))


def finite_weights(nodes, moments):
    matrix = arb_mat([[x**j for x in nodes] for j in range(3)])
    tilted = matrix.solve(arb_mat([[v] for v in moments[:3]])).entries()
    weights = [w*prod(1+j*j*x for j in range(1,5)) for w,x in zip(tilted,nodes,strict=True)]
    return tilted, weights


def gaussian(end, moments):
    tilted = [end*moments[i]-moments[i+1] for i in range(4)]
    gram = arb_mat([[tilted[0],tilted[1]],[tilted[1],tilted[2]]])
    p, q = gram.solve(arb_mat([[-tilted[2]],[-tilted[3]]])).entries()
    nodes = [(-q+s*(q*q-4*p).sqrt())/2 for s in (-1,1)]+[end]
    _, weights = finite_weights(nodes,moments)
    return nodes, weights


def endpoints(a, moments):
    b = (moments[3]-a*moments[2])/(moments[2]-a*moments[1])
    nodes = [arb(0), a, b]
    tilted, weights = finite_weights(nodes,moments)
    last = 576*(moments[4]-sum(w*x**4 for w,x in zip(tilted,nodes,strict=True)))
    return [*nodes,None], [*weights,last]


def integrate(rule, function):
    nodes, weights = rule
    return sum(w*function(x) for w,x in zip(weights,nodes,strict=True))


def stationary_rule(rule, moments, lambdas, slopes, lo, hi, orientation):
    def score(v):
        return orientation*integrate(rule(v,moments),lambda x: gradient(x,lambdas,slopes))
    left, right = arb(lo), arb(hi)
    # 分母及平方根在整个起始区间严格可定义，夹逼连续的驻点方程。
    nodes, weights = rule(left.union(right),moments)
    assert all(x.is_finite() for x in nodes if x is not None) and all(w.is_finite() for w in weights)
    assert score(left) < 0 < score(right)
    for _ in range(250):
        mid = (left+right)/2
        value = score(mid)
        if value > 0:
            right = mid
        elif value < 0:
            left = mid
        else:
            break
    assert right-left < arb("1e-15")
    return rule(left.union(right),moments)


def certify_dual(rule, lambdas, slopes, sign):
    nodes, _ = rule
    contacts = [(x,False) for x in nodes]+[(x,True) for x in nodes if x is not None and x>0]
    matrix, target = [], []
    for x, derivative in contacts:
        if derivative:
            matrix.append([0,*[-j*j/(1+j*j*x)**2 for j in range(1,5)],
                           sum(s/(1+z*x)**2 for s,z in zip(slopes,lambdas,strict=True))])
            target.append([-25/(1+25*x)**2])
        else:
            matrix.append([1,*[time_kernel(x,j) for j in range(1,5)],gradient(x,lambdas,slopes)])
            target.append([time_kernel(x,5)])
    coefficients = arb_mat(matrix).solve(arb_mat(target)).entries()
    xvar = arb_poly([0,1])
    rates = [1,4,9,16,25,*lambdas]
    factors = [1+v*xvar for v in rates]

    def without(i):
        return prod(factors[:i]+factors[i+1:])

    numerator = coefficients[0]*prod(factors)-without(4)
    numerator += sum(coefficients[j+1]*without(j) for j in range(4))
    numerator += coefficients[-1]*xvar*sum(s*without(5+j) for j,s in enumerate(slopes))
    finite_contacts = prod(xvar-x for x,_ in contacts if x is not None)
    quotient, remainder = divmod(numerator,finite_contacts)
    assert all(v.contains(0) for v in remainder)
    if nodes[-1] is None:
        # u=0 处的插值方程使八次首项精确为零。
        assert numerator[8].contains(0) and quotient[3].contains(0)
    assert all(sign*quotient[j] > 0 for j in range(3))
    return coefficients


def certify_maximum(rule, invariants, reference):
    # 用比值 R=D(c)/D(c₀) 及 log R 的二阶上界验证集体味作用严格取极大。
    total = arb_poly([0])
    correction = arb(0)
    delta = arb(0,1)-reference
    for x,w in zip(*rule,strict=True):
        e1,e2,e3 = invariants
        determinant = e3 if x is None else 1+e1*x+e2*x*x+e3*x**3
        normalized = (determinant*(1/determinant(reference)))(arb_poly([reference,1]))
        linear_factor = arb_poly(normalized.coeffs()[1:])(delta)
        if x is None:
            row = 1 + delta*linear_factor
            assert row > 0 and row < arb("1.01") and linear_factor > arb(".043")
            correction = w*arb(".043")**2/(2*arb("1.01")**2)
        total += w*normalized
    assert total[1].contains(0)
    assert arb_poly(total.coeffs()[2:])(delta) < correction


def verify():
    ctx.prec, ctx.cap = 768, 3
    times = shared["clock"]["clock"]["times"]()
    moments = [v/200**j for j,v in enumerate(shared["clock"]["clock"]["cleared_moments"](times,4))]
    invariants, reference, lambdas, slopes = shared["flavor_data"]()
    lower = stationary_rule(gaussian,moments,lambdas,slopes,"23834887","23834888",1)
    upper = stationary_rule(endpoints,moments,lambdas,slopes,".00097276998",".00097276999",-1)
    for rule,sign in ((lower,-1),(upper,1)):
        nodes, weights = rule
        finite = [x for x in nodes if x is not None]
        assert all(a<b for a,b in zip(finite,finite[1:])) and all(w>0 for w in weights)
        assert sum(weights).contains(1)
        for j in range(1,5):
            assert integrate(rule,lambda x:time_kernel(x,j)).overlaps(times[j])
        coefficients = certify_dual(rule,lambdas,slopes,sign)
        value = integrate(rule,lambda x:time_kernel(x,5))
        assert value.overlaps(sum(a*b for a,b in zip(coefficients[:5],times[:5],strict=True)))
        certify_maximum(rule,invariants,reference)
        slack = sum(abs(v) for v in coefficients[1:5])*arb(2)**-33 + abs(coefficients[5])*arb(2)**-38
        assert ((value+sign*slack)*2**22).floor().unique_fmpz() == 4040282
        print("下界" if sign < 0 else "上界", value.str(18), "允许扰动后的外扩",slack.str(12))
    assert (times[5]*2**22).floor().unique_fmpz() == 4040283
    print("第五项算术时间的参考整数与共同模型整数",4040283,4040282)


if __name__ == "__main__":
    verify()
