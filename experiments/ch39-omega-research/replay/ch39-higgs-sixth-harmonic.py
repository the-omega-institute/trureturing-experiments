"""Six harmonics impose a two-tail bound on atom mass for the optimal Higgs modes."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb

prior = run_path(str(Path(__file__).with_name("ch39-higgs-arc-weights.py")))
arc, q = prior["arc"], prior["q"]
fifth = arc["fifth"]
powers = [0, 1, 4, 9, 16, 25, 36]


def double(v):
    return [(v, False), (v, True)]


def polynomial(coeffs, value):
    return sum(c*value**k for c,k in zip(coeffs,powers))


def expected(coeffs, moments):
    return sum(a*b for a,b in zip(coeffs,moments))


def groups(weights, sizes, table):
    result = [arb(0) for _ in range(1+max(max(row) for row in table))]
    for weight, row in zip(weights,table):
        assert len(row)==len(sizes)
        for col,n in zip(row,sizes):
            result[col] += weight*q(F(n,sum(sizes)))
    return result


def sixth(nodes, weights):
    x = q(F(1,2**27))
    z = acb((1-x*x)/(1+x*x),2*x/(1+x*x))
    N = 2**53
    ratio = (1-(z**6).real)/(1-z.real)
    return sum(w*(1-ratio*(1-v**q(F(1,N))))**N for w,v in zip(weights,nodes))


def verify():
    weights = prior["mode_weights"]()
    assert min(weights)/4 > q(".01432")
    a,b,c,H,D = map(q,(".488",".8115",".8922",".949",".491"))
    interpolate = fifth["constrained_polynomial"]
    hi = interpolate(double(a)+double(b)+double(c),powers,1)
    interior = interpolate([(D,False),*double(b),*double(c),(H,False)],powers,1)
    lo = interpolate([*double(b),*double(c),(H,False)],powers[:-1],-1)
    assert polynomial(hi,H)>0 and polynomial(lo,D)>0
    old = (
        ([1]*4, ((4,2,2,4),(0,3,4,2),(2,4,4,1)), (501,771,838,882,902)),
        ([1,3,3,6,6], ((2,2,4,4,2),(4,1,0,4,4),(4,4,3,3,4)), (495,564,823,862,901))
    )
    new = (
        ([1]*5, ((4,4,2,3,4),(5,4,0,4,4),(3,4,1,4,4)), (472,765,814,833,894,950)),
        ([1,3,3,6,6], ((2,4,4,2,3),(5,2,0,4,4),(2,2,2,1,4)), (442,756,838,887,903,965))
    )
    two = fifth["two"]
    scale = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/scale)**2/(8*arb.pi())
    for zero in range(11485,11491):
        target = two["harmonics"](acb.zeta_zero(zero).imag,power,190,6)
        moments = fifth["moment_boxes"](target,powers)
        assert expected(hi,moments)/polynomial(hi,H) < q(".01363")
        assert expected(interior,moments) > q(".0000104")
        assert expected(lo,moments)/polynomial(lo,D) < q(".01361")
        for sizes,table,seeds in old+new:
            mixed = groups(weights,sizes,table)
            initial = [q(F(n,1000)) for n in seeds]
            nodes = arc["upper_certificate"](initial,target[:len(seeds)],2**53,mixed)
            if len(seeds)==5:
                shortfall = target[5]-sixth(nodes,mixed)
                lo_err,hi_err = ((".0001567",".0001572") if len(sizes)==4 else (".0001813",".0001816"))
                assert q(lo_err)<shortfall<q(hi_err)
        print(zero,"both fifth-order predictions excluded; six harmonics certified at five and nineteen labels.")


if __name__ == "__main__":
    verify()
