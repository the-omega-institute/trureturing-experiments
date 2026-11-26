"""Read the shell mass product from unresolved spin contributions and split families."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from math import prod

from flint import ctx, fmpq, fmpq_mat, fmpq_poly

CURRENT = run_path(str(Path(__file__).with_name("ch39-complementary-current.py")))
a = CURRENT["a"]


def fermion_weights(slopes):
    for t, b in slopes.items():
        if t >= 5:
            # The three vector contributions are multiples of eight. Here only
            # the three bridge families may also contribute a Dirac onset.
            b = (b+3)%8-3
            assert -3 <= b <= 3
        if b:
            yield t, b


def collision_check():
    # The product identity applies for unequal per-family ratios as well.
    E, D, C = F(2), F(3), F(9,2)
    vectors = [(F(6),-8),(F(21),-8),(F(26),8)]
    external = [(E,1),(D,-2),(C,1)]
    heavy = [F(6), F(21), C]
    light = [F(26), F(6), D]
    slopes = defaultdict(int)
    for t, b in vectors+external+[(t,1) for t in heavy]+[(t,-1) for t in light]:
        slopes[t] += b
    decoded = dict(fermion_weights(slopes))
    expected = defaultdict(int)
    for t, b in external+[(t,1) for t in heavy]+[(t,-1) for t in light]:
        expected[t] += b
    assert decoded == {t:b for t,b in expected.items() if b != 0}
    assert sum(decoded.values()) == 0
    assert prod(t**b for t,b in decoded.items()) == prod(heavy)/prod(light)
    # Reconstruct the complete slope data from 72 moments, including the
    # constant scalar and vector terms and coincident simple/double poles.
    constants = defaultdict(F)
    for t, weight in vectors:
        constants[t] += F(5, 16)*weight
    for i in range(1, 13):
        constants[F(i, 4)] += F(i-6, 24)
    nodes = set(slopes) | set(constants)
    start = 5
    moments = [sum((constants[t]+(j+1)*slopes[t])/t**j for t in nodes)
               for j in range(start, start+72)]
    tail = [fmpq(v.numerator, v.denominator) for v in moments]
    cap = len(tail)//2
    rank = fmpq_mat([[tail[i+j] for j in range(cap)] for i in range(cap)]).rank()
    block = fmpq_mat([[tail[i+j] for j in range(rank)] for i in range(rank)])
    coeff = block.solve(fmpq_mat([[-v] for v in tail[rank:2*rank]])).entries()
    den = fmpq_poly([1,*reversed(coeff)])
    num = (den*fmpq_poly(tail[:rank])).truncate(rank)
    second = den.derivative().derivative()
    recovered = {}
    for factor, multiplicity in den.factor()[1]:
        if multiplicity == 2:
            t = -factor[0]/factor[1]
            b = 2*t**(start-2)*num(t)/second(t)
            assert b.q == 1
            recovered[F(int(t.p), int(t.q))] = int(b)
    assert dict(fermion_weights(recovered)) == decoded
    print("72-moment collision example: rank",rank,"product",prod(t**b for t,b in decoded.items()))
    assert all((b+8*k+3)%8-3 == b for b in range(-3,4) for k in range(-2,3))


def verify():
    ctx.prec = 232
    collision_check()
    epsilon_table = CURRENT["GRID"]["MATCH"]["SOURCE"]["RATIO"]["high_zero_readings"]
    readings = [(n, eps**(-6)) for n, _, eps in epsilon_table()]
    gaps = [1-right/left for (_,left),(_,right) in zip(readings,readings[1:])]
    assert all(v > gaps[-1] for v in gaps[:-1])
    assert a("9.65477e-19") < gaps[-1] < a("9.65478e-19")
    for (_, left), (_, right) in zip(readings, readings[1:]):
        e = a("4.82e-19")
        assert right*(1+e) < left*(1-e)
    # Finite per-node errors, provided the signed integer weights are correct.
    e = a("4.82e-20")
    worst = ((1+e)/(1-e))**5
    assert worst-1 < a("4.821e-19")
    assert all(right*worst < left/worst for (_,left),(_,right) in zip(readings,readings[1:]))
    print("Shell product range", readings[-1][1].str(29), readings[0][1].str(29))
    print("Closest fractional product gap", gaps[-1].str(18))


if __name__ == "__main__":
    verify()
