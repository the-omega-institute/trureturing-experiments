"""Shortest powers of the native shell loop at the requested scalar precision."""
from bisect import bisect_left
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
model = run_path(str(root/"ch39-clock-schur-mediators.py"))
B, root_cell, base = [model[k] for k in ("B","root_cell","base")]
word = base["word"]
unit = word["DEN"]


def ordered(point):
    (l,h),(a,b) = point
    if a >= 0:
        return (0,-h), (0,-l)
    assert b < 0
    return (1,l), (1,h)


def lattice(value):
    return value.lo*unit//1, -((-value.hi*unit)//1)


def first_power(phase, coefficient, bits, block, ceiling):
    tol = B(Q(1,2**bits))/coefficient
    ends = []
    for c in (phase+tol, phase-tol):
        sine = root_cell((B(1)-c*c)/B(17))
        ends.append((lattice(c),lattice(sine)))
    arcs = (ends, [word["conjugate"](v) for v in ends[::-1]])
    babies = sorted((*ordered(word["power"](k)),k) for k in range(block))
    assert all(u[1]<v[0] for u,v in zip(babies,babies[1:]))
    lows = [u[0] for u in babies]
    def position(point):
        lo,hi = ordered(point)
        pos = bisect_left(lows,lo)
        assert pos == len(babies) or hi < babies[pos][0]
        assert pos == 0 or babies[pos-1][1] < lo
        return pos
    for start in range(0,ceiling+1,block):
        rot = word["power"](-start)
        found = []
        for a,b in arcs:
            left = word["multiply"](rot,a)
            right = word["multiply"](rot,b)
            first,last = position(left),position(right)
            if ordered(left)[0] < ordered(right)[0]:
                subset = babies[first:last]
            else:
                subset = babies[first:]+babies[:last]
            found.extend(start+v[2] for v in subset)
        if found:
            k = min(found)
            assert k <= ceiling
            lo,hi = word["power"](k)[0]
            residual = coefficient*(B(Q(lo,unit),Q(hi,unit))-phase)
            assert max(abs(residual.lo),abs(residual.hi)) < Q(1,2**bits)
            return k,max(abs(residual.lo),abs(residual.hi))
    raise AssertionError("No word through the certified search horizon")


def main():
    ids, z, p, remainder, k, _, theta = model["calibration"]()
    m = B(1)/sum((v/x for x,v in zip(z,p)),B(0))
    m2_over_m = sum((v*x for x,v in zip(z,p)),B(0))
    mixed = B(2)*theta*root_cell(p[0])*(root_cell(p[1])+B(k)*root_cell(remainder))
    coefficient = m*mixed
    phase = B(1) - (m2_over_m-m)/coefficient
    for bits,block,ceiling in ((32,4096,14487319),(36,16384,458390681)):
        n,error = first_power(phase,coefficient,bits,block,ceiling)
        assert n==ceiling
        print("At",bits,"bits the least absolute shell power is",n,"; residual",float(error),flush=True)


if __name__=="__main__":
    main()
