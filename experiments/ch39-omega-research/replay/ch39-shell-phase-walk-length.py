"""Shortest shell walks that reproduce a threshold-compatible multiple of a high-zero phase."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
PHASE = run_path(str(ROOT/"ch39-balanced-mass-phase.py"))
a = PHASE["a"]


def walks():
    # Reverse edges conjugate the phase. All 3 diamond components have this gauge.
    edges = [(0,1,0),(1,2,0),(1,3,0),(0,2,1),(0,3,-1)]
    neighbors = defaultdict(list)
    for x,y,p in edges:
        neighbors[x].append((y,p))
        neighbors[y].append((x,-p))
    counts = defaultdict(int)
    for start in range(4):
        states = {(start,0):1}
        for length in range(1,269):
            next_states = defaultdict(int)
            for (node,charge),count in states.items():
                for end,power in neighbors[node]:
                    next_states[end,charge+power] += count
            states = next_states
            for charge in (131,132,133,134):
                counts[length,charge] += 3*states.get((start,charge),0)
    return counts


def verify():
    ctx.prec = 150
    counts = walks()
    for n,kind,r,v,state in PHASE["MATCH"]["matched_points"]():
        _,minimum,_ = PHASE["phase_bounds"](r,1/state[0])
        theta = 2*(1/(2*acb.zeta_zero(n).imag)).atan()
        ratio = 2*(minimum/2).asin()/theta
        m = int(ratio.ceil().unique_fmpz())
        assert m == (131 if n<=11486 else 132)
        # t=1 is feasible for both mass branches; its chord is an interior upper witness.
        d = (a("13/10")/state[0]).sqrt()*r
        for j in (m,m+1,m+2):
            chord = 2*(j*theta/2).sin()
            assert minimum < chord < d*d
        lo = min(length for (length,j),count in counts.items() if j>=m and count)
        clock = min(length for (length,j),count in counts.items() if j>=m and count and length%3==1)
        assert (lo,clock)==((263,265) if m==131 else (264,268))
        first = {j:counts[clock,j] for j in range(m,135) if counts[clock,j]}
        print(n,kind,"positive exponent",m,"first lengths",lo,clock,"positive oriented counts",first)
    b = F(3,13)**3
    q = (1+2*b)/19
    for length in (265,268):
        bound = 19*q**length/(1-q**3)
        assert bound < F(1,10**330)
        print(length,"absolute trace remainder upper bound",a(bound).str(8))


if __name__ == "__main__":
    verify()
