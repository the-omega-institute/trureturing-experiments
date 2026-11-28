"""五背景两条独立线的全部 64980 种选择，允许各自移动比较能标。"""
from itertools import combinations
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def verify():
    five = run_path(str(ROOT/"ch39-five-split-scale.py"))
    data = five["inputs"]()
    kernel = data[-1]
    build = five["load"]("three-background-calibration")["projection_terms"]
    minimum = five["load"]("shared-continuous-position-gap")["projected_minimum"]
    def lower(terms):
        return minimum(kernel["basis"],kernel["S"],terms,arb("1e-7"))
    directions = json.loads((ROOT/"ch39-five-pair-directions.json").read_text())
    pairs = list(combinations(range(361),2))
    for vector in directions:
        const,_,rows = build(*data,vector,None,five["SCALES"])
        norm = sum(abs(arb(n)/10**11) for vals in vector for n in vals[2:4])
        common = [lower(row) for row in rows]
        base = const+sum(common)
        need = {i for pair in pairs for i in pair}
        gain = {i:common[i]-sum(lower([term]) for term in rows[i]) for i in need}
        remaining,eliminated = [],[]
        for i,j in pairs:
            bound = (base-gain[i]-gain[j])/norm
            if bound>0:
                eliminated.append(bound.lower())
            else:
                remaining.append((i,j))
        pairs = remaining
        assert min(eliminated)>arb("0.01001")
        print("筛后余下的线对数、本步排除的规范差下界",len(pairs),min(eliminated).str(12),flush=True)
    print("剩余线对",pairs,flush=True)
    assert pairs == [(0,15),(0,18),(0,285),(0,323)]


if __name__=="__main__":
    verify()
