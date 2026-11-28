"""十九点壳、可改边的平方弦长上限为 724 时，单位荷两钟四十二改边解在剩余全局重相下唯一。"""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb


def contract(vertices, pairs):
    labels = {v: v for v in vertices}
    for i,j in pairs:
        a,b = labels[i], labels[j]
        labels = {k: (min(a,b) if v in (a,b) else v) for k,v in labels.items()}
    return labels


def verify(replay):
    shell = run_path(str(replay/"ch39-arc-relative-determinant.py"))
    edges = shell["edges"]
    inner = set(shell["inner"])
    charges = [int(site in inner) for site in shell["sites"]]
    components = contract(range(19), ((i,j) for i,j,d,_ in edges if d>724))
    assert len(set(components.values())) == 13
    counts = defaultdict(Counter)
    for i, j, chord, _ in edges:
        a, b = int(chord==724), int(chord==804)+charges[i]-charges[j]
        if chord > 724:
            assert (a,b)==(0,0)
        left, right = components[i], components[j]
        if left > right:
            left, right, a, b = right, left, -a, -b
        counts[left,right][a,b] += 1

    labels = set().union(*(c.keys() for c in counts.values()))
    heights = [acb.zeta_zero(n).imag for n in (2345666,964890)]
    angles = [(shell["theta"](d)-shell["theta"](244))*(heights[0]-heights[1])
              for d in (724,804)]
    for (a,b),(c,d) in combinations(labels,2):
        assert abs(acb(0,(a-c)*angles[0]+(b-d)*angles[1]).exp().arg()) > arb(".09")

    zero = (0,0)
    unchanged = sum(c[zero] for c in counts.values())
    gains = {edge: max(c.values())-c[zero] for edge,c in counts.items()
             if edge[0]!=edge[1] and max(c.values())>c[zero]}
    third = {(3,11):15,(2,12):13,(0,18):9}
    assert gains == {edge:1 for edge in third}
    strict = {edge for edge,c in counts.items() if edge[0]!=edge[1]
              and c[zero]>max((n for label,n in c.items() if label!=zero),default=0)}
    reserved = set(third)
    for (i,j), k in third.items():
        legs = {tuple(sorted((i,k))),tuple(sorted((k,j)))}
        assert len(legs)==2 and legs <= strict and legs.isdisjoint(reserved)
        reserved.update(legs)
    forced = strict-reserved
    assert len(edges)==135 and len(edges)-unchanged==42 and len(forced)==51
    assert len(set(contract(set(components.values()), forced).values()))==1
    print("分量数 13，严格零差的额外分量对 51，相等图连通，最少改边数 42")


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
