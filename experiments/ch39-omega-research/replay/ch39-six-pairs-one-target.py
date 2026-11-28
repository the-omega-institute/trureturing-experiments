"""固定 2345712 的六种两线选择与源重谱的行列式排序。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
SOURCE = 2345698
TARGET = 2345712
PAIRS = tuple(combinations((16,285,301,323),2))


def verify():
    local = run_path(str(ROOT/"ch39-two-source-share-crossings.py"))
    global_check = run_path(str(ROOT/"ch39-two-share-global-forward.py"))["check_pair"]
    seeds = [{"heavy":n,"lines":pair} for n in (SOURCE,TARGET) for pair in PAIRS]
    bounds = {(16,285):(47241054,47241056), (285,301):(335320127,335320129),
              (16,301):(2605790594,2605790606)}
    ranks = {}
    energies = {}
    for n,pair,root,energy,minimum,failed,read,*_ in local["certified_pairs"]((SOURCE,TARGET),seeds):
        if n==SOURCE:
            _,jac = read(arb(1),arb(0),arb(0))
            ranks[pair] = abs(jac.det())
            continue
        visited,small = global_check(n,pair,root,read,(TARGET,),denominator=20)
        print(pair,"全域解数、数学根能标、最小逆耦合、越界条件、标尺盒数、小标尺 K 上界",
              int(not failed),energy.str(17),minimum.str(13),failed,visited,small.str(9),flush=True)
        assert small<arb(".751464")
        if pair in bounds:
            lo,hi = bounds[pair]
            assert not failed and minimum>arb("19.83") and lo<energy<hi
            energies[pair] = energy
        else:
            assert "能标上限" in failed and "正性" in failed
            assert ("第二份额上限" in failed) == (pair==(301,323))
    ranked = sorted(ranks,key=lambda k:ranks[k].mid(),reverse=True)
    assert all(ranks[i]>ranks[j] for i,j in zip(ranked,ranked[1:]))
    print("源重谱下的雅可比绝对行列式从大到小",[(k,ranks[k].str(11)) for k in ranked],flush=True)
    assert ranked == [(16,301),(285,301),(16,285),(16,323),(301,323),(285,323)]
    rank_ratio = ranks[ranked[0]]/ranks[ranked[1]]
    spread = energies[16,301]/energies[16,285]
    assert arb("4.34026")<rank_ratio<arb("4.34028")
    assert arb("55.1594")<spread<arb("55.1595")
    print("源端最大与次大行列式之比、最高与最低物理能标之比",rank_ratio.str(12),spread.str(12),flush=True)


if __name__=="__main__":
    verify()
