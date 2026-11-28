"""三条谱线各自选位、其余 358 条共享插值时，检验全部例外选择的色逆耦合。"""
from collections import Counter
from functools import cache
from itertools import product
from math import comb
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb


ROOT = Path(__file__).resolve().parent
COARSE = 12
EXPECTED = {(15,285,323): 18, (16,17,18): 54, (285,304,323): 11, (285,323,342): 9}


def small_scale(edges, rows, heat, cut, unit, lower):
    pi = arb.pi()
    k0 = sum(alpha for _,alpha,_ in rows)/(2*pi)
    assert k0 < 0
    endpoints = [[w*a for a,_,w in edges], [w*b for _,b,w in edges]]
    minima = [v.min(w) for v,w in zip(*endpoints)]
    error = sum((-(w*alpha)).max(0)*cut*b*(1-heat(2,z/(cut*b)))
                for _,b,w in edges for _,alpha,z in rows)/(2*pi)
    bounds = []
    for endpoint in endpoints:
        corrections = sorted(lower(m-e) for m,e in zip(minima,endpoint))
        first_moment = sum(endpoint)+arb(sum(corrections[:3]))/unit
        bounds.append((k0*cut*first_moment).max(0)+error)
    assert max(v.upper() for v in bounds) < arb(".977")
    print("小标尺牛顿上界，两种公共端点",*[v.str(12) for v in bounds])


def coarse_boxes(tables, offset, limit):
    n = len(tables[0])
    low, high, colour = tables
    totals = (low.sum(axis=0), high.sum(axis=0), offset+colour.sum(axis=0))
    corrections = (low.min(axis=1)[:,None]-low, high.max(axis=1)[:,None]-high,
                   colour.max(axis=1)[:,None]-colour)
    counts = [0,0,0]
    survivors = []
    for i in range(n-2):
        j,k = np.triu_indices(n-i-1,k=1)
        j,k = j+i+1,k+i+1
        counts[0] += len(j)
        loose = [total+d[i]+d[j]+d[k] for total,d in zip(totals,corrections)]
        pairs,common = np.nonzero((loose[0]<=0)&(loose[1]>=1)&(loose[2]>limit))
        counts[1] += len(np.unique(pairs))
        counts[2] += len(common)
        # 成批展开三个独立的十二分区，各轴以同一批次的例外线号取数。
        for start in range(0,len(common),256):
            group = pairs[start:start+256]
            t = common[start:start+256]
            j_batch,k_batch = j[group],k[group]
            possible = True
            for field,(table,total) in enumerate(zip(tables,totals)):
                shared = total[t]-table[i,t]-table[j_batch,t]-table[k_batch,t]
                value = (shared[:,None,None,None]+table[i,None,:,None,None]
                         + table[j_batch,None,:,None]+table[k_batch,None,None,:])
                if field == 0:
                    possible = value <= 0
                elif field == 1:
                    possible &= value >= 1
                else:
                    possible &= value > limit
            for batch,ti,tj,tk in np.argwhere(possible):
                lines = (i,int(j_batch[batch]),int(k_batch[batch]))
                survivors.append((lines,(int(t[batch]),int(ti),int(tj),int(tk))))
    assert counts == [comb(n,3),65336,776946]
    assert Counter(lines for lines,_ in survivors) == EXPECTED
    assert all(cell[0] == 0 for _,cell in survivors)
    print("三元组总数、初筛后组数及公共区间数",counts,"粗筛余盒",len(survivors),flush=True)
    return survivors


def refine(survivors, edges, at, offset, limit):
    @cache
    def line(i, part, count):
        return at(edges[i],arb(part)/count,arb(part+1)/count)

    @cache
    def common(lines, part, count):
        rows = [line(i,part,count) for i in range(len(edges)) if i not in lines]
        return tuple(sum(axis) for axis in zip(*rows))

    counts = Counter()
    permitted = Counter()
    worst = None
    remaining = 0
    for lines, cell in survivors:
        for shifts in product(range(2),repeat=4):
            t, *positions = [2*k+s for k,s in zip(cell,shifts)]
            columns = [common(lines,t,2*COARSE)]
            columns.extend(line(i,pos,2*COARSE) for i,pos in zip(lines,positions))
            lo,hi,colour = map(sum,zip(*columns))
            counts[lines] += 1
            if lo<=0 and hi>=1:
                permitted[lines] += 1
                worst = max(worst,offset+colour) if worst is not None else offset+colour
                if colour+offset > limit:
                    remaining += 1
    assert counts == {lines: n*16 for lines,n in EXPECTED.items()}
    assert permitted == {(15,285,323):103, (16,17,18):431, (285,304,323):42, (285,323,342):41}
    assert worst == -130751447 and remaining == 0
    print("四组粗盒数、二分后真空候选数",[(lines,EXPECTED[lines],permitted[lines]) for lines in EXPECTED])
    print("二分总数、过真空的子盒、其中最大色支上界分子、尚余盒数",
          sum(counts.values()),sum(permitted.values()),worst,remaining)


def verify():
    def load(name):
        return run_path(str(ROOT/f"ch39-{name}.py"))

    edges,rows,heat,_ = load("two-copy-gram-cutoff-two-exceptions")["model_inputs"]()
    two = load("two-exception-colour-positivity")
    cut = arb(1)/625
    unit = two["UNIT"]
    assert len(edges) == 361 and two["PARTS"] == COARSE
    assert all(0<a and (a<b or a==b) and b<=1 for a,b,_ in edges)
    small_scale(edges,rows,heat,cut,unit,two["lower"])
    offset,at = two["line_bounds"](rows,cut)
    tables = two["interval_tables"](edges,at,COARSE)
    limit = -8*unit
    largest = max(abs(int(z)) for table in tables for z in table.flat)
    assert (len(edges)+10)*largest+abs(offset)+abs(limit) < 2**62
    assert arb("1e-122")/cut**2 < arb(1)/unit
    assert 0<7*(arb("1.22089e19")/10**16).log()/(2*arb.pi())<8
    survivors = coarse_boxes(tables,offset,limit)
    refine(survivors,edges,at,offset,limit)


if __name__ == "__main__":
    verify()
