"""实际正零点稀疏检出名单、所有更短名单的精确消零覆盖，以及任意单点的局限。"""

from json import loads
from math import comb
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb, arb_mat, ctx


INDICES = (2,5,10,12,14,21,23,30)
MULTIPLIERS = (-9925065,-16987316,14692113,7443624,6886669,-12476412,16269156,-15319644)


def single_angle(root, types):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    geometry = load("mincut-real-triangle-tomography")
    model = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    angles = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    x = arb("-.01046")
    y = x-angles["reference_angles"](flavor)[0]
    weights = load("height-independent-flavor-vacuum")["minimum_order_weights"](
        model,model["terms"](shell),x,y)
    angles = [geometry["unit"](pair) for clock,pair,_ in types if clock==(1,0)]
    assert len(angles)==14
    # 只用 x 钟的十四组。各行在一个网格点上提出一组正负候选，再以统一斜率界覆盖整圆。
    grid = 8192
    nodes = np.arange(grid)*2*np.pi/grid
    phases = [float(z.arg()) for z in angles]
    trial = sum(float(w)*np.cos(q*(nodes[:,None]+phases)) for q,w in weights.items())
    cs = [(2*arb.pi()*k/grid).cos() for k in range(grid)]
    sn = [(2*arb.pi()*k/grid).sin() for k in range(grid)]
    coeff = [[(q,w*(z**q).real,w*(z**q).imag) for q,w in weights.items()] for z in angles]
    def score(k,j):
        return sum(a*cs[q*k%grid]-b*sn[q*k%grid] for q,a,b in coeff[j])
    positives = [score(k,j) for k,j in enumerate(trial.argmax(axis=1))]
    negatives = [-score(k,j) for k,j in enumerate(trial.argmin(axis=1))]
    margin = min(z.lower() for z in positives+negatives)-arb.pi()/grid*sum(q*w for q,w in weights.items())
    cap = arb(23)/25*arb(2)**-20/14**2
    lower = (1-cap)*margin
    assert lower>arb(".12")
    print("任一角和允许高度上异号两组的读数绝对值下界", lower.str(16))


def masks_of_size(k):
    # 高低两段各十六位。使用整型位运算生成所有名单，精确检验包含关系。
    bits = [np.array([v for v in range(1<<16) if v.bit_count()==h],dtype=np.uint32)
            for h in range(k+1)]
    pieces = [((bits[k-h][:,None]<<16) | bits[h][None,:]).ravel() for h in range(k+1)]
    values = np.concatenate(pieces)
    assert len(values)==comb(32,k)
    return values


def verify(root):
    types,_,rows = run_path(str(root/"ch39-proton-detection-prefix.py"))["samples"](root,32)
    ctx.prec = 256
    denom = sum(map(abs,MULTIPLIERS))
    scores = [sum(z*rows[k-1][j] for z,k in zip(MULTIPLIERS,INDICES))/denom for j in range(32)]
    gap = min(z.lower() for z in scores)
    assert gap>arb(".00223109")
    print("检出编号及最坏最大模的保证下界",INDICES,gap.str(16))
    # 对这份名单，最坏面概率恰在八个轨道上取正，八个带符号读数同时取等。
    support = (6,7,10,11,17,23,25,26)
    signed = [[(1 if v>0 else -1)*x for x in rows[n-1]] for n,v in zip(INDICES,MULTIPLIERS)]
    basis = arb_mat([[r[j] for j in support]+[-1] for r in signed]+[[1]*8+[0]])
    p = basis.solve(arb_mat([[0]]*8+[[1]])).entries()
    dual = basis.transpose().solve(arb_mat([[0]]*8+[[-1]])).entries()
    assert min(p)>0 and min(dual[:8])>0 and (sum(dual[:8])-1).contains(0)
    assert (p[-1]+dual[-1]).contains(0)
    for j in range(32):
        if j not in support:
            assert sum(y*r[j] for r,y in zip(signed,dual))>p[-1]
    print("该八项名单最坏最大绝对读数的可达值",p[-1].str(19))
    single_angle(root,types)
    records = loads(Path(__file__).with_name("ch39-sparse-proton-cover.json").read_text())
    k = len(INDICES)-1
    uncovered = masks_of_size(k)
    minimum = arb(1)
    for nums,cols in records:
        assert 0<nums<2**32 and 0<cols<2**32
        Z = [i for i in range(32) if (nums>>i)&1]
        J = [i for i in range(32) if (cols>>i)&1]
        assert len(J)==len(Z)+1
        M = arb_mat([[1]*len(J)]+[[rows[i][j] for j in J] for i in Z])
        masses = M.solve(arb_mat([[1]]+[[0]]*len(Z))).entries()
        assert min(masses)>0
        minimum=minimum.min(min(masses))
        uncovered = uncovered[(uncovered & np.uint32(nums))!=uncovered]
    assert len(uncovered)==0
    print("全部较短名单的长度、名单数、相消分布证书数、见证最小正组质量",
          k,comb(32,k),len(records),minimum.lower().str(15))


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
