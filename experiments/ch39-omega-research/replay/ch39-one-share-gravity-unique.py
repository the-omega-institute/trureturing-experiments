"""第五壳第16线的既定两个端点：全部物理标尺与份额上的引力解唯一性。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def divided(coeff,s,degree):
    # 在源标尺的基底上表示目标响应除以 s^degree，先合并确定的标尺幂。
    entries = coeff.entries()
    scaled = [value*s**(2-i-degree) for i,value in enumerate(entries[:-1])]
    scaled[2] += entries[-1]*s.log()/s**degree
    return arb_mat([[*scaled, entries[-1]/s**degree]])


def verify():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,change,share = repair["setup"]()
    gauge = load("one-share-repaired-gauge")
    targets = (*gauge["TARGETS"], *gauge["FORWARD"])
    _,lib,bgs = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(repair["SOURCE"],*targets),tight_neutral=True)
    model = lib["chiral"]["source"]["mean"]["inputs"]()[0]
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    scale = arb(7)/1000
    weights = sum(abs(w) for _,_,w in edges)
    n = kernel["ORDER"]+1
    factor = 2/(4*arb.pi()*factorial(n))
    min_x = min(a.lower() for a,_,_ in edges)*scale/10
    assert all(max(z.upper() for _,_,z in bg) < min(arb(1),2*min_x) for bg in bgs)
    physical_d = (-share).union(1-share)
    envelope = [(row[0], row[1], max(bg[i][2].upper() for bg in bgs[1:]))
                for i,row in enumerate(bgs[0])]
    def loss(bg,s):
        result = arb(0)
        for _,b,w in edges:
            x = scale*s*b
            result += abs(w)*x*sum(abs(a)*(1-lib["expint"](2,z/x)) for _,a,z in bg)
        return result/(2*arb.pi())
    # E2 的无质量项在源的已知 K=1 上归一，再加入第16线改变的单一份额。
    massless = sum(a for _,a,_ in bgs[0])/(2*arb.pi())
    slope = massless*change[1,0]*kernel["S"]
    small = (1+loss(bgs[0],1)+abs(physical_d*slope))/10+loss(envelope,arb(1)/10)
    print("s≤0.1 时牛顿矩上界",small.upper().str(10),flush=True)
    assert small < arb(".95317")
    neighborhood = arb(0,arb(".002"))
    local_s = arb(1)+neighborhood
    contraction_bounds = []
    for number,rows in zip(targets,bgs[1:],strict=True):
        series = kernel["series"](rows,fields)
        coeff = [arb_mat([[int(i==k) for i in range(5)]])*series for k in (0,1)]
        tails = [factor*sum(abs(2*a) for _,a,_ in rows),factor*sum(abs(g) for g,_,_ in rows)]
        local_args = (coeff,moments,change,[t*weights for t in tails],[2*t*abs(edges[16][2]) for t in tails])
        eps = arb("1e-122")
        # z/(s x)<2 支持整个范围。总误差对全部物理 d 一致。
        key = 1
        while scale*key/10 < 1:
            box = (arb(key)/10).union(arb(2*key)/10)
            key *= 2
            kar,war = [divided(row,box,i+1) for i,row in enumerate(coeff)]
            werr = (tails[1]*weights*2**n/box**2).upper()
            dw = (war*change)[0,0]
            assert not dw.contains(0)
            d = (eps/box**2-(war*moments)[0,0]+arb(0,werr))/dw
            kerr = (tails[0]*weights*2**n/box).upper()
            inv_s = (kar*(moments+d*change))[0,0]+arb(0,kerr)
            assert inv_s>0
            possible = 1/inv_s
            if possible.overlaps(box):
                assert local_s.contains(possible) and neighborhood.contains(d)
        # 区间内只用雅可比收缩给至多一根；存在性用同样的局部根证书重放。
        assert max(z.upper() for _,_,z in rows)<arb(".0069")*min(a.lower() for a,_,_ in edges)
        root,_ = repair["repair_box"](local_args,share)
        assert local_s.contains(root[0]) and neighborhood.contains(root[1])
        _,jac = repair["bounds"](local_s,neighborhood,*local_args)
        _,middle = repair["bounds"](arb(1),arb(0),*local_args)
        inv = arb_mat(2,2,[v.mid() for v in middle.inv().entries()])
        assert not inv.det().contains(0)
        contraction = max(repair["slopes"](arb_mat([[1,0],[0,1]])-inv*jac))
        assert contraction<1
        contraction_bounds.append(contraction)
    assert max(contraction_bounds) < arb(".00802")
    print("全物理域唯一的修复根数、小盒统一收缩常数上界",
          len(targets),max(contraction_bounds).str(10),flush=True)


if __name__=="__main__":
    verify()
