"""用剩余谱线的最大贡献逐组排除四个单位置例外。"""
import ctypes
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
from runpy import run_path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
from flint import arb


ROOT = Path(__file__).resolve().parent
SEGMENTS = 48
MULTIPLIERS = (-256,-128,-96,-80,-76,-72,-64,-48,-32,-16,-8,-4,-1,0,4,16,64,256)
UNIT = 2**16


def read_limits(tables, offset, extra=()):
    low, high, colour = tables
    # 先编码 W/c² 的两侧必要界，再编码 4U₃+mW/c² 的逐线必要上界。
    fields = [-low, high]
    offsets = [0, -1]
    for multiplier in MULTIPLIERS:
        vacuum = high if multiplier >= 0 else low
        fields.append(4*colour+multiplier*vacuum)
        offsets.append(4*(offset+8*UNIT)+max(0,-multiplier)-1)
    for field, constant in extra:
        fields.append(field)
        offsets.append(constant)
    largest = max(abs(int(x)) for table in fields for x in table.flat)
    assert (len(colour)+20)*largest + max(abs(z) for z in offsets) < 2**62
    values = np.ascontiguousarray(np.array(fields).transpose(1,2,0))
    totals = np.ascontiguousarray(values.sum(axis=0)+np.array(offsets))
    return values, totals


def pick_lines(values, totals):
    n, segments, fields = values.shape
    maxima = values.max(axis=1)
    candidates = []
    for t in range(segments):
        corrections = maxima-values[:,t]
        # 各不等式单独排序，只影响遍历顺序；所有剪枝均在整数域完成。
        dimension = np.maximum(1,np.abs(totals[t]))
        order = np.argsort(-(corrections/dimension).max(axis=1))
        gains = corrections[order]
        suffix = np.zeros((n+1,5,fields),dtype=np.int64)
        for j in reversed(range(n)):
            for k in range(1,5):
                suffix[j,k] = np.maximum(suffix[j+1,k],gains[j]+suffix[j+1,k-1])

        def visit(start, slots, budget, chosen):
            stop = n-slots+1
            potential = budget + gains[start:stop] + suffix[start+1:stop+1,slots-1]
            for index in np.flatnonzero(np.all(potential>=0,axis=1))+start:
                lines = (*chosen,int(order[index]))
                if slots == 1:
                    candidates.append((*sorted(lines),t))
                else:
                    visit(index+1,slots-1,budget+gains[index],lines)

        visit(0,4,totals[t],())
    return candidates


def pick_cells(program, values, totals, candidates):
    candidates = np.ascontiguousarray(candidates,dtype=np.int32)
    segments,fields = values.shape[1:]
    maxima = np.ascontiguousarray(values.max(axis=1))
    reductions = np.ascontiguousarray(values-maxima[:,None,:])
    assert reductions.max() == 0
    keep = np.zeros(len(candidates),dtype=np.int32)
    arguments = (values,maxima,reductions,totals,candidates)
    program(*(m.ctypes.data for m in arguments),len(candidates),fields,segments,keep.ctypes.data)
    return candidates[keep!=0]


def exclude_crossings(program, remaining, edges, rows, two):
    gauge = run_path(str(ROOT/"ch39-signed-gravity-unification.py"))
    inputs = run_path(str(ROOT/"ch39-chiral-minimal-neutral-gravity.py"))
    model, _, _, initial, _ = inputs["chiral"]["source"]["mean"]["inputs"]()
    lib = run_path(str(ROOT/"ch39-upward-crossing-cutoff-inverse-moment.py"))
    const, _, _, heavy = lib["inverse_moment_terms"](ROOT,model,initial,rows)
    pi = arb.pi()
    gap = (arb(10**16)/model["mz"]).log()-const
    assert gap > 0
    moment = [sum(b*z**k for b,z in heavy)/(4*pi*k*factorial(k)) for k in range(1,5)]
    remainder = sum(abs(b)*z**5 for b,z in heavy)/(4*pi*5*factorial(5))
    log_colour = -sum(b[2] for b in gauge["gauge_weights"](model))/(4*pi)
    k0 = sum(a for _,a,_ in rows)/(2*pi)
    active = set(int(i) for group in remaining for i in group[:4])
    shared = sorted(set(int(group[4]) for group in remaining))
    grid = [arb(i)/SEGMENTS for i in range(SEGMENTS+1)]
    # E₁(y) ≤ -log y + E₁(1)（y≤1）；对 y≥1 则用热核亏额 ≤1。
    assert arb(1).expint(1) < arb(".22")
    endpoint = Fraction(14,10000)
    fine_left, fine_right = Fraction(".004725"), Fraction(".01063125")
    batches = 0
    while endpoint < 1:
        boundary = min(1,endpoint*Fraction(3,2))
        if fine_left <= endpoint < fine_right:
            boundary = min(fine_right,endpoint*Fraction(21,20))
        low = arb(endpoint.numerator)/endpoint.denominator
        high = arb(boundary.numerator)/boundary.denominator
        offset, at = two["line_bounds"](rows,low)
        up = two["upper"]
        offset += up(log_colour*high.log())
        values = np.zeros((len(edges),SEGMENTS,6),dtype=np.int64)
        for i,(a,b,w) in enumerate(edges):
            gain, loss = arb(0),arb(0)
            # v(1-E₂(z/(cv))) 随 v 增、随 c 减，统一取 b 与 low。
            for _,alpha,z in rows:
                if z.is_zero():
                    continue
                y = z/(low*b)
                err = (y*(1+arb(".22")+(-y.log()).max(0))).min(1)*b/(2*pi)
                gain += (-w*alpha).max(0)*err
                loss += (w*alpha).max(0)*err
            for j in range(SEGMENTS) if i in active else shared:
                left, right = [(1-u)*a+u*b for u in grid[j:j+2]]
                inverse = (1/(high*right)).union(1/(low*left))
                # -Σ b Ein(z/x)/(4π)。截于四次，五阶余项对任意 y≥0 成立。
                response = sum((-1)**k*coef*inverse**k for k,coef in enumerate(moment,1))
                response += arb(0,(remainder*inverse**5).upper())
                kl,kr = k0*w*left,k0*w*right
                values[i,j] = [*at(edges[i],grid[j],grid[j+1]),
                               up(kl.max(kr)+gain),up(-kl.min(kr)+loss),up(w*response/gap)]
        lo,hi,colour,k_plus,k_minus,projection = values.transpose(2,0,1)
        assert 270*max(abs(int(v)) for v in values.flat) < 2**62
        extra = [(k_plus,up(-1/high)),(k_minus,up(1/low)),(projection,-UNIT)]
        terms,totals = read_limits((lo,hi,colour),offset,extra)
        assert not len(pick_cells(program,terms,totals,remaining)), (endpoint,boundary)
        endpoint = boundary
        batches += 1
    assert batches == 32
    print("三个候选的标尺闭区间、逐段公共插值区间、保留的指定能标正交汇",
          batches,len(remaining),0)


def verify():
    def load(name):
        return run_path(str(ROOT/f"ch39-{name}.py"))

    edges,rows,heat,_ = load("two-copy-gram-cutoff-two-exceptions")["model_inputs"]()
    two = load("two-exception-colour-positivity")
    three = load("three-exception-colour-positivity")
    cut = arb(14)/10000
    assert len(edges)==361 and comb(361,4)==695946630
    assert all(0<a and (a<b or a==b) and b<=1 for a,b,_ in edges)
    three["small_scale"](edges,rows,heat,cut,UNIT,two["lower"],exceptions=4)
    assert arb("1e-122")/cut**2 < arb(1)/UNIT
    assert 7*(arb("1.22089e19")/10**16).log()/(2*arb.pi()) < 8
    offset,at = two["line_bounds"](rows,cut)
    tables = two["interval_tables"](edges,at,SEGMENTS)
    assert 270*max(abs(int(x)) for table in tables for x in table.flat)<2**62
    values,totals = read_limits(tables,offset)
    candidates = pick_lines(values,totals)
    assert len(candidates)==20259 and len(set(x[:4] for x in candidates))==2147
    with TemporaryDirectory() as temp:
        shared = Path(temp)/"scan.so"
        source = Path(__file__).with_suffix(".cpp")
        subprocess.run(["c++","-std=c++17","-O3","-shared","-fPIC",str(source),"-o",str(shared)],check=True)
        library = ctypes.CDLL(str(shared))
        program = library.scan
        ptr,integer = ctypes.c_void_p,ctypes.c_int
        program.argtypes = [ptr]*5 + [integer]*3 + [ptr]
        program.restype = None
        remaining = pick_cells(program,values,totals,candidates)
        groups = {tuple(line[:4]) for line in remaining}
        assert groups == {(15,285,304,323),(15,285,323,342),(16,17,18,323)}
        assert len(remaining) == 17
        print("全域四元组、初筛组数及公共段数、再筛组数及公共段数",
              comb(361,4),2147,len(candidates),len(groups),len(remaining),flush=True)
        exclude_crossings(program,remaining,edges,rows,two)


if __name__ == "__main__":
    verify()
