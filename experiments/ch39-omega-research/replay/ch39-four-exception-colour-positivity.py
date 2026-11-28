"""用剩余谱线的最大贡献逐组排除四个单位置例外。"""
import ctypes
from math import comb
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


def read_limits(tables, offset):
    low, high, colour = tables
    # 先编码 W/c² 的两侧必要界，再编码 4U₃+mW/c² 的逐线必要上界。
    fields = [-low, high]
    offsets = [0, -1]
    for multiplier in MULTIPLIERS:
        vacuum = high if multiplier >= 0 else low
        fields.append(4*colour+multiplier*vacuum)
        offsets.append(4*(offset+8*UNIT)+max(0,-multiplier)-1)
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


if __name__ == "__main__":
    verify()
