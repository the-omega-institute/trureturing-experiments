"""以七阶排列多面体的迹截面判定同一平凡块能否同时提供轻上、粲两个低方向。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb


def block_levels(block):
    a,b,c = block
    gap = ((a-b)**2+4*abs(c)**2).sqrt()
    return (a+b-gap)/2,(a+b+gap)/2


def block(coefficients, entries):
    return tuple(sum(w*z for w,z in zip(coefficients,row,strict=True)) for row in entries)


def join(first,second,t):
    return tuple(a+t*(b-a) for a,b in zip(first,second,strict=True))


def clipped_minimum(weights,entries,trace_cap):
    vertices = {}
    kept = set()
    lower = []
    for order in permutations(range(7)):
        mat = block([weights[j] for j in order],entries)
        trace = mat[0]+mat[1]
        vertices[order] = (mat,trace)
        if trace<trace_cap:
            kept.add(order)
            lower.append(block_levels(mat)[0].lower())
        else:
            assert trace>trace_cap
    crossed = 0
    for order in kept:
        first,t1 = vertices[order]
        # 排列多面体的边交换两个相邻秩。每条截边只从被保留的一端计数。
        for j in range(6):
            pair = list(order)
            a,b = pair.index(j),pair.index(j+1)
            pair[a],pair[b] = pair[b],pair[a]
            k = tuple(pair)
            if k in kept:
                continue
            second,t2 = vertices[k]
            ratio = (trace_cap-t1)/(t2-t1)
            assert 0<ratio<1
            lower.append(block_levels(join(first,second,ratio))[0].lower())
            crossed += 1
    assert kept and crossed
    return min(lower),len(kept),crossed


def verify(root):
    mem = run_path(str(root/"ch39-next-shell-partial-coherence-up-threshold.py"))
    A,D,target,error = mem["source_data"](root)
    entries,_ = mem["target_data"](target)
    up,charm = arb(".00000615"),arb(".00311")
    def spectrum(r_text):
        r = arb(r_text)
        eig = ((1-r)*D+r*A).eig(algorithm="rump")
        assert all(v.imag.contains(0) and v.real>0 for v in eig)
        roots = sorted((v.real.sqrt() for v in eig),key=lambda w:w.mid())
        assert all(a<b for a,b in zip(roots,roots[1:])) and roots[0]>error
        lo = [(w-error).lower()**2 for w in roots[:7]]
        hi = [(w+error).upper()**2 for w in roots[:7]]
        assert all(a<b for a,b in zip(lo,lo[1:]))
        return lo,hi
    for r in (".86250",".862503"):
        lo,_ = spectrum(r)
        bound,old,new = clipped_minimum(lo,entries,up+charm)
        print("相干比例、保留原顶点、相交边、截面最小本征值下界",r,old,new,bound.str(14),flush=True)
        assert bound>up
    for r,t in ((".86251",".71008"),(".862505",".71005245")):
        _,hi = spectrum(r)
        first = block([hi[j] for j in (2,5,6,3,4,0,1)],entries)
        second = block([hi[j] for j in (2,6,5,3,4,0,1)],entries)
        light,heavy = block_levels(join(first,second,arb(t)))
        print("相干比例、固定边混合和两个本征值上界",r,t,light.upper().str(14),heavy.upper().str(14),flush=True)
        assert light<up and heavy<charm


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
