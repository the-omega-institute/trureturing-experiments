"""六壳低线选法中一组逐壳可达的固定谱矩矩形不能共用同一个正矩阵。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def dominant_lines(cells):
    allowed = {}
    for k, left, right in cells:
        old = allowed.get(k, (arb(1), arb(0)))
        allowed[k] = (min(left, old[0]), max(right, old[1]))
    return allowed


def propagate(anchor, interval, columns):
    left, right = interval
    result = []
    for overlap in (anchor*columns).entries():
        # a 是两条谱线的重叠平方；端点仅用乘法与非负平方根。
        a = abs(overlap)**2
        a_lo = max(arb(0), a.lower())
        a_hi = min(arb(1), a.upper())
        assert 0 <= left <= right <= 1 and 0 <= a_lo <= a_hi <= 1
        minus = (left*a_lo).sqrt() - ((1-left)*(1-a_lo)).sqrt()
        lower = max(arb(0), minus.lower())**2
        if right < a_lo:
            p, w = right, a_lo
        elif a_hi < left:
            p, w = left, a_hi
        else:
            result.append((lower.lower(), arb(1)))
            continue
        upper = ((p*w).sqrt()+((1-p)*(1-w)).sqrt())**2
        result.append((lower.lower(), upper.upper()))
    return result


def certify(root, levels, blocks, model, frames):
    sieve = run_path(str(root/"ch39-six-cyclic-line-sieve.py"))
    zero_trace = sieve["squared"](blocks[0])
    trace = sum(sieve["squared"](m) for m in blocks[1:])
    bounds = sieve["bounds"](model,levels,zero_trace,trace,(0,0,15,0,15,0))
    data = []
    all_lines = []
    for s,(u,v),f in zip(model,bounds,frames,strict=True):
        others = [i for i in range(19) if i not in s["triv"]]
        weights = [-s["w"][i] for i in others]
        assert all(w>0 for w in weights)
        u += arb(0,(arb(".04")/(s["ck"]*trace**2)).upper())
        v += arb(0,(arb(".045")/(s["cw"]*trace**4)).upper())
        data.append((weights,u,v))
        all_lines.extend([f[i] for i in others])
    lines = acb_mat(all_lines)
    columns = sieve["adj"](lines)
    gram = sieve["adj"](blocks[1])*blocks[1]+sieve["adj"](blocks[2])*blocks[2]
    total = gram[0,0].real + gram[1,1].real
    largest = (1+(1-4*gram.det().real/total**2).sqrt())/2
    coordinate_cells = run_path(str(root/"ch39-six-trivial-fixed-gram.py"))["maximal_coordinate_cells"]
    tables = []
    for shell in (0,5):
        anchors = dominant_lines(coordinate_cells(*data[shell],largest))
        tables.append({j:propagate(acb_mat([all_lines[12*shell+j]]),window,columns)
                      for j,window in anchors.items()})
    assert set(tables[0]) == {1, 2, 9} and set(tables[1]) == {2, 9}
    gap = arb(10)
    for anchor7,anchor12 in product(*tables):
        rows = list(zip(tables[0][anchor7],tables[1][anchor12],strict=True))
        low = [max(a[0],b[0]) for a,b in rows]
        high = [min(a[1],b[1]) for a,b in rows]
        violations = []
        for j,(weights,u,v) in enumerate(data):
            left, right = low[12*j:12*j+12], high[12*j:12*j+12]
            violations.extend((sum(left)-1,1-sum(right)))
            for q, target in ((2,u),(4,v)):
                lo = sum(w*x**q for w,x in zip(weights,left,strict=True))
                hi = sum(w*x**q for w,x in zip(weights,right,strict=True))
                violations.extend((lo-target.upper(),target.lower()-hi))
        best = max(val.lower() for val in violations)
        print("第七、十二壳的补空间谱行序号与最强排除余量",anchor7+1,anchor12+1,best.str(12),flush=True)
        assert best>arb(".12")
        gap = min(gap,best)
    print("六种最大线配对的最小排除余量",gap.str(12),flush=True)


def experiment(root,fifth,columns,levels,check=certify):
    sieve = run_path(str(root/"ch39-six-cyclic-line-sieve.py"))
    opts = run_path(str(root/"ch39-seven-cyclic-weak-horizon.py"))["CANDIDATES"][:6]
    blocks, (_,_,model) = sieve["prepare_shells"](root,fifth,columns,levels,opts,arb(0))
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    chars,_ = run_path(str(root/"ch39-cyclic-two-characters-nine.py"))["characters"](graph["geo"]["sites"])
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    frame = acb_mat([a+b for a,b in zip(chars[1].tolist(),chars[2].tolist(),strict=True)])
    numbers=range(2_306_888_100_092,2_306_888_100_098)
    fields=run_path(str(root/"ch39-three-index-free-atoms.py"))["broken_fields"](root,numbers)[2]
    frames=[]
    for number,field,s in zip(numbers,fields,model,strict=True):
        ek,ew=sieve["thermal_defect"](field,levels[-1])
        cost=-sum(w for i,w in enumerate(s["w"]) if i not in s["triv"])
        assert ek*cost<arb(".04") and ew*cost<arb(".045")
        _,basis=spectrum(graph,6,number)
        frames.append((sieve["adj"](basis)*frame).tolist())
    check(root,levels,blocks,model,frames)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
