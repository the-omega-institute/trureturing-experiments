"""六壳固定补空间谱的矩松弛及全部可行对角的有理见证。"""
from itertools import product
import json
from pathlib import Path
from runpy import run_path

from flint import arb


def maximal_coordinate_cells(weights, second, fourth, largest_eigenvalue):
    # 任取一条最大占据线 k，v <= p_k² u。向外舍入地覆盖 p_k 的候选区间。
    if second.upper() <= 0 or fourth.upper() <= 0:
        return
    lower = (max(arb(0), fourth.lower()) / second.upper()).sqrt().lower()
    for k, weight in enumerate(weights):
        others = weights[:k] + weights[k+1:]
        holder2 = sum(1/x for x in others)**-1
        holder4 = sum(x**(-arb(1)/3) for x in others)**-3
        heaviest = max(x.upper() for x in others)
        total = sum(others)
        for cell in range(256):
            left = max(arb(cell)/256, lower)
            right = min(arb(cell+1)/256, largest_eigenvalue.upper())
            if left > right:
                continue
            if weight*left**2 + holder2*(1-right)**2 > second.upper():
                continue
            if weight*left**4 + holder4*(1-right)**4 > fourth.upper():
                continue
            rem2 = (second.upper() - weight*left**2).upper()
            rem4 = (fourth.upper() - weight*left**4).upper()
            if rem2 < 0 or rem4 < 0:
                continue
            mass = (1-left).upper()
            if weight*right**4 + min(heaviest*mass**4, mass**2*rem2) < fourth.lower():
                continue
            if weight*right**2 + min(heaviest*mass**2, (total*rem4).sqrt()) < second.lower():
                continue
            yield k, left, right


def certify(root, levels, blocks, model, records):
    sieve = run_path(str(root/"ch39-six-cyclic-line-sieve.py"))
    seven = run_path(str(root/"ch39-seven-cyclic-weak-horizon.py"))
    candidates = seven["CANDIDATES"][:6]
    zero_trace = sieve["squared"](blocks[0])
    comp_gram = sum((sieve["adj"](block)*block for block in blocks[1:]),
                    sieve["adj"](blocks[1])*blocks[1]*0)
    comp_trace = comp_gram[0,0].real + comp_gram[1,1].real
    trace = sum(sieve["squared"](block) for block in blocks[1:])
    largest = (1+(1-4*comp_gram.det().real/comp_trace**2).sqrt())/2
    assert arb(".99828") < largest < arb(".99830")
    rows = {tuple(r["lines"]): r["diagonals"] for r in records}
    assert len(rows) == len(records) == 15
    survivors = set()
    margin = arb(100)
    max_coordinate = arb(0)
    for ids in product(*candidates):
        targets = sieve["bounds"](model, levels, zero_trace, trace, ids)
        intervals = []
        for (u, v), shell in zip(targets, model, strict=True):
            weights = [-weight for k,weight in enumerate(shell["w"]) if k not in shell["triv"]]
            assert all(w>0 for w in weights)
            # 沿用六壳两条谱线筛选所用的完整误差预算。
            u += arb(0, (arb(".04")/(shell["ck"]*trace**2)).upper())
            v += arb(0, (arb(".045")/(shell["cw"]*trace**4)).upper())
            if max(arb(0), u.lower())**2 > sum(weights)*v.upper():
                break
            if any(bound.upper() < sum(w**(-arb(1)/(q-1)) for w in weights)**(1-q)
                   for q,bound in ((2,u),(4,v))):
                break
            if next(maximal_coordinate_cells(weights, u, v, largest), None) is None:
                break
            intervals.append((weights, u, v))
        else:
            key = tuple(k+1 for k in ids)
            survivors.add(key)
            for integers,(weights,u,v) in zip(rows[key], intervals, strict=True):
                assert len(integers)==12 and all(type(n) is int and n>=0 for n in integers)
                assert sum(integers)==10**9
                coords = [arb(n)/10**9 for n in integers]
                assert max(coords) < largest
                sums = [sum(w*p**q for w,p in zip(weights,coords,strict=True)) for q in (2,4)]
                # 外包矩形的端点就是本松弛的目标；两矩的证书放在各自矩形内部。
                gap = min((sums[0]-u.lower()).lower(), (u.upper()-sums[0]).lower(),
                          (sums[1]-v.lower()).lower(), (v.upper()-sums[1]).lower())
                margin = min(margin,gap)
                max_coordinate = max(max_coordinate, max(p.upper() for p in coords))
    assert survivors == set(rows) and margin > arb(".00525") and max_coordinate < arb(".850239")
    print("可行的低线组合、六壳逐壳见证数、最小矩区间余量、最大占据",
          len(survivors), 6*len(survivors), margin.str(12), max_coordinate.str(12),flush=True)


def experiment(root, fifth, columns, levels, records):
    sieve = run_path(str(root/"ch39-six-cyclic-line-sieve.py"))
    seven = run_path(str(root/"ch39-seven-cyclic-weak-horizon.py"))
    blocks, (_, _, model) = sieve["prepare_shells"](root,fifth,columns,levels,seven["CANDIDATES"][:6],arb(0))
    fields = run_path(str(root/"ch39-three-index-free-atoms.py"))["broken_fields"](
        root, range(2_306_888_100_092,2_306_888_100_098))[2]
    for shell,field in zip(model,fields,strict=True):
        ek,ew = sieve["thermal_defect"](field,levels[-1])
        budget = -sum(w for i,w in enumerate(shell["w"]) if i not in shell["triv"])
        assert ek*budget < arb(".04") and ew*budget < arb(".045")
    certify(root,levels,blocks,model,records)


if __name__=="__main__":
    root = Path(__file__).resolve().parent
    records = json.loads((root/"ch39-six-trivial-fixed-gram.json").read_text())
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"),records)
