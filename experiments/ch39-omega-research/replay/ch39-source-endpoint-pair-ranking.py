"""在共同源标尺与源分布上，比较所有谱线两端单位份额转移的三残差响应体积。"""
from itertools import combinations
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def source_columns():
    repair = run_path(str(ROOT / "ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges, kernel, moments, _ = repair["source_data"]()
    _, lib, (rows,) = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT, [repair["SOURCE"]], tight_neutral=True)
    model, _, _, _, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0), *load("signed-gravity-unification")["gauge_weights"](model)]
    beta = model["beta"]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    direction = [1-t,t,-1]
    projection = arb_mat([[1,0,0,0,0], [0,1,0,0,0], [0,0,*direction]])
    coeff = projection * kernel["series"](rows,fields)
    scale = arb(7)/1000
    assert max(z.upper() for _,_,z in rows) < min(arb(1),scale*min(a.lower() for a,_,_ in edges))
    variation = sum(abs(w) for _,_,w in edges)
    factor = 2/(4*arb.pi()*factorial(kernel["ORDER"]+1))
    size = [sum(abs(2*a) for _,a,_ in rows), sum(abs(g) for g,_,_ in rows),
            sum(abs(h*b) for weights in fields for h,b in zip(direction,weights,strict=True))]
    ds = []
    for k in range(3):
        row = arb_mat([[coeff[k,j] for j in range(coeff.ncols())]])
        deriv = repair["dilate"](row,arb(1),True)*moments
        ds.append(deriv[0,0]+arb(0,(2*factor*size[k]*variation).upper()))
    columns = []
    for a,b,w in edges:
        endpoint_change = w*(kernel["basis"](scale*b)-kernel["basis"](scale*a))
        values = (coeff*endpoint_change).entries()
        columns.append([v+arb(0,(2*factor*z*abs(w)).upper())
                        for v,z in zip(values,size,strict=True)])
    return ds, columns


def scores(ds,columns):
    values = {}
    for i,j in combinations(range(len(columns)),2):
        x,y = columns[i],columns[j]
        det = (ds[0]*(x[1]*y[2]-x[2]*y[1])
               +ds[1]*(x[2]*y[0]-x[0]*y[2])
               +ds[2]*(x[0]*y[1]-x[1]*y[0]))
        values[i,j] = abs(det)
    return values


def rank(pair, values):
    subject = values[pair]
    rest = [v for key,v in values.items() if key!=pair]
    return 1+sum(v>subject for v in rest), 1+sum(not (v<subject) for v in rest)


def uniform_bound(columns):
    """外积坐标把指定线对的行列式写成三个竞争线对的固定线性组合。"""
    basis = arb_mat([columns[i] for i in (0,18,300)]).transpose()
    assert not basis.det().contains(0)
    coefficients = [basis.solve(arb_mat([columns[i]]).transpose()).entries() for i in (16,301)]
    x,y = coefficients
    minors = [x[i]*y[j]-x[j]*y[i] for i,j in combinations(range(3),2)]
    bound = sum(abs(z) for z in minors)
    print("三个竞争线对的系数", [a.str(12) for a in minors])
    print("对任意标尺响应的统一上界系数",bound.str(12))
    assert bound < arb(".024107")
    return bound


def verify():
    ds, columns = source_columns()
    uniform_bound(columns)
    values = scores(ds,columns)
    pair = (16,301)
    by_midpoint = sorted(values,key=lambda p: float(values[p].mid()),reverse=True)
    assert len(values) == 64980
    assert by_midpoint[0] == (300,323) and rank((300,323),values) == (1,1)
    assert rank(pair,values) == (112,112)
    four = {p: values[p] for p in combinations((16,285,301,323), 2)}
    print("原四线也改用端点后的六对排序", sorted(four, key=lambda p: -float(four[p].mid())))
    assert rank(pair,four) == (4,4)
    print("16、301 在六对端点方向中的名次", rank(pair,four))
    print("首名与指定线对的响应体积比", (values[300,323]/values[pair]).str(12))
    for candidate in [*by_midpoint[:5],pair]:
        print("线对、名次范围、响应体积", candidate, rank(candidate,values),
              values[candidate].str(12))


if __name__ == "__main__":
    verify()
