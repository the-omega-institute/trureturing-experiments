"""对已预登记的第十四壳，逐个读取七壳八站见证的冻结参数盒。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    number = 2_306_888_100_099
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    adj = lib("eight-optimal-classification")["adj"]
    _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
    weights = [z.real for z in native["readout"](
        basis, range(19), [14, 8, 7, 6], native["gram_cholesky"](root)[3], 3
    )]
    rows = adj(basis)
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [number])
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    small = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    floor = levels[0]
    mp, mz = arb("1.22089e19"), arb("91.1876")
    p0, p1 = (mz/mp)**2, (arb("1.176397e14")/mp)**2
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    records = {}

    def read(sites, image, low_rows, inverse, frame):
        amplitudes = rows*image
        # 第十四壳第一线属于平凡字符。先作固定系数行乘法，再作用于七个小行，
        # 避免先算十九个大振幅再相减而引入的区间依赖。
        trivial = rows*frame
        assert sum(abs(z)**2 for z in trivial.tolist()[0]) > arb(".5")
        low = (trivial*inverse)*low_rows
        for h in range(2):
            assert amplitudes[0,h].overlaps(low[0,h])
            amplitudes[0,h]=low[0,h]
        cuts = []
        for row in amplitudes.tolist():
            lo = floor.lower()+sum(max(arb(0),abs(z).lower())**2 for z in row)
            hi = floor.upper()+sum(abs(z).upper()**2 for z in row)
            cuts.append(lo.union(hi))
        assert all(a.is_finite() and a.lower()**2>p0 for a in cuts)
        values, _ = gravity(amplitudes, [], (floor, floor+1, floor+1), weights, heat, field, True)
        newton, vacuum = values[0,0]+1, values[1,0]+arb("1e-122")
        # 带质量低能级数和高能余项界按同一个冻结截止盒逐线使用。
        correction = arb(0)
        for w, a in zip(weights, cuts, strict=True):
            if a.upper()*mp < arb(100000):
                correction += w*(small((a*mp/mz)**2)-ref)
            else:
                correction += w*(large(a.lower()**2, p0, p1)-ref)
        weak = ref+correction
        assert all(v.is_finite() for v in (newton,vacuum,weak))
        assert sites not in records
        records[sites] = (newton,vacuum,weak)
        print(sites, *(v.str(14) for v in records[sites]), flush=True)

    lib("eight-optimal-classification")["experiment"](root, fifth, columns, levels, observe_projected=read)
    assert len(records)==32
    assert all(arb("3.3408")<k<arb("5.4011") and arb("-2.2166")<w<arb("-.9363")
               and arb("-.074265")<d<arb("-.074249") for k,w,d in records.values())
    for j,name in enumerate(("K", "W", "D")):
        column = [record[j] for record in records.values()]
        print("第十四壳", name, "全盒范围",
              min(v.lower() for v in column).str(14), max(v.upper() for v in column).str(14), flush=True)
    print("D 严格负、含零、严格正的盒数",
          sum(row[2]<0 for row in records.values()), sum(row[2].contains(0) for row in records.values()),
          sum(row[2]>0 for row in records.values()), flush=True)
    return records


if __name__ == "__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
