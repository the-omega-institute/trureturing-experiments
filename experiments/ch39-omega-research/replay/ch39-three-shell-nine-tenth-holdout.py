"""按预登记用第十壳的原生谱基读取同一九式闭盒内的荷六算子。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    source, levels = lib("common-three-shell-nine-endpoints")["verify"](root)
    shell, native = lib("three-ray-chiral-boundary"), lib("native-rh-ray-pvm")
    number = 2_306_888_100_095
    nodes,basis = lib("adjacent-pvm-transport")["spectrum"](shell,6,number)
    _,other,pivots,_ = next(r for r in native["shell_weights"](root,number) if r[0]==6)
    assert all(a.overlaps(b) for a,b in zip(nodes,other,strict=True))
    weights = [w.real for w in pivots[2]]
    moments,heat,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    ctx.prec = 256
    _,column = lib("fifth-seventh-two-pair-thermal-rotation")["thermal_columns"](field,moments,heat)
    p = basis.conjugate().transpose()*source
    low = levels[0]
    amps = [low+sum((levels[j+1]-low)*abs(p[i,j])**2 for j in range(2)) for i in range(19)]
    readings = [column(a) for a in amps]
    k,w = [sum(weight*row[j] for weight,row in zip(weights,readings,strict=True)) for j in range(2)]
    minimum = min(a.lower() for a in amps)
    p0 = (arb("91.1876")/arb("1.22089e19"))**2
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    assert minimum**2 > p0 and sum(weights).contains(1)
    ref = arb("-.074256577194663365545")
    remainder = lib("full-neutral-free-cutoff")["weak_column"](minimum**2,p0,p1)-ref
    d = ref + sum(weight*remainder for weight in weights)
    for label,value in zip(("牛顿项","真空项","弱差","最小振幅下界","第十六线振幅"),(k,w,d,minimum,amps[15]),strict=True):
        print("第十壳",label,value.lower().str(17),value.upper().str(17),flush=True)
    assert arb(".9756574") < k < arb(".9756590")
    assert arb(".00117051") < w < arb(".00117060")
    assert arb("-.07425657719466338") < d < arb("-.07425657719466335")
    assert minimum > arb("3.6663e-7") and amps[15] < arb("3.6669e-7")
    return k,w,d,minimum,amps[15]


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
