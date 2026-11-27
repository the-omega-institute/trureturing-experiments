"""冻结的六壳十八式见证在第十三壳原生第三枢轴上的读数。"""
from pathlib import Path
from runpy import run_path
from flint import arb,ctx

def readout(root,levels,weighted):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_098
    geo=lib("three-ray-chiral-boundary")
    native=lib("native-rh-ray-pvm")
    adj=lib("common-three-shell-nine-endpoints")["adj"]
    gram=native["gram_cholesky"](root)
    _, basis=lib("adjacent-pvm-transport")["spectrum"](geo,6,number)
    weights=[z.real for z in native["readout"](basis,range(19),[14,8,7,6],gram[3],3)]
    _,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[number])
    floor=levels[0]; rows=adj(basis)*weighted
    gravity=lib("shortest-common-complex-gate")["gravity_system"]
    values,_=gravity(rows,[],(floor,floor+1,floor+1),weights,heat,field,True)
    K=values[0,0]+1;W=values[1,0]+arb("1e-122")
    mp,mz=arb("1.22089e19"),arb("91.1876")
    p0,p1=(mz/mp)**2,(arb("1.176397e14")/mp)**2
    ref=arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    large=lib("full-neutral-free-cutoff")["weak_column"]
    minimum=min((floor+sum(abs(z)**2 for z in row)).lower() for row in rows.tolist())
    assert minimum**2>p0
    D=ref+arb(0,((large(minimum**2,p0,p1)-ref).abs_upper()*sum(abs(w).upper() for w in weights)).upper())
    print("第十三壳盒上 K、W、D 及最小归一化截止下界",*[v.str(18) for v in (K,W,D,minimum)],flush=True)
    return K,W,D

if __name__=="__main__":
    root=Path(__file__).resolve().parent
    src=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    a,b,levels=src["verify"](root,return_operator="spectral")
    witness=run_path(str(root/"ch39-six-unrestricted-eighteen.py"))
    readout(root,levels,witness["experiment"](root,a,b,levels))
