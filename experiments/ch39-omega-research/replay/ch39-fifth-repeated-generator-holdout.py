"""三壳校准后，冻结第三步生成元并在 N+5 留出检验六个质量窗及荷六三项。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def mass_window(lower,upper,left,right):
    if upper<left or lower>right:
        return "排除"
    if lower<right and upper>left:
        return "可达"
    return "未判定"


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    z,err,flavor,incoming,(low,h,k) = lib("fourth-full-ckm-gravity")["verify"](root,return_transport=True)
    n = lib("rank-two-generator-up-quark")["N"]+5
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    keep = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert keep==(14,8,7,6)
    print("第五壳零点编号与冻结第三步的算子误差上界",n,err.upper().str(13),flush=True)
    quarks=[]
    for q,power,op,(src,weights),bounds in zip((-4,2),(4,2),(adj(z**4),z**2),flavor,windows,strict=True):
        _,dest=spectrum(shell,q,n)
        overlap=adj(dest)*op*src
        amplitudes=[sum(w*abs(overlap[i,j])**2 for j,w in enumerate(weights)) for i in range(19)]
        sites=acb_mat([[dest[i,j] for j in range(19)] for i in keep])
        block=acb_mat([[sites[i,j]*amplitudes[j] for j in range(19)] for i in range(4)])*adj(sites)
        vals=block.eig(algorithm="rump")
        assert all(v.imag.contains(0) for v in vals)
        spread=max(w.upper() for w in weights)-min(w.lower() for w in weights)
        wiggle=arb(0,(power*err*spread).upper())
        levels=sorted((v.real+wiggle for v in vals),key=lambda v:v.mid())
        assert all(a<b for a,b in zip(levels,levels[1:]))
        status=[mass_window(levels[j],levels[j+1],arb(l),arb(r)) for j,(l,r) in enumerate(bounds)]
        print("荷、四站有序谱、任意三平面三个质量窗判定",q,[v.str(15) for v in levels],status,flush=True)
        quarks.append((levels,status))

    native=lib("native-rh-ray-pvm")
    _,pivot_levels,cross,_=next(row for row in native["shell_weights"](root,n) if row[0]==6)
    moments,heat,fields=lib("three-index-free-atoms")["broken_fields"](root,[n])
    ctx.prec=384
    eig,dest=spectrum(shell,6,n)
    assert all(a.overlaps(b) for a,b in zip(eig,pivot_levels,strict=True))
    ov=adj(dest)*z**6*incoming
    wiggle=arb(0,(6*err*(k-low)).upper())
    amplitudes=[low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2+wiggle for i in range(19)]
    assert all(v>0 for v in amplitudes)
    p0,p1=[(arb(v)/arb("1.22089e19"))**2 for v in ("91.1876","1.176397e14")]
    weak=lib("full-neutral-free-cutoff")["weak_column"]
    thermal=[moments(fields[0],heat,f*f) for f in amplitudes]
    columns=(*zip(*thermal),[weak(f*f,p0,p1) for f in amplitudes])
    weights=[c.real for c in cross[0]]
    assert sum(weights).contains(1)
    gravity=[sum(c*a for c,a in zip(weights,col,strict=True)) for col in columns]
    status=["排除" if not v.overlaps(t) else "未判定" for v,t in zip(gravity,(arb(1),arb("1e-122"),arb(0)),strict=True)]
    print("荷六最小振幅下界、K、W、D 与逐项目标判定",min(v.lower() for v in amplitudes).str(15),
          *[v.str(16) for v in gravity],status,flush=True)
    return quarks,gravity,status


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
