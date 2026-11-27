"""按 cb7c616e 的公开登记检验 N+5，包含原训练盒内完全相干的第四步。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number = 2_306_888_100_090
    z,dz,sectors,incoming,(low,h,k) = lib("fourth-third-pivot-joint")["verify"](root,return_transport=True)
    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    sites = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert sites==(14,8,7,6)
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    readings = []
    for iq,(charge,power,op) in enumerate(((-4,4,adj(z**4)),(2,2,z**2))):
        src,_,train,_,_ = sectors[iq]
        dst = spec(shell,charge,number)[1]
        ov = adj(dst)*op*src
        weights = [sum(abs(ov[i,j])**2*train[j] for j in range(19)) for i in range(19)]
        four = acb_mat([[dst[i,j] for j in range(19)] for i in sites])
        raw = lib("common-unitary-quark-right-trace")["weighted"](four,weights)
        err = (power*dz*(max(v.upper() for v in train)-min(v.lower() for v in train))).upper()
        levels = sorted((v.real for v in raw.eig(algorithm="rump")),key=lambda v:v.mid())
        levels = [v+arb(0,err) for v in levels]
        assert all(a<b for a,b in zip(levels,levels[1:]))
        status = []
        for j,(left,right) in enumerate(windows[iq]):
            lo,hi = arb(left),arb(right)
            if levels[j]>hi or levels[j+1]<lo:
                status.append("排除")
            elif levels[j]<lo and levels[j+1]>hi:
                status.append("覆盖全窗")
            else:
                status.append("未分离")
        published = [
            ((".000111260",".000111264"),(".007898666",".007898668"),(".405367789",".405367793"),(".911712786",".911712790")),
            ((".0002221057",".0002221061"),(".0005228074",".0005228078"),(".002650397",".002650399"),(".0891766721",".0891766725"))
        ]
        assert all(arb(lo)<v<arb(hi) for v,(lo,hi) in zip(levels,published[iq],strict=True))
        assert status==["排除","排除","覆盖全窗"]
        print("荷、四站本征值、三个单窗",charge,*[v.str(13) for v in levels],status,flush=True)
        readings.append((levels,status))

    native = lib("native-rh-ray-pvm")
    _,nodes,pivots,_ = next(row for row in native["shell_weights"](root,number) if row[0]==6)
    c = [w.real for w in pivots[2]]
    assert sum(c).contains(1)
    spectrum,dst = spec(shell,6,number)
    assert all(a.overlaps(b) for a,b in zip(nodes,spectrum,strict=True))
    ov = adj(dst)*z**6*incoming
    delta = (6*dz*(k-low)).upper()
    bands = []
    for i in range(19):
        amp = low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2
        bands.append((amp+arb(0,delta)).intersection(low.union(k))**2)
    assert all(x>0 for x in bands)
    moments,heat,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=number)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec = 192
    p0,p1 = data["p0"],data["p1"]
    fields,_ = neutral["response_rows"](data,masses,signs,(arb("246.22")/arb("1.22089e19"))**2/2,arb(".42995401872"))
    rule = gauge["angular_rule"](11,30)
    a1,a2,a3 = data["initial"]
    r = [arb(0),arb(0),arb(0),5*a1/3+a2-arb(8)*a3/3,a2-a3]
    log = (p1/p0).log()
    for x,co in zip(bands,c,strict=True):
        ka,wa = moments(field,heat,x)
        assert x>256*p0
        weak = neutral["weak_column"](x,p0,p1)
        rp,rm = arb(0),arb(0)
        for kind,mass,(qq,mixed,colour) in fields:
            ea = gauge["angular_exposure"](kind,mass,x.lower(),p0,p1,log,rule)
            eb = gauge["angular_exposure"](kind,mass,x.upper(),p0,p1,log,rule)
            ev = ea.union(eb)/(4*arb.pi())
            rp -= (qq-arb(8)*colour/3)*ev
            rm -= (mixed-colour)*ev
        for j,v in enumerate((ka,wa,weak,rp,rm)):
            r[j] += co*v
    goals = (arb(1),arb("1e-122"),arb(0),arb(0),arb(0))
    result = ["排除" if not a.overlaps(b) else "未分离" for a,b in zip(r,goals,strict=True)]
    assert result==["排除"]*5 and min(v.lower() for v in bands)>arb("3.36198e-6") and dz<arb("1.279e-10")
    bounds = (("-116.822847","-116.822845"),("40.069334","40.069336"),
              ("-.07425657719467","-.07425657719465"),("5.04","5.16"),("2.73","2.87"))
    assert all(arb(lo)<v<arb(hi) for v,(lo,hi) in zip(r,bounds,strict=True))
    print("截止平方最小下界、生成元算子误差",min(b.lower() for b in bands).str(13),dz.str(13),flush=True)
    print("荷六五项及逐项目标判定",*[b.str(13) for b in r],result,flush=True)
    return readings,r


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
