"""在固定能标与第四壳首枢轴交叉权下，用 K 与 D 界定任意荷六截止的中性规范残差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    n = 2_306_888_100_085
    native = lib("native-rh-ray-pvm")
    number = n+4
    _,_,pivots,_ = next(row for row in native["shell_weights"](root,number) if row[0]==6)
    c = [z.real for z in pivots[0]]
    assert sum(c).contains(1) and c[15]>0>c[17]
    _,ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    assert all(x is ordered[0] for x in ordered[:17])
    xmin,xmax = ordered[0],ordered[-1]
    _,heat,(fields,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=number)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec = 192
    p0,p1 = data["p0"],data["p1"]
    assert xmin>p0 and xmax<arb(".571")
    # 所有真实重排均满足这个逐层凸包。牛顿列以无质量斜率和全区间余项夹住。
    tv = sum(map(abs,c))
    beta = sum(a for _,a,_ in fields)/(2*arb.pi())
    assert beta<0
    newton_error = tv*xmax*sum(abs(a)*(1-heat(2,z/xmax)) for _,a,z in fields)/(2*arb.pi())
    neg_other = sum(-w for i,w in enumerate(c) if i!=17 and w<0)
    x18 = (((1-newton_error)/abs(beta)-neg_other*xmax)/abs(c[17])).lower()
    # 弱差为零时，第16层落在小截止区。其余十八层仍用全区间。
    low_weak = native["weak_series"](56)
    weak_band,large_band = lib("outer-shared-distance-kernel")["weak_bands"](low_weak,neutral["weak_column"],p0,p1)
    weak_upper = sum(w*(large_band if i==15 else weak_band) for i,w in enumerate(c)).upper()
    x16 = 256*p0
    assert newton_error<arb(".001737") and x18>arb(".04907203")
    assert x16.overlaps((arb("1459.0016")/arb("1.22089e19"))**2)
    assert weak_upper<arb("-.0152222")
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    species,_ = neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    log = (p1/p0).log()
    rule = gauge["angular_rule"](9,25)
    # 每个物种的有限动量角积分随截止平方单调。按系数符号合并端点包围。
    edges = [xmin,4*p0,16*p0,64*p0,x16,arb("1e-27"),arb("1e-20"),arb("1e-15"),p1,x18,xmax]
    assert all(l<r for l,r in zip(edges,edges[1:]))
    exposures = [[gauge["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*arb.pi())
                  for kind,z,_ in species] for x in edges]
    intervals = []
    start1,start2,start3 = data["initial"]
    starts = [5*start1/3+start2-8*start3/3,start2-start3]
    for j in range(2):
        coef = [(qq-arb(8)*colour/3 if j==0 else mixed-colour) for _,_,(qq,mixed,colour) in species]
        def band(a,b):
            return -sum(w*l.union(r) for w,l,r in zip(coef,exposures[a],exposures[b],strict=True))
        pieces = [band(i,i+1) for i in range(len(edges)-1)]
        def hull(parts):
            result = parts[0]
            for part in parts[1:]:
                result = result.union(part)
            return result
        all_x,low16,high18 = hull(pieces),hull(pieces[:4]),pieces[-1]
        res = starts[j]+sum(w*(low16 if i==15 else high18 if i==17 else all_x) for i,w in enumerate(c))
        intervals.append(res)
        print("规范列、全区间、第16层、第18层、首枢轴总残差",j+1,[(b.lower().str(12),b.upper().str(12)) for b in (all_x,low16,high18,res)],flush=True)
    print("K总误差、x18下界、x16阈值、第16层大截止时弱差上界",
          *[b.str(13) for b in (newton_error,x18,x16,weak_upper)],flush=True)
    assert arb("429.896")<intervals[0]<arb("701.904")
    assert arb("63.805")<intervals[1]<arb("137.580")
    return c,xmin,xmax,x18,p0,species,starts


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
