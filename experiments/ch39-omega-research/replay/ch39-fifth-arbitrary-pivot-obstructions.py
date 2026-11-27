"""原能标下用荷六源谱迹、牛顿真空两矩和角核，排除第五壳前两个枢轴的一切酉输运。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def inputs(root, *, number=2_306_888_100_090):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    n = 2_306_888_100_085
    _,source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    native = lib("native-rh-ray-pvm")
    _,_,pivots,_ = next(row for row in native["shell_weights"](root,number) if row[0]==6)
    weights = [[z.real for z in pivot] for pivot in pivots]
    thermal = lib("three-index-free-atoms")["broken_fields"](root,[number])
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=number)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec = 192
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    rows,_ = neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    return lib,source,weights,(thermal[0],thermal[1],thermal[2][0]),data,rows


def high_cutoff(weights,source,thermal,pivot):
    _,heat,field = thermal
    max_x = source[-1]
    tv = sum(map(abs,weights))
    pi = arb.pi()
    k = sum(a for _,a,_ in field)/(2*pi)
    w = -sum(g for g,_,_ in field)/(8*pi)
    ek = tv*max_x*sum(abs(a)*(1-heat(2,z/max_x)) for _,a,z in field)/(2*pi)
    ew = tv*max_x**2*sum(abs(g)*(arb(1)/2-heat(3,z/max_x)) for g,_,z in field)/(4*pi)
    assert k<0<w and weights[17]<0
    if pivot==1:
        other = sum(-v for i,v in enumerate(weights) if i!=17 and v<0)
        bound = (((1-ek)/abs(k)-other*max_x)/abs(weights[17])).lower()
        assert bound>arb(".04")
        return bound
    assert [i for i,c in enumerate(weights) if c>0]==[0,9,15]
    mplus = min(w.lower() for w in weights if w>0)
    mminus = max((-w).upper() for i,w in enumerate(weights) if i!=17 and w<0)
    trace = sum(x.sqrt() for x in source).upper()
    bound = arb("1e-5")
    delta = (1-ek)/abs(k)-abs(weights[17])*bound
    error = (arb("1e-122")+ew)/w
    assert mminus*trace**2<5 and delta>arb(".35")
    def exclude(pl,pr,ql,qr):
        if pr<mplus*ql**2 or (pl+delta)**3>(pr*qr**2+error)*(mminus*(trace-ql))**2:
            return
        assert max(pr-pl,qr-ql)>arb("1e-6")
        if pr-pl>qr-ql:
            mid = (pl+pr)/2
            exclude(pl,mid,ql,qr)
            exclude(mid,pr,ql,qr)
        else:
            mid = (ql+qr)/2
            exclude(pl,pr,ql,mid)
            exclude(pl,pr,mid,qr)
    exclude(arb(0),arb(5),arb(0),trace)
    return bound


def verify(root, *, number=2_306_888_100_090):
    lib,source,weights,thermal,data,rows = inputs(root,number=number)
    p0,p1 = data["p0"],data["p1"]
    shared = lib("fourth-variable-scale-neutral-gauge")
    gauge = lib("all-field-signed-crossing")
    rule = gauge["angular_rule"](9,25)
    series = lib("native-rh-ray-pvm")["weak_series"](56)
    neutral = lib("full-neutral-free-cutoff")
    a1,a2,a3 = data["initial"]
    starts = [5*a1/3+a2-arb(8)*a3/3,a2-a3]
    coeffs = [(qq-arb(8)*col/3,mix-col) for _,_,(qq,mix,col) in rows]
    photon_limits=[]
    for pivot,small,log_edge in ((1,(15,),8),(2,(0,15),13)):
        c = weights[pivot-1]
        bound = high_cutoff(c,source,thermal,pivot)
        limit = 2**log_edge*p0
        large = neutral["weak_column"](limit,p0,p1)
        whole = large
        for e in range(log_edge):
            step = arb(2)**(e-6)
            for j in range(64):
                whole = whole.union(series(((64+j)*step).union((65+j)*step)))
        d = sum(v*(large if i in small else whole) for i,v in enumerate(c))
        assert source[0]>p0 and d<0
        edges = sorted([source[0],limit,bound,source[-1]]+
                       [p0*2**e for e in range(4,101,5) if e!=log_edge and p0*2**e<source[-1]])
        columns = shared["angular_endpoint"](gauge,rule,rows,p0,p1,edges)
        result = []
        for j in range(2):
            bands = [-sum(b[j]*a.union(z) for b,a,z in zip(coeffs,left,right,strict=True))
                     for left,right in zip(columns,columns[1:])]
            all_x = shared["hull"](bands)
            low = shared["hull"](bands[:edges.index(limit)])
            high = shared["hull"](bands[edges.index(bound):])
            result.append([starts[j]+sum(v*(low if i==s else high if i==17 else all_x) for i,v in enumerate(c))
                           for s in small])
        print("枢轴、x18下界、无小截止时D上界、各分支RP与RM下界",pivot,
              bound.str(12),d.upper().str(12),[[x.lower().str(12) for x in result[j]] for j in range(2)],flush=True)
        assert all(v>0 for v in result[0])
        photon_limits.append([v.lower() for v in result[0]])
        if number==2_306_888_100_090:
            if pivot==1:
                assert bound>arb(".0485981") and d<arb("-.0151457")
                assert result[0][0]>arb("435.50") and result[1][0]>arb("68.85")
            else:
                assert d<arb("-.0057625")
                assert result[0][0]>arb("233.91") and result[0][1]>arb("3760.31")

    return photon_limits


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
