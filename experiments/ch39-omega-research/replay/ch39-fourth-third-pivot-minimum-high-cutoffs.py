"""第四壳第三枢轴的五式端点解至少需要三个高截止，且三个可以达到。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


ANCHORS = (
    "8.327631195715618e-11", "1.296765884551616e-17", "6.039630895538626e-25",
    ".228018575663224", "5.43492569418729e-17", "4.505216600638419e-28",
    "6.918875631155653e-15", "3e-7", ".3545731159040746", "3.307777910935902e-34",
    "3e-7", "3e-7", "8.716261913347304e-11", "3e-7", "3e-7", "6.152688440550019e-24",
    ".0119591223613301", "2.631223560752675e-22", "2.022618420309069e-33",
)


def exclude_pairs(c,source,heat,field,cut):
    pi=arb.pi()
    xmax=source[-1]**2
    total=sum(source)
    slope_k=sum(a for _,a,_ in field)/(2*pi)
    slope_w=-sum(g for g,_,_ in field)/(8*pi)
    assert slope_k<0<slope_w
    def remainder(x):
        return (x*sum(abs(a)*(1-heat(2,z/x)) for _,a,z in field)/(2*pi),
                x*x*sum(abs(g)*(arb(1)/2-heat(3,z/x)) for g,_,z in field)/(4*pi))
    # x^d(E_{d+1}(0)-E_{d+1}(z/x)) 对 x 递增，由积分中 1-(1+v)e^{-v}>=0 得出。
    high_k,high_w=remainder(xmax)
    low_k,low_w=remainder(cut)
    counts=0
    for i,j in combinations(range(19),2):
        high=abs(c[i])+abs(c[j])
        low=sum(abs(w) for k,w in enumerate(c) if k not in (i,j))
        error_k=high*high_k+low*(abs(slope_k)*cut+low_k)
        error_w=high*high_w+low*(abs(slope_w)*cut**2+low_w)
        sk=arb(1)+arb(0,error_k.upper())
        sw=arb("1e-122")+arb(0,error_w.upper())
        trace=(total-17*cut.sqrt()).union(total)
        # 两个保留位置已覆盖真实高位，其余十七个振幅每个至多 sqrt(cut)。
        # 用保迹限制第二个振幅，两个热矩直接在第一个振幅的一维区间上排除。
        start=(trace.lower()-source[-1].upper()).lower()
        assert start>0
        stack=[(start,source[-1].upper())]
        while stack:
            lo,hi=stack.pop()
            counts+=1
            a=lo.union(hi)
            b=(trace-a).intersection(start.union(source[-1]))
            def power(v,n):
                return (v.lower()**n).union(v.upper()**n)
            b4=(sw/slope_w-c[i]*power(a,4))/c[j]
            if not b4.overlaps(power(b,4)):
                continue
            b=b4.intersection(power(b,4)).sqrt().sqrt()
            k=slope_k*(c[i]*power(a,2)+c[j]*power(b,2))
            w=slope_w*(c[i]*power(a,4)+c[j]*power(b,4))
            if not k.overlaps(sk) or not w.overlaps(sw):
                continue
            assert hi-lo>arb("1e-7"),(i,j,a,b,k,w)
            mid=(lo+hi)/2
            stack.extend(((lo,mid),(mid,hi)))
    print("阈值、全部层对数、振幅区间数",cut.str(8),171,counts,flush=True)


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_089
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    source=[v.sqrt() for v in ordered]
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,number) if row[0]==6)
    c=[v.real for v in pivots[2]]
    moments,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[number])
    cut=arb("1e-6")
    exclude_pairs(c,source,heat,field,cut)
    exclude_pairs(c,source,heat,field,arb("9e-5"))
    gauge=lib("all-field-signed-crossing")
    data=gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=192
    p0,p1=data["p0"],data["p1"]
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    neutral=lib("full-neutral-free-cutoff")
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    norm=4*arb.pi()
    rows=[(k,z,((q-arb(8)*col/3)/norm,(mix-col)/norm,col/norm)) for k,z,(q,mix,col) in species]
    a1,a2,a3=data["initial"]
    sites=list(map(arb,ANCHORS))
    selected=(17,0,3,8,16)
    parameter=9
    width=arb("3e-8")
    radius=arb("1e-5")
    rule=gauge["angular_rule"](19,54)
    means=lib("fourth-third-pivot-cyclic-five")["angular_means"]
    exact=[]
    for x in sites:
        values=[gauge["angular_exposure"](k,z,x,p0,p1,(p1/p0).log(),rule) for k,z,_ in rows]
        gauge_cols=[-sum(b[j]*v for (_,_,b),v in zip(rows,values,strict=True)) for j in range(3)]
        exact.append([*gauge_cols[:2],*moments(field,heat,x),gauge_cols[2]])
    def jac(wide):
        cols=[]
        for i in selected:
            x=sites[i]*arb(0,wide.upper()).exp()
            diff=[u-v for u,v in zip(means(p0/x),means(p1/x),strict=True)]
            g=[-sum(b[j]*(-z/x).exp()*diff[k] for k,z,b in rows) for j in range(2)]
            k=x*sum(a*(-z/x).exp() for _,a,z in field)/(2*arb.pi())
            w=-x*x*sum(s*(-z/x).exp() for s,_,z in field)/norm
            cols.append([c[i]*v for v in (*g,k,w)]+[x.sqrt()/2])
        return arb_mat(cols).transpose()
    constants=(5*a1/3+a2-8*a3/3,a2-a3,-1,-arb("1e-122"))
    residual=[v+sum(w*r[j] for w,r in zip(c,exact,strict=True)) for j,v in enumerate(constants)]
    residual.append(sum(x.sqrt() for x in sites)-sum(source))
    low=sites[parameter]*arb(0,width.upper()).exp()
    gauge_norms=[sum(abs(b[j]) for _,_,b in rows) for j in range(3)]
    dparam=[abs(c[parameter])*v for v in (*gauge_norms[:2],low*sum(abs(a) for _,a,_ in field)/(2*arb.pi()),
              low**2*sum(abs(s) for s,_,_ in field)/norm)]+[low.sqrt()/2]
    residual=arb_mat([[v+arb(0,(width*d).upper())] for v,d in zip(residual,dparam,strict=True)])
    pre=jac(arb(0)).mid().inv().mid()
    assert not pre.det().contains(0)
    identity=arb_mat([[int(i==j) for j in range(5)] for i in range(5)])
    defect=identity-pre*jac(radius)
    offset=max(abs(v).upper() for v in (pre*residual).entries())/radius
    slope=max(sum(map(abs,row)).upper() for row in defect.tolist())
    assert offset+slope<1
    boxes=[v*arb(0,radius.upper()).exp() if i in selected else v for i,v in enumerate(sites)]
    boxes[parameter]=low
    weak=native["weak_series"](56)
    def weak_at(x):
        return weak(x/p0) if x<1e7*p0 else neutral["weak_column"](x,p0,p1)
    # 如果第一、十、十六层全在 2.3 TeV 以上，只剩第四层正权能在低端贡献弱差。
    ref=arb("-.074256577194663365545")
    positive=[i for i,v in enumerate(c) if v>0]
    assert positive==[0,3,9,15] and all(c[i]<0 for i in range(19) if i not in positive)
    def weak_range(start):
        end=arb(10)**7
        lower=start
        bounds=neutral["weak_column"](end*p0,p0,p1)-ref
        while lower<end:
            upper=min((lower*arb("1.0008")).upper(),end)
            bounds=bounds.union(weak(lower.union(upper))-ref)
            lower=upper
        return bounds
    all_weak=weak_range((ordered[0]/p0).lower())
    large_weak=weak_range((arb(2300)/arb("91.1876"))**2)
    cap=ref*sum(c)+c[3]*all_weak.upper()+sum(c[i]*large_weak.upper() for i in (0,9,15))
    cap+=sum(c[i]*all_weak.lower() for i in range(19) if i not in positive)
    print("假设三个大正权截止全在2.3 TeV以上时的弱差上界",cap.upper().str(12),flush=True)
    assert cap<0
    rest=sum(c[i]*weak_at(x) for i,x in enumerate(boxes) if i!=parameter)
    ends=[rest+c[parameter]*weak_at(sites[parameter]*(s*width).exp()) for s in (-1,1)]
    assert ends[0]>0>ends[1]
    gaps=lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,[x.sqrt() for x in boxes])
    assert [i+1 for i,x in enumerate(boxes) if x>cut]==[4,9,17]
    assert all(x>cut or x<cut for x in boxes)
    colour=a3+sum(w*r[4] for w,r in zip(c,exact,strict=True))
    variation=width*abs(c[parameter])+radius*sum(abs(c[i]) for i in selected)
    colour+=arb(0,(variation*gauge_norms[2]).upper())
    assert colour>0
    print("三高解的自映位移、导数余量、两端弱差、主序最小余量、端点色逆耦合",
          *(v.str(10) for v in (offset,slope,*ends,min(gaps),colour)),flush=True)
    running=lib("full-neutral-running")
    for mu in ("1e4","1e9"):
        p=(arb(mu)/arb("1.22089e19"))**2
        variation=sum(abs(w)*abs((xb/x).log()).upper()*min(arb(1),((p-p0)/(4*xb.lower())).upper())
                      for w,x,xb in zip(c,sites,boxes,strict=True))
        vals=running["direct"](data,sites,c,species,mu,gauge,gauge["angular_rule"](13,37))
        colour=vals[2]+arb(0,(gauge_norms[2]*variation).upper())
        print("检验能标与色逆耦合上下界",mu,colour.lower().str(11),colour.upper().str(11),flush=True)
        if mu=="1e9":
            assert colour<0
        else:
            assert colour>0
    return data,species,c,sites,boxes,lib


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
