"""第四壳原生第三枢轴的稀疏低截止端点解与全程正动能。"""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, ctx

ROWS=(0, 2, 9, 13, 14, 15, 16, 17, 18)
ACTIVE=(7, 5, 0)
ANCHORS=('7.778373229889432e-35', '8.823322531355292e-25', '8.823318506239489e-25', '8.823399941841849e-25', '4.51250576265497e-34', '4.513467084206126e-34', '4.512764659148236e-34', '2.179594500793099e-34', '3.126525805585156e-30')

def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_089
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    source=[z.sqrt() for z in ordered]
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,number) if row[0]==6)
    c=[z.real for z in pivots[2]]
    assert sum(c).contains(1)
    gauge=lib("all-field-signed-crossing")
    data=gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=192
    p0,p1=data["p0"],data["p1"]
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    neutral=lib("full-neutral-free-cutoff")
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    norm=4*arb.pi()
    coeff=[[(arb(q)-arb(8)*col/3)/norm,(m-col)/norm] for _,_,(q,m,col) in species]
    a1,a2,a3=data["initial"]
    starts=(5*a1/3+a2-8*a3/3,a2-a3)
    low=list(map(arb,ANCHORS))
    active=[ROWS[i] for i in ACTIVE]
    def fill(v):
        return ((sum(source)-sum(x.sqrt() for x in v))/(19-len(v)))**2
    high=fill(low)
    sites=[high]*19
    for i,x in zip(ROWS,low,strict=True):
        sites[i]=x
    log=(p1/p0).log()
    rule=gauge["angular_rule"](19,54)
    weak=native["weak_series"](56)
    reference=arb("-.074256577194663365545")
    def weak_at(x):
        return weak(x/p0) if x<1e5*p0 else neutral["weak_column"](x,p0,p1)
    columns=[]
    for x in sites:
        vals=[gauge["angular_exposure"](k,z,x,p0,p1,log,rule) for k,z,_ in species]
        columns.append([-sum(b[j]*v for b,v in zip(coeff,vals,strict=True)) for j in range(2)]+[weak_at(x)])
    residual=[s+sum(w*r[j] for w,r in zip(c,columns,strict=True)) for j,s in enumerate((*starts,0))]
    radius=arb("6e-7")
    # |a_n| <= 10*4^n/(n*n!)，对数坐标三阶导数的绝对和 <= 10*(16+4)*exp(4)。
    third_limit=200*arb(4).exp()
    step=arb("1e-8")
    avg=lib("fourth-third-pivot-cyclic-five")["angular_means"]
    # 弱差以原56阶多项式为光滑替身，真实读数和替身的统一误差留在自映预算内。
    def jac(xwide):
        mat=[]
        for i in active:
            xcenter=sites[i]
            x=xcenter*arb(0,xwide.upper()).exp()
            diff=[u-v for u,v in zip(avg(p0/x),avg(p1/x),strict=True)]
            dg=[-sum(b[j]*(-z/x).exp()*diff[k] for b,(k,z,_) in zip(coeff,species,strict=True)) for j in range(2)]
            if xcenter<1e5*p0:
                upper=weak(xcenter*step.exp()/p0);lower=weak(xcenter*(-step).exp()/p0)
                hd=arb("1e-5")
                assert xcenter*(-hd-xwide-step).exp()>p0
                second=(weak(xcenter*hd.exp()/p0)-2*weak(xcenter/p0)+weak(xcenter*(-hd).exp()/p0))/hd**2
                slope_limit=abs(second)+third_limit*(hd+xwide+step)
                dw=(upper-lower)/(2*step)+arb(0,(slope_limit*(step+xwide)).upper())
            else:
                dw=arb(0)
            mat.append([c[i]*v for v in (*dg,dw)])
        return arb_mat(mat).transpose()
    pre=jac(arb(0)).mid().inv().mid()
    assert not pre.det().contains(0)
    d=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])-pre*jac(radius)
    # 比较真实弱差和光滑替身在中心及盒内的变化，两端各计一次解析余项。
    noise=arb(0)
    for i in active:
        x=sites[i]*arb(0,radius.upper()).exp()
        if sites[i]<1e5*p0:
            noise+=abs(c[i])*(arb("1e-19")+2*x/p1)
        else:
            noise+=abs(c[i])*abs(weak_at(x)-reference)
    residual[2]+=arb(0,(2*noise).upper())
    low_box=[x*arb(0,radius.upper()).exp() if i in ACTIVE else x for i,x in enumerate(low)]
    high_box=fill(low_box)
    high_weight=sum(c[i] for i in range(19) if i not in ROWS)
    high_delta=abs((high_box/high).log())
    for j in range(2):
        residual[j]+=arb(0,(abs(high_weight)*sum(abs(b[j]) for b in coeff)*(p1-p0)*high_delta/(4*high_box.lower())).upper())
    residual[2]+=arb(0,(2*abs(high_weight)*(arb("2.398e-21")+arb("7.468e-36")/high_box.lower())).upper())
    f=pre*arb_mat([[v] for v in residual])
    ratio=max((abs(v)/radius+sum(map(abs,row))).upper() for v,row in zip(f.entries(),d.tolist(),strict=True))
    print("三项端点条件的自映比",ratio.str(12),flush=True)
    print("位移、余项",[(abs(v)/radius).str(7) for v in f.entries()],[sum(map(abs,row)).str(7) for row in d.tolist()],flush=True)
    assert ratio<1
    result=[high_box]*19
    for i,x in zip(ROWS,low_box,strict=True):
        result[i]=x
    assert all(x<arb("1e-6") for x in low_box) and high_box>arb("1e-6")
    gaps=lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,[v.sqrt() for v in result])
    # 对有符号权按逆截止分部求和；同号质量插入的平滑校正关于 1/x 的 Lipschitz 界与原逆矩系数相同。
    w,z,t=[q*rho for q in (arb(".63811")**2/2,(arb(".63811")**2+3*arb(".467453")**2/5)/2,arb(".8616")**2)]
    invcoef=(2*w+z)/(24*arb.pi())+(w+z)/(48*arb.pi())+t/(16*arb.pi())+t/(24*arb.pi())
    centres=[*low,high]
    boxes=[*low_box,high_box]
    weights=[*[c[i] for i in ROWS],high_weight]
    lib("neutral-positive-profile")["verify"]()
    ctx.prec=192
    order=sorted(range(len(centres)),key=lambda i:centres[i].mid(),reverse=True)
    last=arb(0); tot=arb(1); variation=arb(0)
    for i in order:
        here=1/boxes[i]
        assert here>last
        variation+=(here-last)*max(arb(0),(-tot).upper())
        tot-=weights[i]
        last=here
    weak_bound=lib("full-neutral-running")["weak_bound"](data,[p0],[arb(0)],lib("neutral-axial-top"),lib("neutral-bosonic"))+invcoef*variation
    print("主序最小余量、整段负弱差上界",min(gaps).str(12),weak_bound.upper().str(12),flush=True)
    return data,species,weights,centres,boxes,weak_bound,lib


def certify_path(root):
    data,species,weights,centres,boxes,weak,lib=verify(root)
    path=lib("fourth-third-pivot-neutral-positive-path")["path"]
    middle=arb("1e4")
    pm=(middle/arb("1.22089e19"))**2
    first=dict(data,p1=pm)
    path(first,species,weights,centres,boxes,weak,lib)
    # 第二段以精确匹配的高能端点为弱差参考。小截止直接用两端 E1 界，
    # 大截止用无截止端点矩和有限截止的逆矩余项。
    running=lib("full-neutral-running")
    axial,boson=lib("neutral-axial-top"),lib("neutral-bosonic")
    endpoint=running["weak_bound"]({"p0":pm},[pm],[arb(0)],axial,boson)
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    w,z,t=[q*rho for q in (arb(".63811")**2/2,(arb(".63811")**2+3*arb(".467453")**2/5)/2,arb(".8616")**2)]
    coeff=(2*w+z)/(24*arb.pi())+(w+z)/(48*arb.pi())+t/(16*arb.pi())+t/(24*arb.pi())
    small=[i for i,x in enumerate(boxes) if 4*x.upper()/pm<coeff/x.upper()]
    large=[i for i in range(len(boxes)) if i not in small]
    weak_high=abs(sum(weights[i] for i in large))*endpoint
    weak_high+=sum(4*abs(weights[i])*boxes[i].upper()/pm for i in small)
    weak_high+=sum(abs(weights[i])*coeff/boxes[i].lower() for i in large)
    gauge=lib("all-field-signed-crossing")
    projected=[(k,z,(m,q-m,col)) for k,z,(q,m,col) in species]
    a1,a2,a3=data["initial"]
    initial=[a2,5*a1/3,a3]
    rule=gauge["angular_rule"](16,45)
    for weight,x,box in zip(weights,centres,boxes,strict=True):
        dx=abs((box/x).log()).upper()
        for kind,z,beta in projected:
            h=gauge["angular_exposure"](kind,z,x,data["p0"],pm,(pm/data["p0"]).log(),rule)
            h+=arb(0,(dx*min(arb(1),((pm-data["p0"])/(4*box.lower())).upper())).upper())
            for j in range(3):
                initial[j]-=weight*beta[j]*h/(4*arb.pi())
    second=dict(data,p0=pm,initial=(3*initial[1].lower()/5,initial[0].lower(),initial[2].lower()))
    print("后段弱差绝对上界",weak_high.upper().str(12),flush=True)
    return path(second,species,weights,centres,boxes,weak_high,lib)


if __name__=="__main__":
    certify_path(Path(__file__).resolve().parent)
