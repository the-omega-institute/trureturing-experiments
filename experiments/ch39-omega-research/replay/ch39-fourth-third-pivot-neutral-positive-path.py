"""第四壳第三枢轴的三个端点条件与整段正中性、正色动能。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb, arb_mat, ctx


ANCHORS = ('1.583410934279855e-33', '2.171540424676002e-24', '8.115519225921293e-22', '0.04176967028436505', '0.04834721511054314', '0.04274205267028878', '2.372137868434049e-13', '3.154201720421295e-24', '1.363269289925297e-24', '4.945395480321781e-26', '4.811272484241265e-25', '0.1165733635593383', '0.04522493825058871', '4.180456380381771e-25', '9.075074283126069e-35', '7.859573448899221e-35', '9.074996611485406e-35', '9.074350332281186e-35', '2.221187606209506e-30')
ACTIVE = (15,0,11,17)


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
    sites=list(map(arb,ANCHORS))
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
    residual+=[sum(x.sqrt() for x in sites)-sum(source)]
    radius=arb("6e-7")
    # |a_n| <= 10*4^n/(n*n!)，对数坐标三阶导数的绝对和 <= 10*(16+4)*exp(4)。
    third_limit=200*arb(4).exp()
    step=arb("1e-8")
    avg=lib("fourth-third-pivot-cyclic-five")["angular_means"]
    # 弱差以原56阶多项式为光滑替身，真实读数和替身的统一误差留在自映预算内。
    def jac(xwide):
        mat=[]
        for i in ACTIVE:
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
            mat.append([c[i]*v for v in (*dg,dw)]+[x.sqrt()/2])
        return arb_mat(mat).transpose()
    pre=jac(arb(0)).mid().inv().mid()
    assert not pre.det().contains(0)
    d=arb_mat([[int(i==j) for j in range(4)] for i in range(4)])-pre*jac(radius)
    # 比较真实弱差和光滑替身在中心及盒内的变化，两端各计一次解析余项。
    noise=arb(0)
    for i in ACTIVE:
        x=sites[i]*arb(0,radius.upper()).exp()
        if sites[i]<1e5*p0:
            noise+=abs(c[i])*(arb("1e-19")+2*x/p1)
        else:
            noise+=abs(c[i])*abs(weak_at(x)-reference)
    residual[2]+=arb(0,(2*noise).upper())
    f=pre*arb_mat([[v] for v in residual])
    ratio=max((abs(v)/radius+sum(map(abs,row))).upper() for v,row in zip(f.entries(),d.tolist(),strict=True))
    print("端点与迹的四维自映比",ratio.str(12),flush=True)
    print("位移、余项",[(abs(v)/radius).str(7) for v in f.entries()],[sum(map(abs,row)).str(7) for row in d.tolist()],flush=True)
    assert ratio<1
    result=[v*arb(0,radius.upper()).exp() if i in ACTIVE else v for i,v in enumerate(sites)]
    gaps=lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,[v.sqrt() for v in result])
    # 对有符号权按逆截止分部求和；同号质量插入的平滑校正关于 1/x 的 Lipschitz 界与原逆矩系数相同。
    w,z,t=[q*rho for q in (arb(".63811")**2/2,(arb(".63811")**2+3*arb(".467453")**2/5)/2,arb(".8616")**2)]
    invcoef=(2*w+z)/(24*arb.pi())+(w+z)/(48*arb.pi())+t/(16*arb.pi())+t/(24*arb.pi())
    order=sorted(range(19),key=lambda i:sites[i].mid(),reverse=True)
    last=arb(0); tot=arb(1); variation=arb(0)
    for i in order:
        here=1/result[i]
        assert here>last
        variation+=(here-last)*abs(tot)
        tot-=c[i]
        last=here
    weak_bound=lib("full-neutral-running")["weak_bound"](data,[p0],[arb(0)],lib("neutral-axial-top"),lib("neutral-bosonic"))+invcoef*variation
    print("主序最小余量、整段弱差绝对上界",min(gaps).str(12),weak_bound.upper().str(12),flush=True)
    moments,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[number])
    thermal=[moments(field,heat,x) for x in result]
    kval,wval=[sum(weight*col[i] for weight,col in zip(c,thermal,strict=True)) for i in (0,1)]
    print("本解盒的 K、W",kval.str(12),wval.str(12),flush=True)
    return data,species,c,sites,result,weak_bound,lib


def path(data,species,weights,sites,box,weak,lib):
    gauge=lib("all-field-signed-crossing")
    p0,p1=data["p0"],data["p1"]
    a1,a2,a3=data["initial"]
    nodes,quadrature=gauge["angular_rule"](12,33)
    log=(p1/p0).log()
    eps=arb(2)**-52
    subdivisions=256
    values=[]
    for power in range(-205,6):
        base=arb(2)**power
        for j in range(subdivisions):
            point=base*(subdivisions+j)/subdivisions
            exact=point.expint(1)
            v=float(exact.mid())
            assert abs(exact-v)<arb("1e-12")
            values.append(v)
    assert (-arb(64)).exp()/64<arb("1e-12")
    values.append(0.)
    table=np.array(values)
    # 表格参数是 dyadic；线性插值二阶误差及 IEEE 四则运算的统一绝对误差。
    e_err=1/(8*arb(subdivisions)**2)+arb("3e-11")
    def E(y):
        assert np.all(np.isfinite(y)) and np.all(y>2.**-205)
        y=np.minimum(y,64.)
        m,k=np.frexp(y)
        pos=(2*m-1)*subdivisions
        j=pos.astype(np.int64)
        f=pos-j
        base=(k.astype(np.int64)-1+205)*subdivisions+j
        return table[base]+f*(table[np.minimum(base+1,len(table)-1)]-table[base])
    axes=3
    C=[arb(0)]*axes
    discarded=[arb(0)]*axes
    variation=[arb(0)]*axes
    coorderr=[arb(0)]*axes
    mass,angle,coeff=[],[],[]
    for c,x,xb in zip(weights,sites,box,strict=True):
        for kind,z,(qq,mm,cc) in species:
            bs=[-c*b/(4*arb.pi()) for b in (mm,qq-mm,cc)]
            for i,b in enumerate(bs):
                C[i]+=abs(b)
                variation[i]+=abs(b)*abs((xb/x).log()).upper()
            if z/x>100:
                for i,b in enumerate(bs):
                    discarded[i]+=abs(b)*log*(-arb(100)).exp()
                continue
            v=z/x
            vz=float(v.mid())
            assert v.is_zero() or abs(vz/v-1)<arb('1e-13')
            for h,omega in nodes:
                ph=h/x
                vh=float(ph.mid())
                assert abs(vh/ph-1)<arb('1e-13')
                angle.append(vh)
                mass.append(vz)
                point=[]
                for i,b in enumerate(bs):
                    b=b*omega[kind]
                    bf=float(b.mid())
                    coorderr[i]+=abs(b-bf)
                    point.append(bf)
                coeff.append(point)
    angle=np.array(angle);mass=np.array(mass)
    beta=np.array(coeff).T.copy()
    def pfloat(p):
        v=float(p.mid())
        assert abs(v/p-1)<arb("1e-13")
        return v
    zero=E(mass+pfloat(p0)*angle)
    intervals=[]
    for i,norm in enumerate(C):
        # 保守取逐项顺序求和的 IEEE 误差，即使 NumPy 采用了分块求和也成立。
        n=len(angle)+4
        assert n*eps<arb("1e-6")
        num=(n*eps/(1-n*eps))*(norm+coorderr[i])*(log+1)
        err=norm*quadrature*log+discarded[i]+variation[i]+coorderr[i]*(log+1)+2*(norm+coorderr[i])*e_err+num
        intervals.append(arb(0,err.upper()))
    initial=[a2,5*a1/3,a3]
    def end(j,d):
        if j==0:
            return initial
        p=p0*(log*arb(j)/2**d).exp()
        e=zero-E(mass+pfloat(p)*angle)
        cols=np.sum(beta*e,axis=1,dtype=np.float64)
        return [v+arb(float(k))+err for v,k,err in zip(initial,cols,intervals,strict=True)]
    # 对 t=log p 的二阶导数为 r(1-r)e^-y-r(ph/x)e^-y，r=ph/(z+ph)，绝对值 <=1/e。
    curv=[v/arb(1).exp() for v in C]
    stack=[(0,0,end(0,0),end(1,0))]
    leaves=0;floors=[arb(10)**7]*5
    while stack:
        j,d,left,right=stack.pop()
        chord=(log/2**d)**2/8
        l=[min(a.lower(),b.lower())-c*chord for a,b,c in zip(left,right,curv,strict=True)]
        mix,diff,col=l
        det=mix*diff-(mix+diff)*weak
        if mix>weak and diff>weak and col>2 and det>3:
            eigen=1/(1/mix+1/diff)-weak
            floors=[min(a,b.lower()) for a,b in zip(floors,[mix,diff,col,det,eigen],strict=True)]
            leaves+=1
        else:
            assert d<15,(j,d,l,det)
            mid=end(2*j+1,d+1)
            stack.extend(((2*j,d+1,left,mid),(2*j+1,d+1,mid,right)))
    print("连续能标闭段数、混合、QQ减混合、色、行列式、调和均值减弱差下界",
          leaves,*[v.str(12) for v in floors],flush=True)
    return leaves,floors


if __name__=="__main__":
    path(*verify(Path(__file__).resolve().parent))
