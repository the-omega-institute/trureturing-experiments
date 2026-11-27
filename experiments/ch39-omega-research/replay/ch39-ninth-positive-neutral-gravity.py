"""固定能标下，第九壳第三枢轴引力、弱差与两个中性规范残差的同谱正色解。"""
from inspect import getclosurevars
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


LOGS = (
    "-11.9594880255", "-2.9339", "-2.3336", "-26.3366",
    "-26.4466", "-2.5296", "-1.6744", "-2.7852",
    "-2.7651", "-39.2456", "-24.001", "-2.932",
    "-1.1765295239519", "-27.0467", "-1.5623280809274", "-36.923659144157",
    "-33.419307412191", "-2.6329832118945", "-20.8775",
)
ACTIVE = (0,12,14,15,16,17)


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number = 2_306_888_100_094
    _,sq = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=2_306_888_100_086)[2]
    floor,heavy,top = [sq[i].sqrt() for i in (0,17,18)]
    assert 0 < floor < heavy < top
    native = lib("native-rh-ray-pvm")
    _,_,piv,_ = next(row for row in native["shell_weights"](root,number) if row[0]==6)
    weights = [v.real for v in piv[2]]
    assert sum(weights).contains(1)
    gauge = lib("all-field-signed-crossing")
    data = gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses = lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    moments,heat,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    ctx.prec = 192
    mp = arb("1.22089e19")
    rho = (arb("246.22")/mp)**2/2
    neutral = lib("full-neutral-free-cutoff")
    species,_ = neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    p0,p1 = data["p0"],data["p1"]
    scales = (p1/p0).log()
    rule = gauge["angular_rule"](18,48)
    means = lib("fourth-third-pivot-cyclic-five")["angular_means"]
    slope,column = lib("fifth-seventh-two-pair-thermal-rotation")["thermal_columns"](field,moments,heat)
    finite = native["weak_series"](56)
    poly = getclosurevars(finite).nonlocals["coefficients"]
    ref = arb("-.074256577194663365545")
    terms = [(k,m,[(v/(4*arb.pi())) for v in (q-arb(8)*c/3,t-c,c)]) for k,m,(q,t,c) in species]
    a1,a2,a3 = data["initial"]
    starts = [-(17*floor+heavy+top),-arb(1),-arb("1e-122"),ref,5*a1/3+a2-8*a3/3,a2-a3,a3]
    logs = list(map(arb,LOGS))
    tiny = [z.exp()**2/p0 < arb("1e9") for z in logs]
    def at(i,z):
        a=z.exp();x=a*a
        weak = sum(c*(p0/x)**k for k,c in enumerate(poly,1)) if tiny[i] else arb(0)
        gr = [arb(0)]*3
        for k,m,c in terms:
            exp = gauge["angular_exposure"](k,m,x,p0,p1,scales,rule)
            for j in range(3):
                gr[j] -= c[j]*exp
        return [a,*[weights[i]*v for v in (*column(a),weak,*gr)]]
    def jac(zs):
        dirs=[]
        for i,z in zip(ACTIVE,zs,strict=True):
            a=z.exp();x=a*a
            dk,dw=slope(a)
            weak = -2*sum(k*c*(p0/x)**k for k,c in enumerate(poly,1)) if tiny[i] else arb(0)
            ang=[u-v for u,v in zip(means(p0/x),means(p1/x),strict=True)]
            g=[arb(0)]*3
            for k,m,c in terms:
                v=2*(-m/x).exp()*ang[k]
                for j in range(3):
                    g[j]-=c[j]*v
            dirs.append([a,*[weights[i]*v for v in (a*dk,a*dw,weak,*g)]])
        return arb_mat(dirs).transpose()
    def square(j):
        return arb_mat(j.tolist()[:6])
    fixed = list(starts)
    for i in range(19):
        if i not in ACTIVE:
            vals=at(i,logs[i])
            fixed=[a+b for a,b in zip(fixed,vals,strict=True)]
    def read(zs):
        row=list(fixed)
        for i,z in zip(ACTIVE,zs,strict=True):
            row=[a+b for a,b in zip(row,at(i,z),strict=True)]
        return row
    center=[logs[i] for i in ACTIVE]
    for _ in range(3):
        f=read(center)
        step=square(jac(center)).mid().solve(arb_mat([[v.mid()] for v in f[:6]])).entries()
        center=[(z-d).mid() for z,d in zip(center,step,strict=True)]
    f=read(center)
    pre=square(jac(center)).mid().inv().mid()
    assert not pre.det().contains(0)
    radius=arb(".000001")
    box=[z+arb(0,radius.upper()) for z in center]
    all_logs=list(logs)
    for i,z in zip(ACTIVE,box,strict=True):
        all_logs[i]=z
    noise=arb(0)
    for i,z in enumerate(all_logs):
        a=z.exp();x=a**2
        assert floor<a<top
        assert x>p0
        if tiny[i]:
            err=arb("1e-19")+2*x.upper()/p1
        else:
            err=(neutral["weak_column"](x.lower(),p0,p1)-ref).abs_upper()
        noise+=abs(weights[i]).upper()*err
    f[3]+=arb(0,noise.upper())
    jwide=jac(box)
    eye=arb_mat([[int(i==j) for j in ACTIVE] for i in ACTIVE])
    delta=eye-pre*square(jwide)
    disp=(pre*arb_mat([[v] for v in f[:6]])).entries()
    ratio=max(abs(d).upper()/radius+sum(v.abs_upper() for v in r) for d,r in zip(disp,delta.tolist(),strict=True))
    color=f[6]+sum(jwide[6,k]*(z-c) for k,(z,c) in enumerate(zip(box,center,strict=True)))
    print("六坐标中心",*[v.mid().str(18) for v in center],flush=True)
    print("弱核余量、自映比、端点色逆耦合",*[v.str(15) for v in (noise,ratio,color)],flush=True)
    assert noise<arb("2.67e-12") and ratio<arb(".000263")
    assert arb("32.01")<color<arb("32.03")
    return all_logs,color


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
