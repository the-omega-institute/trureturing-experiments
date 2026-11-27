"""同一第三步以荷二谱线相位保持三个下型中心，并以三次幂精确校准第四壳荷六。"""
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


PHASE_LINES = (18,0,6,15,9,8)


def adj(a):
    return a.conjugate().transpose()


def one_rotation(seed,projection,eye):
    complement = eye-projection
    diag = projection*seed*projection+complement*seed*complement
    error = eye-adj(diag)*diag
    eps = sum(abs(v).upper()**2 for v in error.entries()).sqrt().upper()
    assert eps<arb(".001")
    inverse_root = acb_mat(eye)
    term = acb_mat(eye)
    for j in range(1,11):
        term *= error
        inverse_root += term*(arb(comb(2*j,j))/4**j)
    tail = arb(0,(eps**11/(1-eps)).upper())
    inverse_root += acb_mat([[acb(tail,tail) for j in range(19)] for i in range(19)])
    block = diag*inverse_root
    assert (adj(block)*block).contains(eye)
    off = complement*adj(block)*seed*projection
    u = adj(acb_mat([[off[3,j] for j in range(19)]]))
    u *= 1/(adj(u)*u)[0,0].real.sqrt()
    b = off*u
    b *= 1/(adj(b)*b)[0,0].real.sqrt()
    sn,cs = arb(".008").sin_cos()
    turn = eye+(u*adj(u)+b*adj(b))*(cs-1)+(b*adj(u)-u*adj(b))*sn
    base = block*turn
    assert (adj(base)*base).contains(eye) and sn**2<arb(".000064")
    print("七维空间的单条泄漏概率及修正酉种子的 Frobenius 距离",
          (sn**2).str(13),sum(abs(v).upper()**2 for v in (base-seed).entries()).sqrt().str(10),flush=True)
    return base


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    x,y,incoming,dest,(low,h,k),c,(moments,heat,field) = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    data = lib("common-embedding-central-right-ckm")
    _,labels,_ = data["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    source,t = [lib("adjacent-pvm-transport")["spectrum"](lib("three-ray-chiral-boundary"),2,n+i)[1] for i in (1,4)]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"](root,rank["eye"],"ch39-fourth-down-gravity-cayley.tsv")
    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    rotation = acb_mat([[int(i==cycle[j]) for j in range(19)] for i in range(19)])
    eye = rank["eye"](19)
    projection = (eye+rotation+rotation**2)*(arb(1)/3)
    base = one_rotation(seed,projection,eye)
    down = rank["read"](source,t,labels[n+1,2],base*y*x)
    four = acb_mat([[t[i,j] for j in range(19)] for i in (14,8,7,6)])
    raw = lib("common-unitary-quark-right-trace")["weighted"](four,down)
    vals = sorted((v.real for v in raw.eig(algorithm="rump")),key=lambda z:z.mid())
    mass = list(map(arb,(".0000135",".000268",".01401")))
    assert all(vals[j]<v<vals[j+1] for j,v in enumerate(mass))
    probabilities = [prod(v-m for m in mass)/prod(v-w for j,w in enumerate(vals) if j!=i)
                     for i,v in enumerate(vals)]
    assert all(p>0 for p in probabilities) and sum(probabilities).contains(1)
    print("四站下型本征值与法向最小概率",*[v.str(13) for v in vals],
          min(p.lower() for p in probabilities).str(10),flush=True)

    # Z=T D T*B，故 Tg*Z³I=C D V D V D H，其中 I 为两列荷六相干输入。
    read = dest*t
    v = adj(t)*base*t
    a = adj(t)*base*incoming
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    high_weak = lib("full-neutral-free-cutoff")["weak_column"]
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    def rough(s):
        return c[15]*weak(s)+(1-c[15])*arb("-.074256577194663365545")
    left,right = arb("1.28"),arb(2)
    assert rough(left)>0>rough(right)
    for _ in range(42):
        m = (left+right)/2
        if rough(m)>0:
            left=m
        else:
            assert rough(m)<0
            right=m
    left-=arb("1e-11")
    right+=arb("1e-11")
    time = left.union(right)

    def system(point,s):
        phases = [acb(1)]*19
        for line,z in zip(PHASE_LINES,point,strict=True):
            phases[line] = acb(0,z).exp()
        vd = acb_mat([[v[i,j]*phases[j] for j in range(19)] for i in range(19)])
        cd = acb_mat([[read[i,j]*phases[j] for j in range(19)] for i in range(19)])
        first = vd*a
        second = vd*first
        out = cd*second
        mid,last = cd*vd,cd*vd*vd
        do = [acb_mat([[acb(0,1)*(cd[i,l]*second[l,j]+mid[i,l]*first[l,j]+last[i,l]*a[l,j])
                        for j in (0,1)] for i in range(19)]) for l in PHASE_LINES]
        cutoff = mz*s.sqrt()/mp
        eta = ((cutoff-low)/(h-low)).sqrt()
        z = out[15,0]-eta
        values = [z.real,z.imag,out[15,1].real,out[15,1].imag]
        jac = [[getattr(d[15,j],part) for d in do] for j in (0,1) for part in ("real","imag")]
        f = [low+(h-low)*abs(out[i,0])**2+(k-low)*abs(out[i,1])**2 for i in range(19)]
        f[15] = cutoff
        errors = (abs(c[15])*cutoff**2*sum(abs(a) for _,a,_ in field)/(2*arb.pi()),
                  abs(c[15])*cutoff**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi()))
        thermal = [arb(0,e.upper())-target for e,target in zip(errors,(arb(1),arb("1e-122")),strict=True)]
        slopes = [[],[]]
        for i,u in enumerate(f):
            if i==15:
                for row in slopes:
                    row.append(arb(0))
                continue
            assert u>arb("1e-6")
            for j,z in enumerate(moments(field,heat,u*u)):
                thermal[j] += c[i]*z
            decay = [(-z/u**2).exp() for _,_,z in field]
            slopes[0].append(c[i]*u*sum(a*e for (_,a,_),e in zip(field,decay,strict=True))/arb.pi())
            slopes[1].append(-c[i]*u**3*sum(g*e for (g,_,_),e in zip(field,decay,strict=True))/(2*arb.pi()))
        df = arb_mat([[2*((h-low)*(out[i,0].conjugate()*d[i,0]).real+
                          (k-low)*(out[i,1].conjugate()*d[i,1]).real) for d in do] for i in range(19)])
        values.extend(thermal)
        jac.extend((arb_mat(slopes)*df).tolist())
        return arb_mat([[v] for v in values]),arb_mat(jac),f

    center = [arb(0)]*6
    for _ in range(3):
        val,jac,_ = system(center,time.mid())
        correction = jac.mid().solve(val.mid())
        center = [(v-z).mid() for v,z in zip(center,correction.entries(),strict=True)]
    width = arb("1e-10")
    box = [v+arb(0,width.upper()) for v in center]
    val,jac,_ = system(center,time)
    pre = jac.mid().inv()
    assert not pre.det().contains(0)
    _,wide,amps = system(box,time)
    defect = arb_mat([[int(i==j) for i in range(6)] for j in range(6)])-pre*wide
    q = max(sum(abs(v).upper() for v in row) for row in defect.tolist())
    offset = max(abs(v).upper() for v in (pre*val).entries())/width
    assert q+offset<1
    other = sum(c[i]*high_weak(a*a,p0,p1) for i,a in enumerate(amps) if i!=15)
    ends = [c[15]*weak(s)+other for s in (left,right)]
    assert ends[0]>0>ends[1]
    print("六相位中心最大模、位移占比、导数行和、两个弱差符号",
          *[v.str(12) for v in (max(map(abs,center)),offset,q,*ends)],flush=True)
    print("第16层截止 GeV、第18和19层振幅",*[v.str(13) for v in (mz*time.sqrt(),amps[17],amps[18])],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
