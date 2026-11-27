"""从第三枢轴五式的十九个振幅构造实 Schur–Horn 双列，用两次反射延拓并检验其主立方根。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def inn(a,b):
    return (adj(a)*b)[0,0]


def col(a,j):
    return acb_mat([[a[i,j]] for i in range(a.nrows())])


def horn_columns(amps,low,h,k):
    small = h-low
    large = k-low
    d = [v-low for v in amps]
    assert all(v>0 for v in d)
    i,j = 11,7
    others = [l for l in range(19) if l not in (i,j)]
    assert all(d[l]<d[i] for l in others) and d[i]<d[j]
    rest = small-sum(d[l] for l in others)
    assert 0<rest<d[i]<large
    sn = ((d[i]-rest)/(large-rest)).sqrt()
    cs = ((large-d[i])/(large-rest)).sqrt()
    rem = (rest/small).sqrt()
    us,vs = [arb(0)]*19,[arb(0)]*19
    for l in others:
        us[l] = (d[l]/small).sqrt()
    us[i],us[j],vs[i],vs[j] = rem*cs,-rem*sn,sn,cs
    columns = acb_mat([[a,b] for a,b in zip(us,vs,strict=True)])
    assert (adj(columns)*columns).contains(acb_mat([[1,0],[0,1]]))
    assert all((low+small*a*a+large*b*b).overlaps(x) for x,a,b in zip(amps,us,vs,strict=True))
    return columns,rest,sn


def reflect(s,t):
    a = inn(s,t)
    assert abs(a)>0
    direction = (s+t*(a.conjugate()/abs(a)))*(1/(2+2*abs(a)).sqrt())
    return direction


def root_rule(incoming,target,amps,levels):
    eye = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    unit = acb_mat([[target[i,j]*target[0,j].conjugate()/abs(target[0,j]) for j in range(19)] for i in range(19)])
    assert min(abs(target[0,j]).lower() for j in range(19))>0
    columns,rest,sn = horn_columns(amps,*levels)
    physical = unit*columns
    first = reflect(col(incoming,0),col(physical,0))
    h1 = eye-first*adj(first)*2
    second = reflect(h1*col(incoming,1),col(physical,1))
    h2 = eye-second*adj(second)*2
    b = inn(first,second)
    r = abs(b)
    assert arb(".037")<r<arb(".039")
    e = (second*(b.conjugate()/r)-first*r)*(1/(1-r*r).sqrt())
    theta = 2*r.acos()
    assert arb("1.021")<theta/3<arb("1.023")
    sr,cr = (theta/3).sin_cos()
    z = eye+(first*adj(first)+e*adj(e))*(cr-1)+(e*adj(first)-first*adj(e))*sr
    assert (adj(z)*z).contains(eye)
    assert all(v.contains(0) for v in (z**3-h2*h1).entries())
    print("剩余较低高源迹、Schur–Horn 正弦、反射法向内积模、主根角",
          *[s.str(12) for s in (rest,sn,r,theta/3)],flush=True)
    return z


def original_plane(lib):
    emb = lib("common-embedding-central-right-ckm")
    _,_,intrinsic,rotation = emb["common_embedding"](*lib("independent-unitary-quark-right-holdout")["GEOMETRY"])
    quat = lib("fourth-full-ckm-gravity")["QUATERNION"]
    norm = arb(sum(v*v for v in quat)).sqrt()
    qx,qy,qz,qw = [arb(v)/norm for v in quat]
    plane_turn = acb_mat([[acb(qw,-qz),acb(-qy,-qx),0,0],[acb(qy,-qx),acb(qw,qz),0,0],
                         [0,0,1,0],[0,0,0,1]])
    s,t,u = [1/arb(n).sqrt() for n in (3,2,6)]
    four_basis = acb_mat([[1,0,0,0],[0,s,t,u],[0,s,-t,u],[0,s,0,-2*u]])
    return four_basis*plane_turn*rotation*intrinsic


def midpoint(m):
    return acb_mat([[acb(v.real.mid(),v.imag.mid()) for v in row] for row in m.tolist()])


def frobenius(m):
    return sum(abs(v).upper()**2 for v in m.entries()).sqrt().upper()


def ckm_on_plane(lib,raws,errors):
    frame = original_plane(lib)
    assert (adj(frame)*frame).contains(acb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
    eig = lib("full-ckm-generator-next-shell")["sorted_eigen"]
    pairs = [eig(adj(frame)*raw*frame) for raw in raws]
    angles = []
    for (levels,_),err in zip(pairs,errors,strict=True):
        gaps = [min(abs(m-w).lower() for j,w in enumerate(levels) if i!=j) for i,m in enumerate(levels)]
        assert min(gaps)>2*err
        angles.append([arb(2).sqrt()*err/(gap-err) for gap in gaps])
    up,down = (p[1] for p in pairs)
    row = []
    for i,j in ((0,1),(1,2),(0,2),(1,1)):
        error = arb(0,(angles[0][i]+angles[1][j]).upper())
        row.append(inn(up[i],down[j])+acb(error,error))
    signed = (row[0]*row[1]*row[2].conjugate()*row[3].conjugate()).imag
    vals = [abs(v) for v in row[:3]]
    print("原固定三平面 |Vus|、|Vcb|、|Vub|、带号 J",
          *[v.str(11) for v in (*vals,signed)],flush=True)
    assert arb(".68")<vals[0]<arb(".72") and vals[1]<arb(".053") and vals[2]<arb(".047")
    assert abs(signed)<arb(".001")
    return vals,signed


def quark_reading(lib,z,x,y,sectors):
    centre = midpoint(z)
    dz = frobenius(z-centre)
    raws,errors = [],[]
    for (src,dst,train),(charge,power,prop) in zip(sectors,((-4,2,adj(y*y)*adj(x*x)),(2,1,y*x)),strict=True):
        op = adj(centre*centre) if charge==-4 else centre
        weights = lib("rank-two-generator-up-quark")["read"](src,dst,train,op*prop)
        mat = acb_mat([[dst[i,j] for j in range(19)] for i in (14,8,7,6)])
        raw = lib("common-unitary-quark-right-trace")["weighted"](mat,weights)
        mid = midpoint(raw)
        mid = (mid+adj(mid))*(arb(1)/2)
        eta = (1+dz)**power-1
        err = (frobenius(raw-mid)+max(w.upper() for w in train)*eta*(2+eta)).upper()
        vals = sorted((v.real for v in mid.eig(algorithm="rump")),key=lambda v:v.mid())
        vals = [v+arb(0,err) for v in vals]
        assert all(a<b for a,b in zip(vals,vals[1:]))
        bounds = [(".040",".042"),(".241",".243"),(".362",".365"),(".411",".414")] if charge==-4 else [
            (".0050",".0052"),(".0188",".0191"),(".0225",".0228"),(".0925",".0928")]
        assert all(arb(l)<v<arb(r) for (l,r),v in zip(bounds,vals,strict=True))
        print("荷、算子误差、四站本征值",charge,err.str(8),*[v.str(12) for v in vals],flush=True)
        raws.append(mid)
        errors.append(err)
    return ckm_on_plane(lib,raws,errors)


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    amps = lib("fourth-third-pivot-five")["verify"](root,refined_amplitudes=True)
    x,y,incoming,readout,levels,_,_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    ctx.prec = 384
    z = root_rule(incoming,adj(readout),amps,levels)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    sectors = []
    for charge in (-4,2):
        src,dst = [spec(shell,charge,n+j)[1] for j in (1,4)]
        sectors.append((src,dst,labels[n+1,charge]))
    return quark_reading(lib,z,x,y,sectors)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
