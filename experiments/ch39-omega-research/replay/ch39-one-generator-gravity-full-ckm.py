"""同一酉发生元的荷六精确校准、六个 Yukawa 中心与完整带号 CKM。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


DIRECTIONS = ((14,9,"i"), (14,8,"r"), (14,12,"r"), (7,6,"i"), (14,7,"r"),
              (14,14,"i"), (14,6,"i"), (18,8,"r"), (18,6,"r"), (8,0,"i"),
              (6,3,"i"), (15,11,"i"), (18,15,"r"), (5,0,"r"), (7,7,"i"),
              (8,1,"r"), (10,7,"i"), (14,7,"i"))


def adj(a):
    return a.conjugate().transpose()


def targets(root, independent, embedding):
    raw = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    diag = run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))["diag"]
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))["standard_matrix"]
    projection = run_path(str(root/"ch39-three-gram-source-flavor.py"))["projections"]
    mix = acb_mat(ckm(*[arb(t).sin() for t in (".22704",".04275",".00377")],arb("1.139")))
    masses = [list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    a,b = [diag([y/m[-1] for y in m]) for m in masses]
    a = adj(mix)*a*mix
    basis = raw["target_basis"](a,b*arb(".1"))*diag([acb(0,1),1,acb(0,-1)])
    u,v,intrinsic,rotation = embedding["common_embedding"](*independent["GEOMETRY"])
    frame = rotation*intrinsic
    assert (adj(frame)*frame).contains(diag([1]*3))
    coords = independent["coordinates"]
    proj,goals = [],[]
    for mass,mat in zip(masses,(a,b),strict=True):
        desired = adj(basis)*mat*basis*mass[-1]
        block = rotation*embedding["extend_target"](desired,u,v)*adj(rotation)
        compressed = adj(frame)*block*frame
        assert (compressed-desired).contains(acb_mat(3,3))
        levels,p = projection(compressed*compressed)
        assert all(w.overlaps(y*y) for w,y in zip(levels,mass,strict=True))
        proj.append(p)
        goals.append(arb_mat([[z] for z in coords(block)]))
    pu,pd = proj
    assert all((a*b).trace().real.overlaps(abs(mix[i,j])**2)
               for i,a in enumerate(pu) for j,b in enumerate(pd))
    signed = (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    expected = (mix[0,0]*mix[1,1]*mix[0,1].conjugate()*mix[1,0].conjugate()).imag
    assert signed>0 and signed.overlaps(expected)
    print("目标标架的九个模平方及带号 J",[[ (a*b).trace().real.str(10) for b in pd] for a in pu],
          signed.str(13),flush=True)
    return goals


def verify(root):
    get = lambda file: run_path(str(root/file))
    rank = get("ch39-rank-two-generator-up-quark.py")
    n = rank["N"]+2
    emb = get("ch39-common-embedding-central-right-ckm.py")
    independent = get("ch39-independent-unitary-quark-right-holdout.py")
    spec = get("ch39-adjacent-pvm-transport.py")["spectrum"]
    shell = get("ch39-three-ray-chiral-boundary.py")
    _,labels,_ = emb["verify_nodes"](root,((n-1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ordered = get("ch39-adjacent-ordered-charge-six-cutoffs.py")["verify"](root)[2][1]
    assert all(z is ordered[0] for z in ordered[:17])
    low,h,k = [ordered[j].sqrt() for j in (0,17,18)]
    moments,heat,fields = get("ch39-three-index-free-atoms.py")["broken_fields"](root,[n])
    field = fields[0]
    native = get("ch39-native-rh-ray-pvm.py")
    _,_,pivots,_ = next(row for row in native["shell_weights"](root,n) if row[0]==6)
    c = [z.real for z in pivots[0]]
    assert sum(c).contains(1)
    weak = native["weak_series"]()
    weak_high = get("ch39-full-neutral-free-cutoff.py")["weak_column"]
    ctx.prec = 384
    spectra = [[spec(shell,charge,node)[1] for node in (n-1,n)] for charge in (-4,2,6)]
    goals = targets(root,independent,emb)
    s,t,u = [1/arb(j).sqrt() for j in (3,2,6)]
    four_basis = acb_mat([[1,0,0,0],[0,s,t,u],[0,s,-t,u],[0,s,0,-2*u]])
    sectors = []
    for charge,(source,target) in zip((-4,2),spectra[:2],strict=True):
        site = adj(four_basis)*acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
        columns = [acb_mat([[site[i,j]] for i in range(4)]) for j in range(19)]
        compressed_projectors = [v*adj(v) for v in columns]
        for projector in compressed_projectors:
            assert all(projector[i,j].contains(0) for i in (0,1) for j in (2,3))
            assert (projector[2,2]-projector[3,3]).contains(0) and projector[2,3].real.contains(0)
        coords = arb_mat([independent["coordinates"](p) for p in compressed_projectors]).transpose()
        sectors.append((source,adj(target),labels[n-1,charge],coords))
    source,target = spectra[2]
    source = acb_mat([[source[i,17],source[i,18]] for i in range(19)])
    dest = adj(target)
    eye = rank["eye"](19)
    x0 = get("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    basis = []
    for i,j,axis in DIRECTIONS:
        e = acb_mat(19,19)
        z = acb(0,1) if axis=="i" else acb(1)
        e[i,j],e[j,i] = z,-z.conjugate()
        basis.append(e)
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    near,far = arb("1.30"),arb("1.50")
    def rough(t):
        return c[15]*weak(t)+(1-c[15])*arb("-.074256577194663365545")
    assert rough(near)>0>rough(far)
    for _ in range(32):
        mid = (near+far)/2
        if rough(mid)>0:
            near = mid
        else:
            assert rough(mid)<0
            far = mid
    # 第16层的定向重叠方程用任意一个精确选定的单位相位。
    start = (dest*x0*x0*x0*source)[15,0]
    phase = start/abs(start)

    def system(point,time):
        cutoff = (p0*time).sqrt()
        eta = ((cutoff-low)/(h-low)).sqrt()
        H = sum((e*z for e,z in zip(basis,point,strict=True)),acb_mat(19,19))
        inv = (eye-H).inv()
        x = (2*inv-eye)*x0
        dx = [2*inv*e*inv*x0 for e in basis]
        x2 = x*x
        dx2 = [d*x+x*d for d in dx]
        values, rows = [],[]
        for (src,dst,weights,coords),goal,op,dop in zip(
                sectors,goals,(adj(x2),x),([adj(d) for d in dx2],dx),strict=True):
            ov = dst*op*src
            wr = arb_mat([[sum(w*abs(ov[i,j])**2 for j,w in enumerate(weights))] for i in range(19)])
            values.extend((coords*wr-goal).entries())
            jac = []
            for d in dop:
                do = dst*d*src
                jac.append((coords*arb_mat([[2*sum(w*(ov[i,j].conjugate()*do[i,j]).real for j,w in enumerate(weights))]
                                            for i in range(19)])).entries())
            rows.extend(arb_mat(jac).transpose().tolist())
        ov = dest*(x2*x)*source
        dov = [dest*(d2*x+x2*d)*source for d2,d in zip(dx2,dx,strict=True)]
        residual = ov[15,0]-eta*phase
        values.extend((residual.real,residual.imag,ov[15,1].real,ov[15,1].imag))
        rows.extend([[getattr(d[15,j],part) for d in dov] for j in (0,1) for part in ("real","imag")])
        amplitudes = [low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2 for i in range(19)]
        amplitudes[15] = cutoff
        error = (abs(c[15])*cutoff**2*sum(abs(a) for _,a,_ in field)/(2*arb.pi()),
                 abs(c[15])*cutoff**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi()))
        thermal = [arb(0,b.upper())-goal for b,goal in zip(error,(arb(1),arb("1e-122")),strict=True)]
        slopes = [[],[]]
        for i,y in enumerate(amplitudes):
            if i==15:
                for col in slopes:
                    col.append(arb(0))
                continue
            assert y>arb(".00001")
            terms = moments(field,heat,y*y)
            for j,val in enumerate(terms):
                thermal[j] += c[i]*val
            decay = [(-z/y**2).exp() for _,_,z in field]
            slopes[0].append(c[i]*y*sum(a*b for (_,a,_),b in zip(field,decay,strict=True))/arb.pi())
            slopes[1].append(-c[i]*y**3*sum(g*b for (g,_,_),b in zip(field,decay,strict=True))/(2*arb.pi()))
        da = arb_mat([[2*((h-low)*(ov[i,0].conjugate()*d[i,0]).real+
                          (k-low)*(ov[i,1].conjugate()*d[i,1]).real) for d in dov] for i in range(19)])
        values.extend(thermal)
        rows.extend((arb_mat(slopes)*da).tolist())
        return arb_mat([[v] for v in values]),arb_mat(rows),amplitudes

    time = near.union(far)
    center = [arb(0)]*18
    for _ in range(2):
        val,jac,_ = system(center,time.mid())
        step = jac.mid().solve(val.mid())
        center = [(z-d).mid() for z,d in zip(center,step.entries(),strict=True)]
    radius = arb("1e-9")
    box = [z+arb(0,radius.upper()) for z in center]
    val,jac,_ = system(center,time)
    pre = jac.mid().inv()
    assert not pre.det().contains(0)
    _,wide,amplitudes = system(box,time)
    defect = arb_mat([[int(i==j) for j in range(18)] for i in range(18)])-pre*wide
    jac_error = max(sum(abs(y) for y in row).upper() for row in defect.tolist())
    offset = max(abs(y).upper() for y in (pre*val).entries())/radius
    print("十八参数最大中心、半宽、导数行和及中心位移占比",max(abs(y).upper() for y in center).str(12),
          radius.str(5),jac_error.str(10),offset.str(10),flush=True)
    assert jac_error+offset<1
    other = sum(c[j]*weak_high(y*y,p0,p1) for j,y in enumerate(amplitudes) if j!=15)
    ends = [c[15]*weak(t)+other for t in (near,far)]
    assert ends[0]>0>ends[1]
    energy = mz*time.sqrt()
    print("荷六第16层能量及弱差端点",energy.str(16),*[z.str(10) for z in ends],flush=True)
    return center,radius


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
