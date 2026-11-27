"""N+3 单个新生成元以 Y^-2 和 Y 联合延续两荷完整四站块。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


DIRECTIONS = ((3,0,"i"),(11,3,"i"),(18,5,"r"),(9,8,"r"),(16,7,"r"),(17,14,"i"),
              (17,6,"r"),(16,6,"i"),(8,4,"r"),(13,7,"r"),(8,4,"i"),(14,8,"i"))


def adj(a):
    return a.conjugate().transpose()


def outer_sum(a):
    return max(sum(abs(v).upper() for v in row) for row in a.tolist())


def local(seed,point,directions,eye,derivative=True):
    h = acb_mat(19,19)
    for z,(i,j,part) in zip(point,directions,strict=True):
        b = acb(0,z) if part=="i" else acb(z)
        h[i,j] += b
        if i!=j:
            h[j,i] -= b.conjugate()
    inv = (eye-h).inv()
    y = (2*inv-eye)*seed
    deriv = []
    if derivative:
        right = inv*seed
        for i,j,part in directions:
            a = acb(0,1) if part=="i" else acb(1)
            e = acb_mat(19,19)
            e[i,j],e[j,i] = a,-a.conjugate()
            deriv.append(2*inv*e*right)
    return y,deriv


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    original = lib("one-generator-gravity-full-ckm")
    xc,xradius = original["verify"](root,refine=True)
    rank = lib("rank-two-generator-up-quark")
    common = lib("common-embedding-central-right-ckm")
    independent = lib("independent-unitary-quark-right-holdout")
    read = lib("next-coherent-full-ckm")
    n = rank["N"]
    _,trained,_ = common["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    goal = original["targets"](root,independent,common)
    ints = read["QUATERNION"]
    norm = arb(sum(j*j for j in ints)).sqrt()
    x,y,z,w = [arb(j)/norm for j in ints]
    rotation = acb_mat([[acb(w,-z),acb(-y,-x),0,0],[acb(y,-x),acb(w,z),0,0],
                        [0,0,1,0],[0,0,0,1]])
    coords = independent["coordinates"]
    rotated = [rotation*read["block_from_coordinates"](g.entries())*adj(rotation) for g in goal]
    desired = [arb_mat([[b] for b in coords(block)]) for block in rotated]
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    proj = lib("free-four-site-native-positive-spectrum")["native"](root,n+3)
    ctx.prec = 384
    spectra = []
    for q,rows,goal in zip((-4,2),proj,desired,strict=True):
        src,dest = [spec(shell,q,k)[1] for k in (n+1,n+3)]
        cols = arb_mat([coords(a) for a in rows]).transpose()
        spectra.append((src,adj(dest),trained[n+1,q],cols,goal))

    eye = rank["eye"](19)
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    x0 = seed(root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    y0 = seed(root,rank["eye"],"ch39-next-common-ckm-cayley.tsv")
    def system(px,py,axis):
        x,dx = local(x0,px,original["DIRECTIONS"],eye,axis=="x")
        y,dy = local(y0,py,DIRECTIONS,eye,axis=="y")
        xx,yy = adj(x*x),adj(y*y)
        vals,rows = [],[]
        for iq,(src,dst,weights,cols,goal) in enumerate(spectra):
            op = yy*xx if iq==0 else y*x
            if axis=="y":
                dop = [adj(e*y+y*e)*xx for e in dy] if iq==0 else [e*x for e in dy]
            else:
                dop = [yy*adj(e*x+x*e) for e in dx] if iq==0 else [y*e for e in dx]
            ov = dst*op*src
            wr = arb_mat([[sum(v*abs(ov[i,j])**2 for j,v in enumerate(weights))] for i in range(19)])
            vals.extend((cols*wr-goal).entries())
            derivatives = []
            for e in dop:
                do = dst*e*src
                dw = arb_mat([[2*sum(v*(ov[i,j].conjugate()*do[i,j]).real for j,v in enumerate(weights))]
                              for i in range(19)])
                derivatives.append((cols*dw).entries())
            rows.extend(arb_mat(derivatives).transpose().tolist())
        return arb_mat([[v] for v in vals]),arb_mat(rows)

    yc = [arb(0)]*12
    for _ in range(2):
        f,j = system(xc,yc,"y")
        corr = j.mid().solve(f.mid())
        yc = [(v-e).mid() for v,e in zip(yc,corr.entries(),strict=True)]
    f,j = system(xc,yc,"y")
    pre = j.mid().inv()
    assert not pre.det().contains(0)
    xwide = [a+arb(0,xradius.upper()) for a in xc]
    input_radius = max(abs(b-a).upper() for a,b in zip(xc,xwide,strict=True))
    _,dx = system(xwide,yc,"x")
    seed_error = max(abs(z).upper() for z in (pre*f).entries())
    parameter_error = outer_sum(pre*dx)*input_radius
    print("参数中心、数值残差、上游参数误差",max(abs(v).upper() for v in yc).str(13),
          seed_error.str(13),parameter_error.str(13),flush=True)
    radius = 2*(seed_error+parameter_error).upper()
    wide = [a+arb(0,radius.upper()) for a in yc]
    _,dy = system(xwide,wide,"y")
    defect = arb_mat([[int(i==j) for i in range(12)] for j in range(12)])-pre*dy
    offset = (seed_error+parameter_error)/radius
    lip = outer_sum(defect)
    print("下游半宽、中心位移占比、导数行和、总和",*[a.str(13) for a in (radius,offset,lip,lip+offset)],flush=True)
    assert lip+offset<1
    # 另行读出同一对生成元的荷六三次幂，使用 N+3 实际交叉权和破缺背景。
    ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root)[2][1]
    moments,heat,fields = lib("three-index-free-atoms")["broken_fields"](root,[n+3])
    native = lib("native-rh-ray-pvm")
    _,_,cross,_ = next(row for row in native["shell_weights"](root,n+3) if row[0]==6)
    ctx.prec = 384
    x,_ = local(x0,xwide,original["DIRECTIONS"],eye,False)
    y,_ = local(y0,wide,DIRECTIONS,eye,False)
    gs,gt = [spec(shell,6,k)[1] for k in (n+1,n+3)]
    grav = rank["read"](gs,gt,[z.sqrt() for z in ordered],y*y*y*x*x*x)
    assert min(grav)>arb("1e-9")
    p0,p1 = [(arb(z)/arb("1.22089e19"))**2 for z in ("91.1876","1.176397e14")]
    weak = lib("full-neutral-free-cutoff")["weak_column"]
    thermal = [moments(fields[0],heat,f**2) for f in grav]
    columns = (*zip(*thermal),[weak(f**2,p0,p1) for f in grav])
    k,w,d = [sum(a.real*b for a,b in zip(cross[0],col,strict=True)) for col in columns]
    assert d<arb("-.07")
    print("在同一味解盒内的荷六第16层和 K、W、D",*[a.str(12) for a in (grav[15],k,w,d)],flush=True)
    return xc,xradius,yc,radius


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
