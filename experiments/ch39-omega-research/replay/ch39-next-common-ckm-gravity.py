"""从 N+1 经上游完整见证 X，以同一新 Y 同时实现 N+3 的完整味与荷六三读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


DIRECTIONS = ((3,0,"i"),(14,13,"i"),(14,14,"i"),(14,6,"r"),(14,8,"r"),(11,3,"i"),
              (14,7,"r"),(18,5,"r"),(17,6,"r"),(17,14,"i"),(9,6,"i"),(17,7,"r"),
              (16,7,"r"),(6,3,"r"),(8,8,"i"),(10,8,"r"),(14,2,"r"),(16,8,"r"))


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    upstream = lib("one-generator-gravity-full-ckm")
    cx,rx = upstream["verify"](root,refine=True)
    common = lib("next-common-ckm")
    local,norm = common["local"],common["outer_sum"]
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    embed = lib("common-embedding-central-right-ckm")
    independent = lib("independent-unitary-quark-right-holdout")
    reference = lib("next-coherent-full-ckm")
    _,trained,_ = embed["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    goals = upstream["targets"](root,independent,embed)
    ints = reference["QUATERNION"]
    length = arb(sum(i*i for i in ints)).sqrt()
    a,b,c,d = [arb(i)/length for i in ints]
    rot = acb_mat([[acb(d,-c),acb(-b,-a),0,0],[acb(b,-a),acb(d,c),0,0],
                   [0,0,1,0],[0,0,0,1]])
    coords = independent["coordinates"]
    blocks = [rot*reference["block_from_coordinates"](g.entries())*adj(rot) for g in goals]
    rows = lib("free-four-site-native-positive-spectrum")["native"](root,n+3)
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2][1]
    assert all(v is ordered[0] for v in ordered[:17])
    low,h,k = [ordered[j].sqrt() for j in (0,17,18)]
    moments,heat,fields = lib("three-index-free-atoms")["broken_fields"](root,[n+3])
    field = fields[0]
    native = lib("native-rh-ray-pvm")
    _,_,cross,_ = next(row for row in native["shell_weights"](root,n+3) if row[0]==6)
    c = [v.real for v in cross[0]]
    ctx.prec = 384
    weak = native["weak_series"](56)
    weak_high = lib("full-neutral-free-cutoff")["weak_column"]
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    constant = arb("-.074256577194663365545")
    left,right = arb("1.28"),arb(2)
    def approximate_weak(t):
        return c[15]*weak(t)+(1-c[15])*constant
    assert approximate_weak(left)>0>approximate_weak(right)
    for _ in range(42):
        mid = (left+right)/2
        if approximate_weak(mid)>0:
            left = mid
        else:
            assert approximate_weak(mid)<0
            right = mid
    left -= arb("1e-11")
    right += arb("1e-11")
    time = left.union(right)
    spectra = []
    for q,group,block in zip((-4,2),rows,blocks,strict=True):
        src,dst = [spec(shell,q,j)[1] for j in (n+1,n+3)]
        cols = arb_mat([coords(row) for row in group]).transpose()
        target = arb_mat([[v] for v in coords(block)])
        spectra.append((src,adj(dst),cols,trained[n+1,q],target))
    gs,gt = [spec(shell,6,j)[1] for j in (n+1,n+3)]
    gs = acb_mat([[gs[i,17],gs[i,18]] for i in range(19)])
    gt = adj(gt)
    eye = rank["eye"](19)
    get_seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    x0 = get_seed(root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    y0 = get_seed(root,rank["eye"],"ch39-next-common-ckm-gravity-cayley.tsv")
    xm,_ = local(x0,cx,upstream["DIRECTIONS"],eye,False)
    link = (gt*y0*y0*y0*xm*xm*xm*gs)[15,0]
    phase = link/abs(link)

    def system(sx,sy,t,axis):
        x,dx = local(x0,sx,upstream["DIRECTIONS"],eye,axis=="x")
        y,dy = local(y0,sy,DIRECTIONS,eye,axis=="y")
        x2,y2 = x*x,y*y
        val,jac = [],[]
        for iq,(src,dst,cols,weights,goal) in enumerate(spectra):
            op = (adj(y2)*adj(x2) if iq==0 else y*x)
            if axis=="x":
                deriv = ([adj(y2)*adj(e*x+x*e) for e in dx] if iq==0 else [y*e for e in dx])
            else:
                deriv = ([adj(e*y+y*e)*adj(x2) for e in dy] if iq==0 else [e*x for e in dy])
            ov = dst*op*src
            mass = arb_mat([[sum(w*abs(ov[i,j])**2 for j,w in enumerate(weights))] for i in range(19)])
            val.extend((cols*mass-goal).entries())
            derivatives = []
            for e in deriv:
                dov = dst*e*src
                dw = arb_mat([[2*sum(w*(ov[i,j].conjugate()*dov[i,j]).real for j,w in enumerate(weights))]
                              for i in range(19)])
                derivatives.append((cols*dw).entries())
            jac.extend(arb_mat(derivatives).transpose().tolist())
        x3,y3 = x2*x,y2*y
        ov = gt*y3*x3*gs
        if axis=="x":
            dov = [gt*y3*(e*x2+x*e*x+x2*e)*gs for e in dx]
        else:
            dov = [gt*(e*y2+y*e*y+y2*e)*x3*gs for e in dy]
        cutoff = mz*t.sqrt()/mp
        eta = ((cutoff-low)/(h-low)).sqrt()
        res = ov[15,0]-eta*phase
        val.extend((res.real,res.imag,ov[15,1].real,ov[15,1].imag))
        jac.extend([[getattr(d[15,j],part) for d in dov] for j in (0,1) for part in ("real","imag")])
        f = [low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2 for i in range(19)]
        f[15] = cutoff
        err = (abs(c[15])*cutoff**2*sum(abs(a) for _,a,_ in field)/(2*arb.pi()),
               abs(c[15])*cutoff**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi()))
        thermal = [arb(0,e.upper())-goal for e,goal in zip(err,(arb(1),arb("1e-122")),strict=True)]
        slopes = [[],[]]
        for i,v in enumerate(f):
            if i==15:
                for col in slopes:
                    col.append(arb(0))
                continue
            assert v>arb("1e-6")
            for j,z in enumerate(moments(field,heat,v*v)):
                thermal[j] += c[i]*z
            decay = [(-z/v**2).exp() for _,_,z in field]
            slopes[0].append(c[i]*v*sum(a*b for (_,a,_),b in zip(field,decay,strict=True))/arb.pi())
            slopes[1].append(-c[i]*v**3*sum(g*b for (g,_,_),b in zip(field,decay,strict=True))/(2*arb.pi()))
        df = arb_mat([[2*((h-low)*(ov[i,0].conjugate()*d[i,0]).real+
                          (k-low)*(ov[i,1].conjugate()*d[i,1]).real) for d in dov] for i in range(19)])
        val.extend(thermal)
        jac.extend((arb_mat(slopes)*df).tolist())
        return arb_mat([[v] for v in val]),arb_mat(jac),f

    cy = [arb(0)]*18
    for _ in range(3):
        val,jac,_ = system(cx,cy,time.mid(),"y")
        correction = jac.mid().solve(val.mid())
        cy = [(v-c).mid() for v,c in zip(cy,correction.entries(),strict=True)]
    val,jac,_ = system(cx,cy,time,"y")
    pre = jac.mid().inv()
    assert not pre.det().contains(0)
    xwide = [v+arb(0,rx.upper()) for v in cx]
    xwidth = max(abs(v-z).upper() for v,z in zip(cx,xwide,strict=True))
    _,dx,_ = system(xwide,cy,time,"x")
    err = max(abs(v).upper() for v in (pre*val).entries())+norm(pre*dx)*xwidth
    radius = 2*err.upper()
    ywide = [v+arb(0,radius.upper()) for v in cy]
    _,dy,amps = system(xwide,ywide,time,"y")
    defect = arb_mat([[int(i==j) for j in range(18)] for i in range(18)])-pre*dy
    lip = norm(defect)
    offset = err/radius
    print("下游中心最大绝对坐标、盒半宽、中心位移占比与导数界",
          *[v.str(13) for v in (max(abs(z).upper() for z in cy),radius,offset,lip)],flush=True)
    assert lip+offset<1
    other = sum(c[i]*weak_high(v*v,p0,p1) for i,v in enumerate(amps) if i!=15)
    ends = [c[15]*weak(t)+other for t in (left,right)]
    assert ends[0]>0>ends[1]
    print("末壳第16层交叉权、截止 GeV、弱差两端",
          c[15].str(12),(mz*time.sqrt()).str(15),*[v.str(12) for v in ends],flush=True)
    print("第18、19层振幅及其余中高层最小振幅",amps[17].str(12),amps[18].str(12),
          min(v.lower() for i,v in enumerate(amps) if i!=15).str(12),flush=True)
    return cx,rx,cy,radius,time


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
