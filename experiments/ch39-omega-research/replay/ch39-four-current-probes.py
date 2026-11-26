"""在连续耦合区间覆盖全部较小零点编号，再以正性排除较大编号。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

three = run_path(str(Path(__file__).with_name("ch39-three-probe-spectral-test.py")))
mat,poles,q,shape,form = (three[key] for key in ("mat","poles","q","shape","form"))
basis = arb_mat([[3,0,q("1/6")],[-1,0,q("1/2")],[0,1,0]])
eye = arb_mat([[1,0,0],[0,1,0],[0,0,1]])
family = arb_mat([[0],[1],[0]])
sources = [arb_mat([[3],[-1],[2]]), three["probe"]]


def exchanges(kbox, geometry, h, k0, tolerance):
    rt = kbox.sqrt()
    cols = arb_mat([[1,0,0],[0,1,0],[0,0,rt]])
    vac = arb_mat([[rt,0,0],[0,rt,0],[0,0,1]])
    matrices = []
    limits = []
    hyp = mat["hypercharge"]*h
    for i,(p,a,pb,ab,refx,refj,coef) in enumerate(geometry):
        transformed = cols*pb*cols
        transformed[2,2] += h*q("25/324")
        if i==0:
            matrices.append(transformed)
            continue
        upd = eye-transformed*vac*ab*vac
        if upd.det().contains(0):
            return None
        g = upd.inv()*transformed
        matrices.append(g)
        w = (eye-a*(p*kbox+hyp)).inv()*a
        col = arb_mat([[p[j,2]] for j in range(3)])
        difference = (k0-kbox)*(coef*(col+p*w*col*kbox))[0,0]
        limits.append(difference+(-1 if i==1 else 1)*tolerance*refj)
    return matrices, limits


def test_box(kbox, geometry, h, k0, bound, cap, source_index):
    # Scale the two hypercharge-orthogonal columns by 1/sqrt(k).
    rt = kbox.sqrt()
    general = source_index==1
    coords = basis.inv()*sources[source_index]
    if general:
        coords[2,0] /= rt
    else:
        assert abs(coords[2,0]) < q("1e-35")
        coords[2,0] = 0
    bundle = exchanges(kbox, geometry, h, k0, q("2e-8"))
    if bundle is None:
        return False
    matrices, limits = bundle
    if not all(x>0 for x in limits):
        return False
    zero,first,second = matrices
    metric = second-first*cap
    if not mat["positive"](metric):
        return False
    inv = metric.inv()
    u,v,w = first*family,second*family,second*coords
    axes = arb_mat([[v[i,0],w[i,0]] for i in range(3)])
    # Fix scalar coefficients at a midpoint. All subsequent tests use the
    # interval matrices, including the error of the projected direction.
    mid = lambda m: arb_mat([[m[i,j].mid() for j in range(m.ncols())] for i in range(m.nrows())])
    trans = mid(axes).transpose()*mid(inv)
    gram = trans*mid(axes)
    if gram.det().contains(0):
        return False
    weights = mid(gram.inv()*trans*mid(u))
    alpha,beta = weights.entries()
    remainder = u-axes*weights
    loss = max(arb(0),form(remainder,inv,remainder).upper()).sqrt()
    slack = (limits[0]/cap).sqrt()-abs(alpha)*limits[1].sqrt()-loss
    if not slack>0 or beta.contains(0):
        return False
    # Every admissible updated exchange at x2 has this much excess on the probe.
    need = (slack/abs(beta))**2
    origin = form(coords,zero,coords)
    lower = form(coords,second,coords)+need
    return origin>0 and lower>origin*geometry[2][4][source_index]*(1+bound)


def cover(lo,hi,geometry,h,k0,bound,cap,source_index):
    stack = [(lo,hi)]
    count = 0
    while stack:
        left,right = stack.pop()
        if test_box(left.union(right),geometry,h,k0,bound,cap,source_index):
            count += 1
            continue
        if right-left < q("1e-9"):
            raise ArithmeticError(f"unresolved coupling interval {left} to {right}")
        middle = ((left+right)/2).mid()
        stack.extend(((left,middle),(middle,right)))
    return count


def setup(logs, time, targets):
    (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
    neutral = poles["reference_readings"](logs)[3]
    e,b,d,k0 = neutral**2/scale,c/rho,eta/rho,1/rho
    h = mat["flow"]["running"]["hypercharge_window"](time+neutral.log()).mid()**2
    ks = [1/(poles["offset"]+goal/(15*scale)) for goal in targets]
    raw = [(mat["massive_exchange"](e,z,b,d),arb_mat(3,3),q(1))]
    raw += [shape["probe"](log_x,e,z,b,d,masses[0]) for log_x in (q(".5"),q(1))]
    geometry = []
    bas = basis.inv()
    for p,a,_ in raw:
        t = p*k0+mat["hypercharge"]*h
        g = (eye-t*a).inv()*t
        origin = raw[0][0]*k0+mat["hypercharge"]*h
        refx = [form(s,g,s)/form(s,origin,s) for s in sources]
        col = arb_mat([[p[j,2]] for j in range(3)])
        coef = col.transpose()*(eye-a*t).inv()*a
        geometry.append((p,a,basis.transpose()*p*basis,bas*a*bas.transpose(),
                         refx,g[2,2]/k0,coef))
    x1,x2 = raw[1][2],raw[2][2]
    threshold = 4*masses[0]**2/e
    cap = (x1-1)/(x2-1)*(threshold+x2)/(threshold+x1)
    return geometry,h,k0,ks,cap


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (1,11430,11432)]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        geometry,h,k0,ks,cap = setup(logs,time,targets)
        bound = q(".00005")
        counts = [cover(ks[0].lower(),q(2),geometry,h,k0,bound,cap,0),
                  cover(q(2),ks[1].upper(),geometry,h,k0,bound,cap,1)]
        p,a,_,_,_,refj,_ = geometry[2]
        t32 = p*ks[2]+mat["hypercharge"]*h
        family32 = ((eye-t32*a).inv()*t32)[2,2]/ks[2]
        assert family32 > refj*(1+q("2e-8"))
        print("case",case,"certified intervals",counts,"next higher label excess (ppm)",
              ((family32/refj-1)*10**6).lower().str(10),flush=True)


if __name__=="__main__":
    verify()
