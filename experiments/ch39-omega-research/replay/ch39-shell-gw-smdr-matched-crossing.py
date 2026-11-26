"""Match an SMDR low-energy flow to the shell's nonlinear two-loop gauge crossing."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb, ctx

folder = Path(__file__).resolve().parent
model = run_path(str(folder/"ch39-shell-gw-nonlinear-crossing.py"))
smdr = run_path(str(folder/"ch39-shell-gw-smdr-quartic.py"))
density = run_path(str(folder/"ch39-shell-gw-positive-density-cp.py"))
flow,fit,q,loop = [model[key] for key in ("flow","fit","q","loop")]
ctx.cap = 5
mz = q("91.1876")


def get_boundary(first_low,first_high):
    # Once the common scale is above the lower end, the small tail is bounded
    # by a verified tube, retaining correlations in the gauge correction.
    beta = smdr["evaluator"](5)
    y = smdr["init"](".967",".55853")
    h = q(1,1536)
    x = (first_low/mz).log()/loop
    N = int(float((x/h).lower()))
    assert N*h<x<(N+1)*h
    for _ in range(N):
        y = smdr["advance"](beta,y,h)
    y = smdr["advance"](beta,y,x-N*h)
    def check_hierarchy(ys):
        top,bottom = ys[3],ys[4]
        assert max((ys[9]/ys[7])**2,(ys[7]/top)**2,(ys[10]/ys[8])**2,(ys[8]/bottom)**2)<q(".004")
    end = (first_high/first_low).log()/loop
    tube = [v+arb(0,".02") for v in y]
    speed = beta(tube)
    assert all(abs(end*v)<q(".019") for v in speed)
    tube = [v+q(0).union(end)*v1 for v,v1 in zip(y,speed)]
    check_hierarchy(tube)
    # Directly integrate the correction to the one-loop reference from first_low.
    gauges = [q(5,3)*y[2]**2, y[1]**2, y[0]**2]
    one = flow["gauge"]((first_low/1000).log()/loop,(flow["start"],flow["beta"]))
    z0 = [4*arb.pi()*(1/g-1/v) for g,v in zip(gauges,one)]
    slope = beta(tube)
    for i,(j,normal) in enumerate(((2,q(5,3)),(1,q(1)),(0,q(1)))):
        dz = 8*arb.pi()*(flow["beta"][i]-slope[j]/(normal*tube[j]**3))
        z0[i] += q(0).union(end)*dz
    return tube,z0


def matrices(y):
    def diag(values):
        M = acb_mat(3,3)
        for i,t in enumerate(values):
            M[i,i] = t
        return M
    V,_,_ = fit["orbit"]["reference"]()
    return [diag((y[9],y[7],y[3])),V*diag((y[10],y[8],y[4])),diag((y[12],y[11],y[5]))]


def check_angle_bounds(ceiling):
    # Use the tabulated CKM angles at 1 TeV and the one-loop mixing equation
    # with the SMDR eigenvalues. Its norm estimates only require broad bounds.
    fine = smdr["evaluator"](5)
    y = smdr["init"](".967",".55853")
    h = q(1,1536)
    interval = (ceiling/mz).log()/loop
    N = int(float((interval/h).upper()))+1
    V = fit["orbit"]["reference"]()[0]
    radius = q(".023")
    for i in range(N):
        dy = fine([t+arb(0,".022") for t in y])
        assert all(abs(h*b)<q(".0219") for b in dy)
        a = [t+q(0).union(h)*v for t,v in zip(y,dy)]
        if (i+1)*h > (q(1000)/mz).log()/loop:
            assert a[3]<q(".87") and a[4]<q(".015") and a[7]<q(".0033")
            assert max((a[9]/a[7])**2,(a[7]/a[3])**2,(a[10]/a[8])**2,(a[8]/a[4])**2)<q(".004")
        y = smdr["advance"](fine,y,h)
    # Each off-diagonal is controlled by the products of the third-row
    # moduli; the opposite direction uses the conservative bottom² bound.
    third = [abs(V[2,j]).abs_upper()+radius for j in (0,1)]+[q(1)]
    coupling = q("1.5")*(1+q(".004"))/(1-q(".004"))
    off = sum((q(".87")**2*third[i]*third[j]+q(".0033")**2)**2
              for i,j in ((0,1),(0,2),(1,2)))
    speed = coupling*((2*off).sqrt()+q(6).sqrt()*q(".015")**2)
    assert speed*(ceiling/1000).log()/loop < radius
    return radius


def check(n):
    Q,m = flow["masses"](n)
    width = q("3e-8" if n==43 else "1e-8")
    e = density["window"](n,width)
    gap = (m[0]**2-125**2)/q("1.2e12")**2
    eps = 4*width/gap
    stretch = ((1+eps)/(1-eps)).sqrt()
    mass = [v*(1/stretch).union(stretch) for v in m]
    s0,s1,t0,t1 = map(q,(".36",".391","-.139","-.112"))
    s,t = s0.union(s1),t0.union(t1)
    first = model["shifted"](mass[0],s)
    angle = check_angle_bounds(first.upper())
    y,z = get_boundary(first.lower(),first.upper())
    data = matrices(y)
    # For each allowed s use its true CKM. Its matrix differs from the
    # reference by at most radius*y_b in Hilbert-Schmidt norm.
    mid = matrices([v.mid() for v in y])
    allowance = angle*y[4]
    trace = model["continuous_light"](n,data,mid,e,12*width,allowance)
    point = (data,z)  # Unitary CKM motion preserves the total Frobenius norm.
    def face(ss,tt):
        return model["run_face"](Q,mass,ss,tt,(s0,s1),point,trace)
    faces = [face(s,t0)[0],face(s,t1)[0],face(s0,t)[1],face(s1,t)[1]]
    print("zero",n,"faces",[(float(v.lower()),float(v.upper())) for v in faces])
    assert faces[0]>0 and faces[1]<0 and faces[2]<0 and faces[3]>0
    _,_,inv,cap = face(s,t)
    assert q("38.82")<inv<q("38.98") and cap<q("2.32")
    if n==43:
        assert q("7.41e10")<first<q("7.66e10") and y[6]>q(".00091")
    else:
        assert q("1.80e11")<first<q("1.861e11") and y[6]<q("-.00025")
    print("first heavy GeV",first.str(9),"light quartic",y[6].str(9),
          "unified inverse",inv.str(8),"Yukawa sum upper",cap.upper().str(8))


if __name__=="__main__":
    for n in (43,195):
        check(n)
