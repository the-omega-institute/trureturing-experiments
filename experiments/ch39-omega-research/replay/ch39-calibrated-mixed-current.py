"""用减法点电流比和四个有限精度质量校准混合谱识别。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

s = run_path(str(Path(__file__).with_name("ch39-three-probe-spectral-test.py")))
mat, poles, shape, q, form = (s[key] for key in ("mat","poles","shape","q","form"))
# Columns are the family source, the mixed source, and the calibration source.
basis = arb_mat([[0,28,1],[0,-9,3],[1,18,0]])
mixed, family = arb_mat([[0],[1],[0]]), arb_mat([[1],[0],[0]])
calibration = arb_mat([[0],[0],[1]])
mapped_y = basis.transpose()*mat["y"]
hypercharge = mapped_y*mapped_y.transpose()
inv_basis = basis.inv()
eye = arb_mat([[1,0,0],[0,1,0],[0,0,1]])


def uncertain(error):
    return 1+arb(0,error.upper())


def dressed(tree, a, k):
    return (eye-tree*a*k).inv()*tree


def projected_tree(energy,z,b,d):
    # Polarize the sum-of-squares Schur complement before interval evaluation.
    weights = (poles["colour"](-energy,mat["unit"])/24,
               poles["right"](-energy,mat["unit"])/4, 1/(b+d*energy))
    vev = (3,-1,2)
    den = z+sum(w*n*n for w,n in zip(weights,vev,strict=True))
    total = arb_mat(3,3)
    for i in range(3):
        row = arb_mat([[basis[i,j]] for j in range(3)])
        total += row*row.transpose()*(z*weights[i])
        for m in range(i):
            wedge = arb_mat([[basis[i,j]*vev[m]-basis[m,j]*vev[i]] for j in range(3)])
            total += wedge*wedge.transpose()*(weights[i]*weights[m])
    return total*(energy/den)-hypercharge*q("1/24")


def collect(logs, time, targets, mass_error):
    (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
    c1,c2,_,n1,_ = poles["reference_readings"](logs)
    b,d,k0 = c/rho,eta/rho,1/rho
    n2 = mat["geometry"]["second_neutral_mass"](scale,z,b,d)
    e = n1**2/scale
    h0 = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
    ref0 = projected_tree(e,z,b,d)+hypercharge*(h0/k0)
    ratio = form(calibration,ref0,calibration)/ref0[0,0]
    observed = [v*uncertain(mass_error) for v in (c1,c2,n1,n2)]
    L,zz,_ = poles["charged_geometry"](*observed[:2])
    bb,dd = mat["geometry"]["neutral_shape"](L,zz,observed[2:])
    ee = observed[2]**2/L
    triplet,_ = mat["geometry"]["physical_masses"](bb,q(1),zz)
    origin = projected_tree(ee,zz,bb,dd)
    t = (ratio*uncertain(q("1e-6"))*origin[0,0]-form(calibration,origin,calibration))
    t /= form(calibration,hypercharge,calibration)
    assert t>0
    origin += hypercharge*t
    return probe_records((e,z,b,d,masses[0],h0,k0,ref0),
                         (ee,zz,bb,dd,triplet,t,origin,L), targets, (q(".5"),q(1)))


def probe_records(reference, fitted, targets, log_ratios, source=mixed):
    e,z,b,d,nominal_triplet,h0,k0,ref0 = reference
    ee,zz,bb,dd,triplet,t,origin,L = fitted
    records = []
    for log_x in log_ratios:
        _,a0,x = shape["probe"](log_x,e,z,b,d,nominal_triplet)
        direct = dressed(projected_tree(e*x,z,b,d)+hypercharge*(h0/k0),inv_basis*a0*inv_basis.transpose(),k0)
        dest = direct[0,0]*origin[0,0]/ref0[0,0]
        _,a,_ = shape["probe"](log_x,ee,zz,bb,dd,triplet)
        tree = projected_tree(ee*x,zz,bb,dd)+hypercharge*t
        a = inv_basis*a*inv_basis.transpose()
        base = dressed(tree,a,k0)
        col = tree*family
        coef = col.transpose()*(eye-a*tree*k0).inv()*a
        cross_target = form(source,direct,source)*form(source,origin,source)/form(source,ref0,source)
        records.append((tree,a,col,coef,dest-base[0,0],dest,cross_target))
    x1,x2 = (log_x.exp() for log_x in log_ratios)
    threshold = 4*triplet**2/ee
    cap = ((x1-1)/(x2-1)*(threshold+x2)/(threshold+x1)).upper()
    ks = [1/(poles["offset"]+g/(15*L)) for g in targets]
    return records, k0, ks, cap, t


def exchange_bounds(interval, records, k0, cap):
    graphs, limits = [], []
    for i,(tree,a,col,coef,offset,target,_) in enumerate(records):
        inv = eye-tree*a*interval
        if inv.det().contains(0):
            return None
        g = inv.inv()*tree
        w = (eye-a*tree*interval).inv()*a
        gap = offset+(k0-interval)*(coef*(col+tree*w*col*interval))[0,0]
        graphs.append(g)
        limits.append(gap+(-1 if i==0 else 1)*q("2e-8")*target)
    if not all(v>0 for v in limits):
        return None
    first,second = graphs
    metric = second-first*cap
    if not mat["positive"](metric):
        return None
    return first, second, limits, metric


def separate(interval, records, k0, cap, required, source=mixed):
    bounds = exchange_bounds(interval,records,k0,cap)
    if bounds is None:
        return False
    first,second,limits,metric = bounds
    inv = metric.inv()
    u,v,w = first*family,second*family,second*source
    axes = arb_mat([[v[i,0],w[i,0]] for i in range(3)])
    mid = lambda m: arb_mat([[m[i,j].mid() for j in range(m.ncols())] for i in range(m.nrows())])
    trans = mid(axes).transpose()*mid(inv)
    gram = trans*mid(axes)
    if gram.det().contains(0):
        return False
    weights = mid(gram.inv()*trans*mid(u))
    alpha,beta = weights.entries()
    r = u-axes*weights
    noise = max(arb(0),form(r,inv,r).upper()).sqrt()
    slack = (limits[0]/cap).sqrt()-abs(alpha)*limits[1].sqrt()-noise
    if not slack>0 or beta.contains(0):
        return False
    return form(source,second,source)+(slack/abs(beta))**2 > required


def identify(records,k0,ks,cap,source=mixed,separator=separate):
    required = (records[1][-1]*(1+q(".00002"))).upper()
    stack = [(q(0),ks[0].upper())]
    count = 0
    while stack:
        left,right = stack.pop()
        try:
            certified = separator(left.union(right),records,k0,cap,required,source)
        except ZeroDivisionError:
            certified = False
        if certified:
            count += 1
        else:
            if right-left < q("1e-9"):
                raise ArithmeticError(f"unresolved coupling {left} to {right}")
            middle = ((left+right)/2).mid()
            stack.extend(((left,middle),(middle,right)))
    tree,a,_,_,_,dest,_ = records[1]
    high = dressed(tree,a,ks[1])[0,0]/dest-1
    assert high > q("2e-8")
    return count, high


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11432)]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        records,k0,ks,cap,t = collect(logs,time,targets,q("1e-10"))
        count, high = identify(records,k0,ks,cap)
        print("case",case,"calibrated gY^2/k",t.str(9),"covered boxes",count,
              "higher-label excess (ppm)", (high*10**6).lower().str(9),flush=True)


if __name__=="__main__":
    verify()
