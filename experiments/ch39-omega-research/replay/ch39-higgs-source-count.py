"""分类一个颜色源加 n 份希格斯方向的两动量标量识别协议。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

scheme = run_path(str(Path(__file__).with_name("ch39-rational-momentum-calibration.py")))
cal, poles, q = (scheme[key] for key in ("cal", "poles", "q"))
form = cal["form"]
j = arb_mat([[0],[0],[1]])


def exclude(interval, rows, k0, cap, required, source):
    bounds = cal["exchange_bounds"](interval,rows,k0,cap)
    if bounds is None:
        return False
    g1, g2, (low, high), metric = bounds
    mixed_room = required-form(source,g2,source)
    if mixed_room < 0:
        return True
    if not mixed_room > 0:
        return False
    inv = metric.inv()
    u,v,z = g1*cal["family"],g2*cal["family"],g2*source
    axes = arb_mat([[v[i,0],z[i,0]] for i in range(3)])
    def mid(a):
        return arb_mat([[a[i,j].mid() for j in range(a.ncols())] for i in range(a.nrows())])
    basis = mid(axes)
    gram = basis.transpose()*mid(inv)*basis
    rhs = basis.transpose()*mid(inv)*mid(u)
    for alpha in (q("39/40"),q("9/10"),q("3/4"),q("1/2")):
        diagonal = arb_mat([[high.upper()/alpha,0],[0,mixed_room.upper()/(1-alpha)]])
        projector = gram+diagonal
        if projector.det().contains(0):
            continue
        coeff = mid(projector.inv()*rhs)
        rem = u-axes*coeff
        if cap*(form(rem,inv,rem)+form(coeff,diagonal,coeff)) < low:
            return True
    return False


def nominal_graphs(logs,time,target):
    (length,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
    e = scheme["anchor"]/length
    k0 = 1/rho
    k = 1/(poles["offset"]+target/(15*length))
    h = scheme["mat"]["flow"]["running"]["hypercharge_window"](time+scheme["anchor"].log()/2).mid()**2
    old, ref = [], []
    for log_x in scheme["logs_x"]:
        p,a,_ = cal["shape"]["probe"](log_x,e,z,c/rho,eta/rho,masses[0])
        bare = p+scheme["mat"]["hypercharge"]*(h/k0)
        assert all(cal["mat"]["positive"](a.inv()-bare*coupling) for coupling in (k,k0))
        old.append(cal["dressed"](bare,a,k))
        ref.append(cal["dressed"](bare,a,k0))
    x1,x2 = (v.exp() for v in scheme["logs_x"])
    gap = 4*masses[0]**2/e
    cap = (x1-1)/(x2-1)*(gap+x2)/(gap+x1)
    return old, ref, cap


def cross(a,b):
    return arb_mat([[a[(i+1)%3,0]*b[(i+2)%3,0]-a[(i+2)%3,0]*b[(i+1)%3,0]]
                    for i in range(3)])


def witness(graphs,source,sign=1,errors=False):
    (g1,g2),(r1,r2),cap = graphs
    ef, em = (q("1.9e-8"),q("1.9e-5")) if errors else (q(0),q(0))
    first = (r1[2,2]-g1[2,2])-ef*r1[2,2]
    second = (r2[2,2]-g2[2,2])+ef*r2[2,2]
    third = form(source,(r2-g2)+em*r2,source)
    if not (first > 0 and second > 0 and third > 0):
        return False
    u,v,z = g1*j,g2*j,g2*source
    # Force the two x2 projections exactly. A free line parameter remains.
    axis = cross(v,z)
    h = max(range(3), key=lambda i: float(abs(axis[i,0].mid())))
    pivot = axis[h,0]
    if pivot.contains(0):
        return False
    i,jj = (h+1)%3,(h+2)%3
    w0 = arb_mat(3,1)
    w0[i,0] = (second.sqrt()*z[jj,0]-sign*third.sqrt()*v[jj,0])/pivot
    w0[jj,0] = (sign*third.sqrt()*v[i,0]-second.sqrt()*z[i,0])/pivot
    axis = axis*(1/pivot)
    metric = g2-cap*g1
    au = (u.transpose()*axis)[0,0]
    bu = (u.transpose()*w0)[0,0]
    a = cap*au**2-first*form(axis,metric,axis)
    b = cap*au*bu-first*form(axis,metric,w0)
    c = cap*bu**2-first*(1+form(w0,metric,w0))
    return b*b-a*c > 0


def cover_tail(graphs):
    # Dividing the source by n does not change a normalized reading.
    stack = [(q(0),q("1/11").upper())]
    count = 0
    while stack:
        left,right = stack.pop()
        source = arb_mat([[3+left.union(right)],[-1],[2]])
        if witness(graphs,source):
            count += 1
        else:
            if right-left < q("1e-10"):
                raise ArithmeticError(f"unresolved inverse source count {left}, {right}")
            middle = ((left+right)/2).mid()
            stack.extend(((left,middle),(middle,right)))
    return count


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11432)]
    source10 = cal["inv_basis"]*arb_mat([[31],[-10],[20]])
    for i, logs in enumerate(poles["retune"]["log_sets"]):
        (records,k0,ks,cap,_),_ = scheme["collect"](logs,time,targets,q("1e-11"),q("1e-10"),source10)
        count,_ = cal["identify"](records,k0,ks,cap,source10,exclude)
        graphs = nominal_graphs(logs,time,targets[0])
        assert cal["mat"]["positive"](graphs[0][1]-graphs[0][0]*graphs[2])
        for n in range(9):
            source = arb_mat([[3*n+1],[-n],[2*n]])
            assert witness(graphs,source,-1 if n==8 else 1,n==8)
        tail = cover_tail(graphs)
        print("case",i,"ten Higgs charges: lower coupling boxes",count,
              "counterexamples at 0..8 and all n>=11; continuous tail boxes",tail,flush=True)


if __name__ == "__main__":
    verify()
