"""在重谱 2345712 的三个原有支点份额上认证连续交汇能区。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
TARGET = 2345712
LINES = (16,285,301)


def setup():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,choices = repair["source_data"]()
    _,lib,(rows,) = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,(TARGET,),tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    beta = [v/(2*arb.pi()) for v in model["beta"]]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    q = [1-t,t,-1]
    projection = arb_mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,*q],
                          [0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    coeff = projection*kernel["series"](rows,fields)
    global_source = load("two-share-global-forward")
    scalar,evaluate = global_source["scalar"],global_source["evaluate"]
    # 下面只在 s>=1 上调用，故 z/(s*x)<1。
    assert max(z.upper() for _,_,z in rows)<min(arb(1),arb(".007")*min(a.lower() for a,_,_ in edges))
    n = kernel["ORDER"]+1
    factor = 2/(4*arb.pi()*factorial(n))
    sizes = [sum(abs(2*a) for _,a,_ in rows),sum(abs(g) for g,_,_ in rows),
             sum(abs(f*b) for bs in fields for f,b in zip(q,bs,strict=True)),
             *[sum(abs(bs[k]) for bs in fields) for k in range(3)]]
    weights = [sum(abs(w) for _,_,w in edges),*[2*abs(edges[i][2]) for i in LINES]]
    measures = [moments,*[choices[i][0] for i in LINES]]
    data = []
    for j,size in enumerate(sizes):
        row = arb_mat([[coeff[j,k] for k in range(coeff.ncols())]])
        data.append([(scalar(row,m),arb(0,(factor*size*w).upper())) for m,w in zip(measures,weights,strict=True)])
    target = [arb(1),arb("1e-122"),-sum(f*a for f,a in zip(q,initial,strict=True))]

    def point(s):
        assert s>=1
        invs,logarithm = 1/s,s.log()
        vals = [[evaluate(co,err,s,invs,logarithm,0) for co,err in row] for row in data]
        matrix = arb_mat([a[1:] for a in vals[:3]])
        rhs = arb_mat([[b-a[0]] for b,a in zip(target,vals[:3],strict=True)])
        if matrix.det().contains(0):
            return None
        moves = matrix.solve(rhs).entries()
        shares = [choices[i][1]+v for i,v in zip(LINES,moves,strict=True)]
        couplings = [a+row[0]+sum(b*v for b,v in zip(row[1:],moves,strict=True)) for a,row in zip(initial,vals[3:],strict=True)]
        log_energy = (couplings[0]-couplings[2])/(beta[0]-beta[2])
        common = couplings[2]-beta[2]*log_energy
        floor = min(v.lower() for v in (*couplings,common))
        return model["mz"]*log_energy.exp(),shares,floor
    return point,model["mz"]


def verify():
    point,mz = setup()
    left,right = arb("1.000605"),arb("1.0042")
    high,_,_ = point(left)
    low,_,_ = point(right)
    print("左端与右端能标（GeV）",high.str(13),low.str(13),flush=True)
    assert arb("3.06225e15")<high<arb("3.06227e15") and arb("70.86765")<low<arb("70.86767")<mz
    # 所有相邻区间均使用同一精确线性系统。可逆性给连续的唯一份额函数。
    pending = [(left,right)]
    certified = []
    while pending:
        a,b = pending.pop()
        span = a.union(b)
        values = point(span)
        if values:
            energy,shares,floor = values
            if all(0<v<1 for v in shares) and floor>arb(".5") and energy<arb("1.22089e19"):
                certified.append((min(v.lower() for v in shares),max(v.upper() for v in shares),floor))
                continue
        assert b-a>arb("1e-12"),(span,values)
        mid = ((a+b)/2).mid()
        pending.extend(((a,mid),(mid,b)))
    print("认证盒数、最小份额下界、最大份额上界、沿程最小逆耦合下界",len(certified),
          min(s[0] for s in certified).str(10),max(s[1] for s in certified).str(10),
          min(s[2] for s in certified).str(10),flush=True)


if __name__=="__main__":
    verify()
