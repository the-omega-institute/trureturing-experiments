"""在重谱 2345712 的三个原有支点份额上认证连续交汇能区。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
TARGET = 2345712
LINES = (16,285,301)
EXCLUDED = object()


def setup():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    global_source = load("two-share-global-forward")
    kernel,moments,choices,edges,fields,bgs,q,v0,total,baseline,increments,loss = global_source["context"]((TARGET,),20)
    small = (baseline+sum(increments[i] for i in LINES))/20+loss
    assert small<arb(".76")
    rows = bgs[TARGET]
    lib = load("chiral-minimal-neutral-gravity")
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    beta = [v/(2*arb.pi()) for v in model["beta"]]
    projection = arb_mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,*q],
                          [0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    coeff = projection*kernel["series"](rows,fields)
    scalar,evaluate = global_source["scalar"],global_source["evaluate"]
    # context 已证 s≥1/20 时 z/(s*x)<4。对 d/d log s 的余项也用此界：
    # 每个遗漏负幂被其幂次乘一次，三个精确核的系数仍由 1/n! 支配。
    n = kernel["ORDER"]+1
    factor = 2*arb(4)**n/(4*arb.pi()*factorial(n))
    sizes = [sum(abs(2*a) for _,a,_ in rows),sum(abs(g) for g,_,_ in rows),
             sum(abs(f*b) for bs in fields for f,b in zip(q,bs,strict=True)),
             *[sum(abs(bs[k]) for bs in fields) for k in range(3)]]
    weights = [total,*[2*abs(edges[i][2]) for i in LINES]]
    measures = [moments,*[choices[i][0] for i in LINES]]
    data = []
    for j,size in enumerate(sizes):
        row = arb_mat([[coeff[j,k] for k in range(coeff.ncols())]])
        terms = []
        for measure,weight in zip(measures,weights,strict=True):
            co = scalar(row,measure)
            powers,log_weight = co
            ds = [(2-k)*p for k,p in enumerate(powers)]
            ds[2] += log_weight
            terms.append((co,(ds,arb(0)),arb(0,(factor*size*weight).upper())))
        data.append(terms)
    target = [arb(1),arb("1e-122"),-v0]
    allowed_moves = [(-choices[i][1]).union(1-choices[i][1]) for i in LINES]

    def point(s, *, derivatives=False, physical_only=False, share_curvature=False):
        invs,logarithm = 1/s,s.log()
        vals = [[evaluate(co,err,s,invs,logarithm,0) for co,ds,err in row] for row in data]
        if physical_only:
            for row,rhs in zip(vals[:3],target,strict=True):
                span = row[0]+sum(v*x for v,x in zip(row[1:],allowed_moves,strict=True))
                if not span.overlaps(rhs):
                    return EXCLUDED
        matrix = arb_mat([a[1:] for a in vals[:3]])
        rhs = arb_mat([[b-a[0]] for b,a in zip(target,vals[:3],strict=True)])
        if matrix.det().contains(0):
            return None
        try:
            moves = matrix.solve(rhs).entries()
        except ZeroDivisionError:
            return None
        shares = [choices[i][1]+v for i,v in zip(LINES,moves,strict=True)]
        if physical_only and any(v<0 or v>1 for v in shares):
            return EXCLUDED
        couplings = [a+row[0]+sum(b*v for b,v in zip(row[1:],moves,strict=True)) for a,row in zip(initial,vals[3:],strict=True)]
        log_energy = (couplings[0]-couplings[2])/(beta[0]-beta[2])
        common = couplings[2]-beta[2]*log_energy
        floor = min(v.lower() for v in (*couplings,common))
        if physical_only and (log_energy<0 or log_energy>(arb("1.22089e19")/model["mz"]).log()
                              or any(v<0 for v in (*couplings,common))):
            return EXCLUDED
        energy = model["mz"]*log_energy.exp()
        if not derivatives and not share_curvature:
            return energy,shares,floor
        ds = [[evaluate(dco,err,s,invs,logarithm,0)*invs for co,dco,err in row] for row in data]
        at_moves = [row[0]+sum(b*v for b,v in zip(row[1:],moves,strict=True)) for row in ds]
        share_slopes = matrix.solve(arb_mat([[-v] for v in at_moves[:3]])).entries()
        if share_curvature:
            # s² f″=(D²-D)f, D=s d/ds。对遗漏的 r^-m，二阶系数最多为 (m+1)/m!。
            # z/x<4 且 m>=49，后继绝对项比小于 1/2。
            accelerations = []
            for row,first in zip(data[:3],ds[:3],strict=True):
                second = []
                for (powers,log_weight),_,err in row:
                    coeff2 = [(2-k)*(1-k)*v for k,v in enumerate(powers)]
                    coeff2[2] -= log_weight
                    second.append(evaluate((coeff2,arb(0)),(n+1)*err,s,invs,logarithm,0)*invs**2)
                fixed = second[0]+sum(a*v for a,v in zip(second[1:],moves,strict=True))
                drift = 2*sum(a*v for a,v in zip(first[1:],share_slopes,strict=True))
                accelerations.append([-fixed-drift])
            return shares,share_slopes,matrix.solve(arb_mat(accelerations)).entries()
        slopes = [a+sum(b*d for b,d in zip(row[1:],share_slopes,strict=True))
                  for a,row in zip(at_moves[3:],vals[3:],strict=True)]
        log_slope = (slopes[0]-slopes[2])/(beta[0]-beta[2])
        return energy,shares,floor,share_slopes,log_slope
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
