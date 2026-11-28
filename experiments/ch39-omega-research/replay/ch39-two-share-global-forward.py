"""第16、301线迁移二十个前向重谱时，普朗克有界标尺上的全部物理解。"""
from functools import cache
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


@cache
def context():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,choices = repair["source_data"]()
    targets = tuple(range(2345988,2346027,2))
    _,lib,bgs = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(repair["SOURCE"],*targets),tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    beta = [v/(2*arb.pi()) for v in model["beta"]]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    q = [1-t,t,-1]
    v0 = sum(a*d for a,d in zip(initial,q,strict=True))
    total = sum(abs(w) for a,b,w in edges)
    # 对小标尺直接将各线实际位置一次矩与无质量热核误差分开。
    scale = arb(7)/1000
    low = arb(1)/10
    src = bgs[0]
    loss_rows = [(g,a,max(bg[i][2].upper() for bg in bgs[1:])) for i,(g,a,z) in enumerate(src)]
    def loss(rows,s):
        return sum(abs(w)*sum(abs(h)*(scale*s*b)*(1-lib["expint"](2,z/(scale*s*b)))
                   for g,h,z in rows) for _,b,w in edges)/(2*arb.pi())
    main = sum(a for g,a,z in src)/(2*arb.pi())
    maximum = 1+loss(src,1)
    for line in (16,301):
        delta,f = choices[line]
        maximum += ((-f).union(1-f)*main*delta[1,0]*kernel["S"]).upper()
    assert maximum>0
    small = low*maximum+loss(loss_rows,low)
    print("标尺比不超过 0.1 时的牛顿矩上界",small.str(12),flush=True)
    assert small<arb(".949705")
    assert all(max(z.upper() for g,a,z in rows)<min(arb(1),2*low*scale*min(a.lower() for a,b,w in edges))
               for rows in bgs)
    return kernel,moments,choices,edges,fields,dict(zip(targets,bgs[1:],strict=True)),q,v0,total


def scalar(row, measure):
    """固定谱矩向量下，预乘幂项和 log(s) 的系数。"""
    entries = row.entries()
    powers = [coefficient*measure[j,0] for j,coefficient in enumerate(entries[:-1])]
    powers[2] += entries[-1]*measure[measure.nrows()-1,0]
    return powers, entries[-1]*measure[2,0]


def evaluate(coefficients, error, s, invs, logarithm, degree):
    """在标尺区间上包住单个读数除以 s^degree 的值。"""
    powers, log_weight = coefficients
    tail = arb(0)
    for coefficient in reversed(powers[3:]):
        tail = (tail+coefficient)*invs
    constant = powers[2]+log_weight*logarithm+tail+error
    return powers[0]*s**(2-degree)+powers[1]*s**(1-degree)+constant*invs**degree


def global_domain(number,root,physical,read):
    kernel,moments,choices,edges,fields,bgs,q,v0,total = context()
    rows = bgs[number]
    series = kernel["series"](rows,fields)
    mat = arb_mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,*q]])*series
    mi,pi = choices[16]
    mj,pj = choices[301]
    physical_d = (-pi).union(1-pi)
    physical_e = (-pj).union(1-pj)
    rows_by_observable = [arb_mat([[mat[k,j] for j in range(mat.ncols())]]) for k in range(3)]
    tab = [[scalar(row,measure) for measure in (moments,mi,mj)] for row in rows_by_observable]
    # 任意标尺 s≥0.1 时 z/x<2。基准矩为原概率分布，差矩为两原子之差。
    n = kernel["ORDER"]+1
    factor = 2*arb(2)**n/(4*arb.pi()*factorial(n))
    sizes = [sum(abs(2*a) for g,a,z in rows),sum(abs(g) for g,a,z in rows),
             sum(sum(abs(t*b) for t,b in zip(q,bs,strict=True)) for bs in fields)]
    errors = [[arb(0,(factor*size*w).upper()) for w in (total,2*abs(edges[16][2]),2*abs(edges[301][2]))]
              for size in sizes]
    # 用同一三元方程认证较宽盒内至多一个数学根。
    radius = arb(".0001")
    domain = [center+arb(0,radius) for center in root]
    _,midjac = read(*root)
    _,interval = read(*domain)
    inv = arb_mat(3,3,[z.mid() for z in midjac.inv().entries()])
    assert not inv.det().contains(0)
    C = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])-inv*interval
    assert max(sum(abs(C[i,j]) for j in range(3)) for i in range(3))<1
    todo = [(arb(1)/10,arb(1000)/7)]
    visited = 0
    while todo:
        lo,hi = todo.pop()
        s = lo.union(hi)
        assert s>0
        invs = 1/s
        logarithm = s.log()
        # 除以相应标尺幂后，依次求 K/s、W/s² 与 D 的常数项和两个份额系数。
        out = [[evaluate(coeff,err,s,invs,logarithm,degree)
                for coeff,err in zip(columns,remainder,strict=True)]
               for columns,remainder,degree in zip(tab,errors,(1,2,0),strict=True)]
        k,a,b = out[0]
        w,c,d = out[1]
        invariant,u,v = out[2]
        det = a*d-b*c
        visited+=1
        if not det.contains(0):
            rhs1 = invs-k
            rhs2 = arb("1e-122")*invs**2-w
            x = (rhs1*d-b*rhs2)/det
            y = (a*rhs2-rhs1*c)/det
            if not physical_d.overlaps(x) or not physical_e.overlaps(y):
                continue
            x = x.intersection(physical_d)
            y = y.intersection(physical_e)
            if not (v0+invariant+u*x+v*y).contains(0):
                continue
            if all(dom.contains(val) for dom,val in zip(domain,(s,x,y),strict=True)):
                continue
        mid = ((lo+hi)/2).mid()
        assert hi-lo>arb("1e-11"),(number,s,det)
        todo.extend(((lo,mid),(mid,hi)))
    print("目标、全域物理解数、已检验标尺盒数",number,int(physical),visited,flush=True)


if __name__=="__main__":
    replay = run_path(str(ROOT/"ch39-two-source-share-crossings.py"))
    replay["verify"](forward=True,check_domain=global_domain)
