"""两组零点共用连续插值位置的任意逐线份额，认证规范投影的下界。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
DIRECTION = (-100000000000,-100000000000,-30037881860,87237760804,
             -71525119424,-100000000000,34465086478,-100000000000)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def projected_minimum(basis, unit, lines, tolerance):
    """区间导数有定号时直接取端点，否则用中值界或二等分。"""
    def response(t, derivative=0):
        result = arb(0)
        for left,span,coeff in lines:
            term = (coeff*basis(left+span*t,derivative))[0,0]
            result += term*(span/unit if derivative else 1)
        return result

    zero, one = arb(0), arb(1)
    f0, f1 = response(zero), response(one)
    best = min(f0.upper(),f1.upper())
    todo = [(zero,one,f0,f1)]
    lower = best
    while todo:
        lo,hi,fl,fr = todo.pop()
        mid = (lo+hi)/2
        der = response(lo.union(hi),1)
        if der >= 0:
            lower = min(lower,fl.lower())
        elif der <= 0:
            lower = min(lower,fr.lower())
        else:
            fm = response(mid)
            best = min(best,fm.upper())
            bound = (fm-abs(der)*(hi-lo)/2).lower()
            if bound >= best-tolerance:
                lower = min(lower,bound)
            else:
                assert hi-lo > arb(2)**-48
                todo.extend(((lo,mid,fl,fm),(mid,hi,fm,fr)))
    return lower


def verify():
    source = load("shared-twelve-node-gauge-gap")
    kernel = load("continuous-split-roots")
    series, basis = kernel["series"],kernel["basis"]
    _,lib,backgrounds = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,source["HEAVY"],tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0)]+load("signed-gravity-unification")["gauge_weights"](model)
    shell = load("native-rh-ray-pvm")["shell_weights"]
    scale = arb(7)/1000
    y = [arb(n)/10**11 for n in DIRECTION]
    limit = (arb("1.22089e19")/model["mz"]).log()
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    constant, error = arb(0), arb(0)
    tables = []
    for n,bg in enumerate(backgrounds):
        ky,wy,u,v = y[4*n:4*n+4]
        direction = [ky,wy,u,v,-u-v]
        coeff = arb_mat([direction])*series(bg,fields)
        edges = source["paired_lines"](ROOT,shell,source["SHELL"][n])
        assert all(0 < a and (a < b or a == b) and b <= 1 for a,b,_ in edges)
        assert len(edges) == 361
        floor = min((scale*a).lower() for a,_,_ in edges)
        zmax = max(z.upper() for _,_,z in bg)
        assert zmax < 1 and zmax/floor < 6
        order = kernel["ORDER"]
        # 三种热核的对数项均已保留。每个 n 阶幂尾由 (z/x)^n/n! 乘场系数控制。
        tail = 2*arb(6)**(order+1)/factorial(order+1)
        assert arb(6)/(order+2) < arb(1)/2
        size = sum(abs(2*ky*a)+abs(wy*g)+sum(abs(d*b) for d,b in
                   zip(direction[2:],bs,strict=True))
                   for (g,a,_),bs in zip(bg,fields,strict=True))/(4*arb.pi())
        error += tail*size*sum(abs(w) for _,_,w in edges)
        tables.append([(scale*a,scale*(b-a),w*coeff) for a,b,w in edges])
        rate = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
        constant += -ky-wy*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
        constant -= (rate*limit).max(0)
    lower = constant-error
    for lines in zip(*tables,strict=True):
        lower += projected_minimum(basis,kernel["S"],lines,arb("1e-5"))
    norm = sum(abs(y[j]) for j in (2,3,6,7))
    print("连续位置逐线投影最小和含基线及两组最坏能标，最大规范差的下界",
          lower.str(14),(lower/norm).str(14))
    assert lower/norm > arb(".8387")


if __name__ == "__main__":
    verify()
