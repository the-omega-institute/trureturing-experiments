"""两组零点共用连续插值位置的任意逐线份额，认证规范投影的下界。"""
from functools import cache
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
    scale_min, scale_max = arb(5)/1000, arb(9)/1000
    curvature = arb(0)
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
        floor = min((scale_min*a).lower() for a,_,_ in edges)
        zmax = max(z.upper() for _,_,z in bg)
        assert zmax < 1 and zmax/floor < 7
        order = kernel["ORDER"]
        # 三种热核的对数项均已保留。每个 n 阶幂尾由 (z/x)^n/n! 乘场系数控制。
        tail = 2*arb(7)**(order+1)/factorial(order+1)
        assert arb(7)/(order+2) < arb(1)/2
        size = sum(abs(2*ky*a)+abs(wy*g)+sum(abs(d*b) for d,b in
                   zip(direction[2:],bs,strict=True))
                   for (g,a,_),bs in zip(bg,fields,strict=True))/(4*arb.pi())
        error += tail*size*sum(abs(w) for _,_,w in edges)
        tables.append([(a,b-a,w*coeff) for a,b,w in edges])
        # s=log(c) 时，幂 r^p 的二阶导数为 p²r^p，对数项二阶导数为零。
        # 只需截断函数的曲率，真实热核的统一尾界在插值后扣除一次。
        unit = kernel["S"]
        for a,b,w in edges:
            high = scale_max*b/unit
            inverse = unit/(scale_min*a)
            bound = 4*abs(coeff[0,0])*high**2+abs(coeff[0,1])*high
            power = arb(1)
            for j in range(1,order+1):
                power *= inverse
                bound += j*j*abs(coeff[0,j+2])*power
            curvature += abs(w)*bound
        rate = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
        constant += -ky-wy*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
        constant -= (rate*limit).max(0)
    norm = sum(abs(y[j]) for j in (2,3,6,7))

    @cache
    def at(micro):
        scale = arb(micro)/10**6
        lower = constant-error
        for row in zip(*tables,strict=True):
            lines = [(scale*a,scale*span,coeff) for a,span,coeff in row]
            lower += projected_minimum(basis,kernel["S"],lines,arb("1e-5"))
        return lower.lower()

    fixed = at(7000)/norm
    print("共同标尺 0.007 的连续位置规范差下界",fixed.str(14),flush=True)
    assert fixed > arb(".8387")
    print("对数标尺截断投影的曲率上界",curvature.upper().str(12),flush=True)
    def band(left,right,minimum):
        work = [(left,right)]
        bounds = []
        while work:
            a,b = work.pop()
            width = (arb(b)/a).log()
            bound = (min(at(a),at(b))-curvature*width**2/8)/norm
            if bound > minimum:
                bounds.append(bound.lower())
            else:
                mid = (a+b)//2
                assert a < mid < b, (a,b,bound)
                work.extend(((a,mid),(mid,b)))
        return min(bounds),len(bounds)

    for left,right,goal in [(5000,6000,".05"),(6000,8000,".45"),(8000,9000,".10")]:
        bound,cells = band(left,right,arb(goal))
        print("共同标尺（百万分之一）、闭段数、规范差最大绝对值下界",
              left,right,cells,bound.str(12),flush=True)


if __name__ == "__main__":
    verify()
