"""低共同标尺下，两个实际背景的引力与连续位置规范交汇。"""
from functools import cache, partial
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
ORDER = 96
# 百万分之一为单位的标尺闭段，随后是投影中的四个规范系数（十亿分之一）。
BANDS = (
    (1400,1900,(-27059416,78208139,30275533,-87341509)),
    (1900,2650,(-53083065,153710094,59724967,-172688750)),
    (2650,3650,(-118644232,343990403,133897408,-387792651)),
    (3650,5000,(-212643424,616650318,240119042,-695767403)),
)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def low_newton(edges, rows, heat):
    # 对任意 0<c≤cut，x(1-E₂(z/x)) 随 x 单调不减。
    cut, pi = arb(14)/10000, arb.pi()
    slope = sum(a for _,a,_ in rows)/(2*pi)
    main = sum((slope*w*a).max(slope*w*b) for a,b,w in edges)
    deficit = sum((-w*alpha).max(0)*(cut*b)*(1-heat(2,z/(cut*b)))
                  for _,b,w in edges for _,alpha,z in rows)/(2*pi)
    bound = cut*main.max(0)+deficit
    assert bound < arb(".985")
    print("任意逐线位置下中央背景的小标尺牛顿读数上界",bound.str(12),flush=True)


def certify(left, right, gauge, data, model, initial, kernel, minimum):
    starting, ending = arb(left)/10**6, arb(right)/10**6
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    length = (arb("1.22089e19")/model["mz"]).log()
    y = [arb(n)/10**9 for n in gauge]
    norm = sum(abs(a) for a in y)
    constant, error, curvature = arb(0), arb(0), arb(0)
    tables = []
    for j,(edges,rows,fields,series) in enumerate(data):
        u,v = y[2*j:2*j+2]
        d = (-1,1,u,v,-u-v)
        coeff = arb_mat([list(d)])*series
        tables.append([(a,b-a,w*coeff) for a,b,w in edges])
        constant += 1-arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
        constant -= (length*(u*(beta[0]-beta[2])+v*(beta[1]-beta[2]))).max(0)
        limit = max(z.upper() for _,_,z in rows)/min((starting*a).lower() for a,_,_ in edges)
        assert limit < 25 and limit/(ORDER+2) < arb(1)/2 and all(z < 1 for _,_,z in rows)
        tail = 2*limit**(ORDER+1)/factorial(ORDER+1)
        heavy = [(abs(2*alpha),abs(g),abs(sum(x*b for x,b in zip(d[2:],bs,strict=True))),z)
                 for (g,alpha,z),bs in zip(rows,fields,strict=True)]
        size = sum(A+D+H for A,D,H,z in heavy)/(4*arb.pi())
        error += tail*size*sum(abs(w) for _,_,w in edges)
        for a,b,w in edges:
            # 对 log c 求两次导数，取 u=z/x，用 u exp(-u)≤min(1/e,u)。
            xlow,xhigh = starting*a,ending*b
            upper = sum(A*xhigh+2*D*xhigh**2+H*(z/xlow).min(1/arb(1).exp())
                        for A,D,H,z in heavy)
            curvature += abs(w)*upper/(4*arb.pi())
    basis = partial(kernel["basis"],order=ORDER)

    @cache
    def at(micro):
        c = arb(micro)/10**6
        value = constant-error
        for line in zip(*tables,strict=True):
            atoms = [(c*a,c*span,coeff) for a,span,coeff in line]
            value += minimum(basis,kernel["S"],atoms,arb("1e-5"))
        return value.lower()

    todo, bounds = [(left,right)], []
    while todo:
        lo,hi = todo.pop()
        lower = (min(at(lo),at(hi))-curvature*(arb(hi)/lo).log()**2/8)/norm
        if lower > arb(".05"):
            bounds.append(lower.lower())
        else:
            mid = (lo+hi)//2
            assert lo < mid < hi, (lo,hi,lower)
            todo.extend(((lo,mid),(mid,hi)))
    assert min(bounds) > arb({1400:"3.24",1900:".44",2650:".95",3650:".52"}[left])
    print("低标尺（百万分之一）、覆盖闭段数、最大规范差下界",
          left,right,len(bounds),min(bounds).str(12),flush=True)


def verify():
    source = load("shared-twelve-node-gauge-gap")
    kernel = load("continuous-split-roots")
    _,lib,backgrounds = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,source["HEAVY"],tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0)]+load("signed-gravity-unification")["gauge_weights"](model)
    shell = load("native-rh-ray-pvm")["shell_weights"]
    data = []
    for code,rows in zip(source["SHELL"],backgrounds,strict=True):
        edges = source["paired_lines"](ROOT,shell,code)
        data.append((edges,rows,fields,kernel["series"](rows,fields,order=ORDER)))
    low_newton(data[0][0],data[0][1],lib["expint"])
    minimum = load("shared-continuous-position-gap")["projected_minimum"]
    for lo,hi,coefficients in BANDS:
        certify(lo,hi,coefficients,data,model,initial,kernel,minimum)


if __name__ == "__main__":
    verify()
