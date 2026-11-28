"""实际第六背景上，检查所有普朗克有界标尺和全部逐线连续分布。"""
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
ORDER = 24
START = 450


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


@cache
def inputs():
    shell = load("native-rh-ray-pvm")["shell_weights"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](ROOT,shell,2306888100090)
    assert len(edges)==361 and all(0<a and (a<b or a==b) and b<=1 for a,b,w in edges)
    _,lib,(rows,) = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,(2345676,),tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    kernel = load("continuous-split-roots")
    matrix = kernel["series"](rows,fields,order=ORDER)
    assert all(z>=0 for g,alpha,z in rows)
    return edges,rows,lib["expint"],model,initial,beta,fields,matrix,kernel["S"]


def small_scale(data):
    edges,rows,heat,*_ = data
    pi,c,l = arb.pi(),arb(START)/10**6,arb(1)/4
    A = sum(alpha for g,alpha,z in rows)/(2*pi)
    B = -sum(g for g,alpha,z in rows)/(8*pi)
    assert A<0<B
    vertex = -A/(2*l*B)
    def quadratic(t):
        return t*(A+l*B*t)
    main = arb(0)
    for a,b,w in edges:
        maximum = (w*quadratic(a)).max(w*quadratic(b))
        if not (vertex<a or vertex>b):
            maximum = maximum.max(w*quadratic(vertex))
        main += maximum
    error_k = sum((-w*alpha).max(0)*c*b*(1-heat(2,z/(c*b)))
                  for a,b,w in edges for g,alpha,z in rows)/(2*pi)
    error_w = sum((w*g).max(0)*b*z for a,b,w in edges for g,alpha,z in rows)/(4*pi)
    bound = c*main.max(0)+error_k+l*error_w
    print("小标尺上 K+W/(4c) 的上界",bound.str(15),flush=True)
    assert bound < arb(".963884")


def certify(part):
    left,right,vector = part
    edges,rows,_,model,initial,beta,fields,matrix,unit = inputs()
    k,h,u,v = [arb(n)/10**11 for n in vector]
    d = [k,h,u,v,-u-v]
    coeff = (arb_mat([d])*matrix).entries()
    low,high = arb(left)/10**6,arb(right)/10**6
    floor = min((low*a).lower() for a,b,w in edges)
    assert max(z.upper() for g,alpha,z in rows)<min(arb(1),2*floor)
    tail = 2*arb(2)**(ORDER+1)/factorial(ORDER+1)
    size = sum(abs(2*k*alpha)+abs(h*g)+sum(abs(t*b) for t,b in zip(d[2:],bs,strict=True))
               for (g,alpha,z),bs in zip(rows,fields,strict=True))/(4*arb.pi())
    rate = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
    limit = (arb("1.22089e19")/model["mz"]).log()
    const = -k-h*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
    const -= (rate*limit).max(0)+tail*size*sum(abs(w) for a,b,w in edges)
    norm = abs(u)+abs(v)
    assert norm>0
    def polynomial(r,derivative=False):
        # 对数导数作用在 r^p 上等于乘以 p。
        inv = 1/r
        s = arb(0)
        for j in range(ORDER,0,-1):
            s = (s+((-j if derivative else 1)*coeff[2+j]))*inv
        if derivative:
            return s+2*coeff[0]*r*r+coeff[1]*r+coeff[-1]
        return s+coeff[0]*r*r+coeff[1]*r+coeff[2]+coeff[-1]*r.log()

    # 每条线的投影只相差带符号线权。一次圈出全段所有驻点，复用到每条线上。
    todo = [(floor/unit,(high/unit).upper())]
    stationary = []
    while todo:
        a,b = todo.pop()
        box = a.union(b)
        if not polynomial(box,True).contains(0):
            continue
        if b-a < arb("1e-11")*a:
            stationary.append((box,polynomial(box)))
        else:
            middle = ((a+b)/2).mid()
            todo.extend(((a,middle),(middle,b)))

    # 对 log c 线性插值时，二阶导的上界对任取的逐线概率测度一致。
    curvature = arb(0)
    for a,b,w in edges:
        r = high*b/unit
        inverse = unit/(low*a)
        loss = arb(0)
        for j in range(ORDER,0,-1):
            loss = (loss+j*j*abs(coeff[j+2]))*inverse
        curvature += abs(w)*(4*abs(coeff[0])*r*r+abs(coeff[1])*r+loss)

    @cache
    def endpoint(micro):
        c = arb(micro)/10**6
        total = const
        for a,b,w in edges:
            l,r = c*a/unit,c*b/unit
            minimum = (w*polynomial(l)).min(w*polynomial(r))
            for root,value in stationary:
                if not (root<l or root>r):
                    minimum = minimum.min(w*value)
            total += minimum
        return total.lower()

    todo = [(left,right)]
    bounds = []
    while todo:
        a,b = todo.pop()
        width = (arb(b)/a).log()
        bound = (min(endpoint(a),endpoint(b))-curvature*width**2/8)/norm
        if bound>arb(".005"):
            bounds.append(bound.lower())
        else:
            mid = (a+b)//2
            assert a<mid<b and min(endpoint(a),endpoint(b))>arb(".005")*norm, (a,b,vector,bound)
            todo.extend(((a,mid),(mid,b)))
    return min(bounds).str(20),len(bounds)


def verify():
    data = inputs()
    small_scale(data)
    parts = json.loads((ROOT/"ch39-sixth-scale-directions.json").read_text())
    assert parts[0][0]==START and parts[-1][1]==10**6
    assert all(a[1]==b[0] for a,b in zip(parts,parts[1:]))
    with ProcessPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(certify,parts))
    lower = min(arb(bound).lower() for bound,n in results)
    print("投影方向数、闭标尺段数、最大绝对规范差之下界",len(parts),sum(n for bound,n in results),lower.str(14))
    assert lower>arb(".00728")


if __name__=="__main__":
    verify()
