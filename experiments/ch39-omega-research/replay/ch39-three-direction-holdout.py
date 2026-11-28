"""三组预登记方向对完整正负荷端点谱的必要条件及独立留出。"""
from concurrent.futures import ProcessPoolExecutor
from functools import cache,partial
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb,arb_mat

ROOT = Path(__file__).resolve().parent
DIRECTIONS = ((-5102579103,-100000000000,120238954,-346875920),
              (-4915099059,-100000000000,286857007,-830616459),
              (-4993520318,-100000000000,83442995,-240153821))


@cache
def sources():
    return run_path(str(ROOT/"ch39-real-minima-window.py"))


@cache
def coefficients(vector):
    _,_,rows,model,initial,beta,fields,kernel,series = sources()["inputs"]()
    k,h,u,v = [arb(n)/10**11 for n in vector]
    y = [k,h,u,v,-u-v]
    c = (arb_mat([y])*series).entries()
    rate = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
    base = (-k-h*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
            -rate*(arb(10**12)/model["mz"]).log())
    size = sum(abs(2*k*a)+abs(h*g)+sum(abs(d*b) for d,b in zip(y[2:],bs,strict=True))
               for (g,a,z),bs in zip(rows,fields,strict=True))/(4*arb.pi())
    return c,base,size*2*arb(9)**49/factorial(49),abs(u)+abs(v)


def lower_score(edges,vector,kernel):
    c,total,tail,norm = coefficients(vector)
    basis = partial(kernel["basis"],order=2)
    short = arb_mat([[*c[:5],c[-1]]])
    minimum = sources()["load"]("shared-continuous-position-gap")["projected_minimum"]
    for a,b,w in edges:
        inv = 200*kernel["S"]/a
        remainder = arb(0)
        for j in range(48,2,-1):
            remainder = (remainder+abs(c[j+2]))*inv
        total -= abs(w)*(tail+remainder*inv**2)
        total += minimum(basis,kernel["S"],[(a/200,(b-a)/200,w*short)],arb("1e-8"))
    return total.lower()/norm


def evaluate(index):
    shell = sources()
    edges,_,_ = shell["shell_data"](index)
    _,_,rows,_,_,_,_,kernel,_ = shell["inputs"]()
    assert max(z.upper() for g,a,z in rows)<min(arb(1),9*min(a.lower()/200 for a,b,w in edges))
    vectors = DIRECTIONS
    if index<128:
        vectors = (*vectors,shell["load"]("four-geometric-scores")["DIRECTION"])
    scores = [lower_score(edges,v,kernel) for v in vectors]
    return [s.str(18) for s in scores]


def verify():
    original = json.loads((ROOT/"ch39-real-minima-window.json").read_text())
    extension = json.loads((ROOT/"ch39-real-minima-window-extension.json").read_text())
    holdout = json.loads((ROOT/"ch39-three-direction-holdout-layouts.json").read_text())
    known = {**original,**extension}
    assert sorted(map(int,known))==list(range(128)) and sorted(map(int,holdout))==list(range(128,192))
    with ProcessPoolExecutor(max_workers=6) as pool:
        bounds = [[arb(v) for v in row] for row in pool.map(evaluate,range(192))]
        # 预报与各壳独立构造或连续域分离证书逐一比对。
        items = [(i,holdout[str(i)]) for i in range(128,192)]
        judged = list(pool.map(independent_classification,items))
    new = {i:positive for i,positive,s,t,lower in judged}
    for start,stop in [(0,128),(128,192)]:
        confusion = {(p,t):0 for p in (True,False) for t in (True,False)}
        margins,misses = [],[]
        for i in range(start,stop):
            bound = max(s.lower() for s in bounds[i][:3])
            predicted = bound<=0
            actual = ("base" in known[str(i)]) if i<128 else new[i]
            confusion[predicted,actual]+=1
            if predicted!=actual:
                misses.append(i)
            if not predicted:
                margins.append(bound)
        print("起止偏移、预测与真实计数、最小排除量、错分",start,stop,confusion,min(margins).str(12),misses)
        if start==0:
            assert confusion=={(True,True):28,(False,False):100,(False,True):0,(True,False):0}
            assert min(margins)>arb(".06725")
        else:
            assert confusion=={(True,True):17,(False,False):46,(False,True):0,(True,False):1}
            assert misses==[162] and min(margins)>arb("1.1314")
    reference_count=sum(row[3]>0 for row in bounds[:128])
    assert reference_count==86
    print("旧方向在前 128 枚中认证的排除数",reference_count)


def independent_classification(item):
    return sources()["classify"](item)


if __name__=="__main__":
    verify()
