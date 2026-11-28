"""预登记的连续三十二枚壳零点，逐一验证两锚点预测与最小异类间隔。"""
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from functools import cache
from io import StringIO
import json
from pathlib import Path
from runpy import run_path

from flint import acb, arb

ROOT = Path(__file__).resolve().parent
START,WEIGHTS,HEAVY = 2306888100085,2306888100089,2345676


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


@cache
def inputs():
    shell = load("native-rh-ray-pvm")["shell_weights"]
    weight = [w for a,b,w in load("shared-twelve-node-gauge-gap")["paired_lines"](ROOT,shell,WEIGHTS)]
    _,rows,_,model,initial,beta,fields,_,_ = load("sixth-all-scales")["inputs"]()
    kernel = load("continuous-split-roots")
    series = kernel["series"](rows,fields)
    sectors = load("three-ray-chiral-boundary")["laplacian"]
    return weight,sectors,rows,model,initial,beta,fields,kernel,series


def classify(item):
    index,certificate = item
    weight,laplacian,rows,model,initial,beta,fields,kernel,series = inputs()
    height = acb.zeta_zero(START+index).imag
    spectra = []
    for charge in (-4,2):
        eig = laplacian(charge,height).eig(multiple=True)
        assert all(z.imag.contains(0) for z in eig)
        levels = sorted((z.real for z in eig),key=float)
        assert len(levels)==19 and levels[0]>0 and all(a<b for a,b in zip(levels,levels[1:]))
        spectra.append([v/levels[-1] for v in levels[:-1]]+[arb(1)])
    minus,plus = spectra
    edges = [(s*t,(s+t)/2,weight[19*i+j]) for i,s in enumerate(minus) for j,t in enumerate(plus)]
    positive = "base" in certificate
    lower = None
    if positive:
        labels = [arb(k)/11000 for row in certificate["base"] for k in row]
        extra = [(None,i,arb(k)/11000) for i,k in certificate["extra"]]
        with redirect_stdout(StringIO()):
            load("one-background-specific-line")["solve_layout"](
                [labels],extra,set(),["11.70"],[START+index],[HEAVY],False,reference_gev=10**12,
                scale_millionths=[5000],edge_tables=[edges])
    else:
        direction = certificate["direction"]+[0]
        bound,_ = load("three-background-calibration")["projection"](
            [(edges,series,rows)],fields,beta,initial,model,kernel,[direction],10**12,set(),
            scale_millionths=[5000])
        norm = sum(abs(arb(v)/10**11) for v in direction[2:4])
        assert bound/norm>0,(index,bound)
        lower = (bound/norm).str(25)
    return index,positive,minus[0].str(30),plus[0].str(30),lower


def verify():
    cert = json.loads((ROOT/"ch39-real-minima-window.json").read_text())
    assert sorted(map(int,cert))==list(range(32))
    with ProcessPoolExecutor(max_workers=6) as pool:
        observations = list(pool.map(classify,[(i,cert[str(i)]) for i in range(32)]))
    vectors = [(arb(s),arb(t)) for i,positive,s,t,lower in observations]
    def distance(i,j):
        a,b = vectors[i]
        c,d = vectors[j]
        return abs(a-c).max(abs(b-d))
    positives = [i for i,positive,*_ in observations if positive]
    negatives = [i for i,positive,*_ in observations if not positive]
    assert positives==[*range(5),*range(25,32)]
    floor = min(arb(lower).lower() for i,positive,s,t,lower in observations if not positive)
    assert floor>arb(".3433")
    print("可行偏移、排除的最小规范差下界",positives,floor.str(14))
    correct = []
    for i in range(6,32):
        a,b = distance(i,4),distance(i,5)
        assert a<b or b<a
        if (a<b)==(i in positives):
            correct.append(i)
    print("预测样本数与命中数",26,len(correct))
    assert len(correct)==26
    for start,expected in [(0,(4,24)),(6,(25,24))]:
        pairs = [(i,j) for i in positives for j in negatives if i>=start and j>=start]
        value = distance(*expected)
        assert expected in pairs and all(value<distance(*p) for p in pairs if p!=expected)
        print("异类最近两枚、距离",*(START+i for i in expected),value.str(17))
    print("锚点间距",distance(4,5).str(16))


if __name__=="__main__":
    verify()
