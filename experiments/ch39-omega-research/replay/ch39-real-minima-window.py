"""两个预登记窗口的连续壳零点，逐一验证两锚点预测与最小异类间隔。"""
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


def shell_data(index):
    weight,laplacian,*_ = inputs()
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
    return edges,minus[0],plus[0]


def classify(item):
    index,certificate = item
    _,_,rows,model,initial,beta,fields,kernel,series = inputs()
    edges,x,y = shell_data(index)
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
    return index,positive,x.str(30),y.str(30),lower


def verify():
    cert = json.loads((ROOT/"ch39-real-minima-window.json").read_text())
    assert sorted(map(int,cert))==list(range(32))
    extension = json.loads((ROOT/"ch39-real-minima-window-extension.json").read_text())
    assert sorted(map(int,extension))==list(range(32,128))
    cert.update(extension)
    with ProcessPoolExecutor(max_workers=6) as pool:
        observations = list(pool.map(classify,[(i,cert[str(i)]) for i in range(128)]))
    vectors = [(arb(s),arb(t)) for i,positive,s,t,lower in observations]
    def distance(i,j):
        a,b = vectors[i]
        c,d = vectors[j]
        return abs(a-c).max(abs(b-d))
    positives = [i for i,positive,*_ in observations if positive]
    negatives = [i for i,positive,*_ in observations if not positive]
    first_positives = [i for i in positives if i<32]
    assert first_positives==[*range(5),*range(25,32)]
    first_floor = min(arb(lower).lower() for i,positive,s,t,lower in observations if not positive and i<32)
    assert first_floor>arb(".3433")
    floor = min(arb(lower).lower() for i,positive,s,t,lower in observations if not positive)
    assert floor>arb(".06387")
    print("全部可行偏移、排除的最小规范差下界",positives,floor.str(14))
    assert len(positives)==28
    def confusion(start,stop):
        result = {(p,r):0 for p in (False,True) for r in (False,True)}
        for i in range(start,stop):
            a,b = distance(i,4),distance(i,5)
            assert a<b or b<a
            result[(a<b,i in positives)] += 1
        return result
    first,second = confusion(6,32),confusion(32,128)
    assert first=={(True,True):7,(False,True):0,(True,False):0,(False,False):19}
    assert second=={(True,True):16,(False,True):0,(True,False):25,(False,False):55}
    print("预登记预测与实测的计数",first,second)
    cases = [(0,32,(4,24)),(6,32,(25,24)),(32,128,(57,58)),(0,128,(57,58))]
    for start,stop,expected in cases:
        pairs = [(i,j) for i in positives for j in negatives if start<=i<stop and start<=j<stop]
        value = distance(*expected)
        assert expected in pairs and all(value<distance(*p) for p in pairs if p!=expected)
        print("窗口、异类最近两枚、距离",start,stop,*(START+i for i in expected),value.str(17))
    print("锚点间距",distance(4,5).str(16))
    print("最近对的最低值、排除壳规范差下界",vectors[57],vectors[58],observations[58][-1])


if __name__=="__main__":
    verify()
