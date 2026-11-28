"""固定重谱与第五壳权，四种扇区组合的同一预计算投影。"""
from functools import partial
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
DIRECTION = (-5166543295,-100000000000,541403203,-1568992051)
FIFTH,SIXTH = 2306888100089,2306888100090


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def score(edges,coeff,c0,tail0,kernel,norm):
    basis = partial(kernel["basis"],order=2)
    short = arb_mat([[*coeff[:5],coeff[-1]]])
    minimum = load("shared-continuous-position-gap")["projected_minimum"]
    scale,unit = arb(1)/200,kernel["S"]
    error = arb(0)
    low,high = c0,c0
    for a,b,w in edges:
        ratio = unit/(scale*a)
        power_error = arb(0)
        for j in range(48,2,-1):
            power_error = (power_error+abs(coeff[j+2]))*ratio
        error += abs(w)*(tail0+power_error*ratio**2)
        weighted = w*short
        low += minimum(basis,unit,[(scale*a,scale*(b-a),weighted)],arb("1e-8"))
        high += min((weighted*basis(scale*(a+(b-a)*i/64)))[0,0].upper() for i in range(65))
    return (low-error).lower()/norm,(high+error).upper()/norm,error.upper()


def verify():
    _,rows,_,model,initial,beta,fields,_,_ = load("sixth-all-scales")["inputs"]()
    kernel = load("continuous-split-roots")
    k,h,u,v = [arb(n)/10**11 for n in DIRECTION]
    direction = [k,h,u,v,-u-v]
    coeff = (arb_mat([direction])*kernel["series"](rows,fields)).entries()
    size = sum(abs(2*k*a)+abs(h*g)+sum(abs(y*b) for y,b in zip(direction[2:],bs,strict=True))
               for (g,a,z),bs in zip(rows,fields,strict=True))/(4*arb.pi())
    tail0 = size*2*arb(2)**49/factorial(49)
    rate = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
    c0 = -k-h*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
    c0 -= (rate*(arb("1.22089e19")/model["mz"]).log()).max(0)
    paired = load("shared-twelve-node-gauge-gap")["paired_lines"]
    shell = load("native-rh-ray-pvm")["shell_weights"]
    for minus in (FIFTH,SIXTH):
        for plus in (FIFTH,SIXTH):
            edges = paired(ROOT,shell,minus,positive_code=plus,weight_code=FIFTH)
            assert max(z.upper() for g,a,z in rows)<min(arb(1),2*min(a.lower()/200 for a,b,w in edges))
            lo,hi,err = score(edges,coeff,c0,tail0,kernel,abs(u)+abs(v))
            print("两荷来源、归一化投影下确界的包络、未除规范系数的尾界",
                  minus,plus,lo.str(12),hi.str(12),err.str(12),flush=True)
            assert lo>0 if minus==plus==SIXTH else hi<0
            assert err<arb(".000261")
    counterexample(edges)


def counterexample(native):
    # 以最大归一化本征值 1 为搭档，乘积端点恰等于另一扇区的本征值。
    minus = [native[19*i+18][0] for i in range(19)]
    plus = [native[342+j][0] for j in range(19)]
    minimum = minus[0]
    assert minus[11]>minimum and minus[-1]==plus[-1]==1
    ratios = [(arb(i)/11)**20 for i in range(1,11)]
    assert 0<ratios[0] and ratios[-1]<1 and all(a<b for a,b in zip(ratios,ratios[1:]))
    minus[1:11] = [minimum+(minus[11]-minimum)*r for r in ratios]
    edges = [(s*t,(s+t)/2,native[19*i+j][2]) for i,s in enumerate(minus) for j,t in enumerate(plus)]
    layout = json.loads((ROOT/"ch39-same-minima-layout.json").read_text())
    labels = [arb(k)/11000 for row in layout["base"] for k in row]
    extra = [(None,i,arb(k)/11000) for i,k in layout["extra"]]
    print("保持的两荷最低本征值及其乘积",minimum.str(12),plus[0].str(12),(minimum*plus[0]).str(12))
    load("one-background-specific-line")["solve_layout"](
        [labels],extra,set(),["11.73"],[SIXTH],[2345676],False,reference_gev=10**12,
        scale_millionths=[5000],edge_tables=[edges])


if __name__=="__main__":
    verify()
