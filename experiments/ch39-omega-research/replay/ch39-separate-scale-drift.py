"""两个背景沿用逐线位置分布时，要使联合读数精确所需的标尺差。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def directions(high_path):
    small = run_path(str(ROOT/"ch39-shared-small-scale-gap.py"))
    floors = ("3.24",".44",".95",".52")
    for (lo,hi,coeffs),margin in zip(small["BANDS"],floors,strict=True):
        u0,v0,u1,v1 = [arb(n)/10**9 for n in coeffs]
        yield lo,hi,[-1,1,u0,v0,-1,1,u1,v1],arb(margin)
    mid = run_path(str(ROOT/"ch39-shared-continuous-position-gap.py"))
    for lo,hi,margin in ((5000,6000,".05"),(6000,8000,".45"),(8000,9000,".10")):
        yield lo,hi,[arb(n)/10**11 for n in mid["DIRECTION"]],arb(margin)
    high = json.loads(high_path.read_text())
    keys = sorted(map(int,high))
    assert keys[0] == 9000 and keys[-1] == 1000000
    cuts = [keys[0],*[(a+b)//2 for a,b in zip(keys,keys[1:])],keys[-1]]
    for k,lo,hi in zip(keys,cuts[:-1],cuts[1:],strict=True):
        yield lo,hi,[arb(n)/10**11 for n in high[str(k)]],arb(".06354")


def rate_bound(y, rows, fields, edges, left, right):
    ky,wy,u,v = y[4:]
    terms = [(2*ky*a,-wy*g,-u*(b[0]-b[2])-v*(b[1]-b[2]),z)
             for (g,a,z),b in zip(rows,fields,strict=True)]
    A,D,H = [sum(t[j] for t in terms) for j in range(3)]
    e1 = sum(abs(a)*z for a,d,h,z in terms)
    e2 = sum(abs(d)*z for a,d,h,z in terms)
    main = H*sum(w for a,b,w in edges)
    error = arb(0)
    for a,b,w in edges:
        low,high = left*a,right*b
        # Ax+Dx² 用两端和二次 Bernstein 中间系数包住。
        coeff = [low*(A+D*low), A*(low+high)/2+D*low*high, high*(A+D*high)]
        main += w*coeff[0].union(coeff[1]).union(coeff[2])
        e3 = sum(abs(h)*(z/low).min(1) for _,_,h,z in terms)
        error += abs(w)*(e1+high*e2+e3)
    return (abs(main)+error)/(4*arb.pi())


def verify():
    shared = run_path(str(ROOT/"ch39-shared-twelve-node-gauge-gap.py"))
    def load(n):
        return run_path(str(ROOT/f"ch39-{n}.py"))
    _,lib,(rows,) = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,shared["HEAVY"][1:],tight_neutral=True)
    model,_,_,_,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0)]+load("signed-gravity-unification")["gauge_weights"](model)
    edges = shared["paired_lines"](ROOT,load("native-rh-ray-pvm")["shell_weights"],shared["SHELL"][1])
    ratio = arb(101)/100
    distances = []
    previous = 1400
    for lo,hi,y,margin in directions(ROOT/"ch39-shared-large-scale-directions.json"):
        assert lo == previous < hi
        previous = hi
        left = arb(lo)/(10**6*ratio)
        right = (arb(hi)*ratio/10**6).min(1)
        speed = rate_bound(y,rows,fields,edges,left,right)
        norm = sum(abs(y[j]) for j in (2,3,6,7))
        distances.append((margin*norm/speed.upper()).lower())
    assert previous == 1000000
    drift = min(ratio.log().lower(),*distances)
    print("两个独立标尺对数之差的绝对值下界",drift.str(12),flush=True)
    assert drift > arb(".00008385")


if __name__ == "__main__":
    verify()
