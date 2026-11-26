"""对 11410、11431、11432 的三段作用，用整数时钟同时区分规范与端口读数。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root = Path(__file__).parent
support = run_path(str(root / "ch39-zero-port-support-switches.py"))
ordered = run_path(str(root / "ch39-ordered-flavor-unitary.py"))
study,box,edges = (support[k] for k in ("study","box","edges"))


def least_clock(lower,upper):
    center = Fraction(str(float((lower+upper)/2)))
    # 近似只用于找候选；最短性再用确切相邻分数证明。
    def inside(f):
        return lower < arb(f.numerator)/f.denominator < upper
    lo,hi = 1,1
    while not inside(center.limit_denominator(hi)):
        hi *= 2
    while lo+1<hi:
        mid=(lo+hi)//2
        if inside(center.limit_denominator(mid)):
            hi=mid
        else:
            lo=mid
    best=center.limit_denominator(hi)
    p,q=best.numerator,best.denominator
    b=pow(p,-1,q)
    a=(p*b-1)//q
    c,d=p-a,q-b
    assert b>0 and d>0 and c*b-a*d==1
    assert arb(a)/b < lower < arb(p)/q < upper < arb(c)/d
    return a,b,p,q,c,d


def binary_clock(lower,upper):
    for bits in range(1,80):
        denominator=2**bits
        below=(lower*denominator).floor().unique_fmpz()
        assert below is not None
        numerator=int(below)+1
        f=arb(numerator)/denominator
        assert 0 not in upper-f
        if f<upper:
            break
    else:
        raise AssertionError("binary clock unresolved")
    assert lower<f<upper
    return bits,numerator,denominator


def verify():
    cost,strength,scale,*_ = support["prev"]["setup"]()
    atoms = {}
    for n in (11431,11432,11415,11430,11372,11304,11410):
        height = acb.zeta_zero(n).imag
        x = arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
        c = arb(0) if n==11431 else cost(x)[0]
        t = study["native"]["plane_response"](height)
        atoms[n] = x,c,t
    r = arb(1)/3
    # 原非相干模型的平分误差边界，归一读数为 X 与 Z=X+sC。
    corner = [v+sgn*arb(2)**-28 for v,sgn in zip(box["centers"],(1,-1,1,-1))]
    coords,costs,rows = [],[],[]
    for n in (11431,11432,11415,11430,11372,11304):
        x,c,t = atoms[n]
        coords.append(x+scale*c/2)
        costs.append(c)
        rows.append(box["flavor_data"](t,corner)[0])
    ws,_ = edges["candidate"](coords,rows,r)
    limit = scale*support["dot"](ws,costs)/2
    left,cl,_ = atoms[11410]
    right,cr,_ = atoms[11432]
    lows,highs = [],[]
    for a,b in ((left,right),(left+scale*cl,right+scale*cr)):
        assert a>b
        lows.append((-limit/r-b)/(a-b))
        highs.append((limit/r-b)/(a-b))
    lower,upper = max(lows),min(highs)
    assert lower<upper and 0<lower and upper<1
    print("Budget per normalized reading",limit.str(19),"admissible external left share",
          lower.str(20),upper.str(20),flush=True)
    a,b,p,q,c,d=least_clock(lower,upper)
    print("Farey parents and shortest",a,b,p,q,c,d)
    bits,numerator,denominator=binary_clock(lower,upper)
    f=arb(numerator)/denominator
    print("First binary external clock",bits,numerator,denominator)
    cube=[v+arb(0,arb(2)**-28) for v in box["centers"]]
    # 在整个可分辨开区间闭包上检查相同乘积的四项风味读数。
    for label,f in (("whole interval",lower.union(upper)),("shortest",arb(p)/q),("binary",f)):
        weights={11410:r*f,11431:1-r,11432:r*(1-f)}
        steps={n:ordered["evolution"](atoms[n][2],w,cube) for n,w in weights.items()}
        obs=ordered["readings"](steps,(11432,11431,11410))
        assert all(abs(o-v)<tol for o,v,tol in zip(obs,ordered["flavor"]["shell"]["centers"],
                                                  [arb(".000055")]*3+[arb(".006")]))
        if label!="whole interval":
            x=r*(right+f*(left-right))
            z=x+scale*r*(cr+f*(cl-cr))
            assert max(abs(x),abs(z))<limit
            print(label,"relative gauge shift",x.str(18),"absolute port error",
                  (abs(z)*strength/scale).str(18))


if __name__=="__main__":
    verify()
