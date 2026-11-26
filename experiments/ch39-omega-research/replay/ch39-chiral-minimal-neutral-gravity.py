"""三枚规范候选同时匹配牛顿与微小正曲率所需的最少额外实单态。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
chiral = run_path(str(ROOT / "ch39-chiral-common-gap-crossing.py"))


def fields(layers):
    scale = arb("2.435e18")/(724**2*arb("1.22089e19"))
    # (玻色减费米偏振, 牛顿热核权, M²/M_P²)，电弱对称原点。
    rows = [(-62,arb(45+4)/6-8,arb(0))]
    for q,dim in ((1,6),(4,3),(2,3),(3,2),(6,1)):
        for h in layers[q][3:]:
            extra = 4 if q==3 else 0
            rows.append((-4*dim+extra,arb(dim)/3+arb(extra)/6,
                         ((32 if q%2==0 else 1)*scale*h)**2))
    rows.extend((36,arb(-6),(32*scale*h)**2) for h in layers[0])
    return rows


def expint(k,v):
    if v.is_zero():
        return arb(1)/(k-1)
    # 核对正参数单调递减。
    return v.upper().expint(k).union(v.lower().expint(k))


def verify_case(rows):
    pi = arb.pi()
    assert sum(g for g,_,_ in rows)==-310
    assert (sum(a for _,a,_ in rows)+arb(103)/6).contains(0)
    moment = sum(g*z for g,_,z in rows)
    absolute = sum(abs(g)*z for g,_,z in rows)
    zmax = max(z.upper() for _,_,z in rows)
    # 若最多 310 份，则任意单态质量下 K=M_P² 均要求 x 超过此粗界。
    ceiling = arb(45+4)/6-8+arb(240)/3+arb(64+310)/6
    xmin = 2*pi/ceiling
    margin = moment-2*(zmax/xmin).sqrt()*absolute
    assert margin>arb(".000044")
    bound = xmin*margin/(4*pi)
    assert bound>arb("1.5e-7")

    # 311 份取相同的 u=M_s²/Q²。归一化 K=x ΣαE₂/(2π)。
    def sums(x,u):
        heat = sum(a*expint(2,z/x) for _,a,z in rows)+arb(311)/6*expint(2,u)
        trace = sum(g*expint(3,z/x) for g,_,z in rows)+311*expint(3,u)
        return x*heat/(2*pi),-x**2*trace/(4*pi)
    u0,u1 = arb(".0015"),arb(".0018")
    ubox = u0.union(u1)
    xa,xb = arb(".17"),arb(".21")
    xbox = xa.union(xb)
    assert sums(xa,ubox)[0]<1<sums(xb,ubox)[0]
    derivative = sum(a*(-z/xbox).exp() for _,a,z in rows)+arb(311)/6*expint(2,ubox)
    assert derivative>30
    def x_for(u):
        lo,hi = xa,xb
        for _ in range(39):
            mid=(lo+hi)/2
            val=sums(mid,u)[0]
            if val<1:
                lo=mid
            elif val>1:
                hi=mid
            else:
                break
        return lo.union(hi)
    A,B = sums(x_for(u0),u0)[1],sums(x_for(u1),u1)[1]
    assert A<0<arb("1e-122")<B
    lo,hi=u0,u1
    for _ in range(31):
        mid=(lo+hi)/2
        v=sums(x_for(mid),mid)[1]
        if v<0:
            lo=mid
        elif v>arb("1e-122"):
            hi=mid
        else:
            break
    us = lo.union(hi)
    # x(u) 递增。E₂'(u)=−E₁(u)。
    xx=x_for(lo).union(x_for(hi))
    Kx = sum(a*(-z/xx).exp() for _,a,z in rows)+arb(311)/6*expint(2,us)
    dx = xx*(arb(311)/6)*expint(1,us)/Kx
    trace = sum(g*expint(3,z/xx) for g,_,z in rows)+311*expint(3,us)
    dt_dx = sum(g*z*expint(2,z/xx) for g,_,z in rows)/xx**2
    dcurv = (xx**2*311*expint(2,us) - dx*(2*xx*trace+xx**2*dt_dx))/(4*pi)
    assert dcurv>0
    delta = arb("1e-122")/dcurv
    print("310 份以下曲率下界",bound.lower().str(9),"311 份解的 u、Q/M_P、M_s/M_P、零到目标的 Δu",
          *[(v.lower().str(11),v.upper().str(11)) for v in (us,xx.sqrt(),(us*xx).sqrt(),delta)])


def verify():
    model,_,_,_,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec=192
    read=chiral["source"]["spectral_supplier"](model)
    for n in (11661,11934,12207):
        rows=read(acb.zeta_zero(n).imag)
        layers={q:chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}
        print("规范候选",n)
        verify_case(fields(layers))


if __name__=="__main__":
    verify()
