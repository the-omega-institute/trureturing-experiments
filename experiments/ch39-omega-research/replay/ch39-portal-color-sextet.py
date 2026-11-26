"""最少夸克标量伙伴的色表示限制与顶夸克质量处的六维色阈值。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
gauge = run_path(str(root / "ch39-chiral-clock-unit-readout.py"))
sm = gauge["sm"]


def color_choices():
    # SU(3) 用 (p,q) 标记。复表示以实标量实现时需加上共轭。
    choices = []
    for p in range(13):
        for q in range(p+1):
            real_dim = (p+1)*(q+1)*(p+q+2)//2 * (1 if p==q else 2)
            if real_dim <= 12:
                casimir = Fraction(p*p+q*q+p*q+3*p+3*q,3)
                choices.append((real_dim,casimir))
    assert sorted(choices)==[(1,0),(6,Fraction(4,3)),(8,3),(12,Fraction(10,3))]
    averages = [(n3,n8,(6*n3*Fraction(4,3)+8*n8*3)/12)
                for n3 in range(3) for n8 in range(2) if 6*n3+8*n8<=12]
    assert max(avg for _,_,avg in averages)==2


def mismatch():
    rel = Fraction(1,1000)
    g1, _, g3 = sm["hull_inputs"](".467453 .63811 1.0583",rel)
    ys = sm["hull_inputs"](".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101",rel)
    hyper = 3*g1**2/5
    terms = [16*g3**2+4*hyper/3]*3+[16*g3**2-2*hyper/3]*3+[6*hyper]*3
    weights = [d*y**4 for y,d in zip(ys,[12]*6+[4]*3)]
    for i in range(6):
        casimir = arb(2) if i==2 else arb(10)/3
        terms[i] -= 6*casimir*g3**2 + 2*hyper/3
    lower = sum(w*(y*y+term) for w,y,term in zip(weights,ys,terms))/(8*arb.pi()**2*sum(weights))
    assert lower>arb(".06")
    print("顶层不取六维色表示时的最小对数斜率下界",lower.lower().str(10))


def all_portals_tangent():
    # 选择对角 CKM，门户和 Yukawa 在此基底的导数也对角。
    rel = Fraction(1,1000)
    g1, _, g3 = sm["hull_inputs"](".467453 .63811 1.0583",rel)
    ys = sm["hull_inputs"](".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101",rel)
    p = [y*y for y in ys]
    n = [12]*6+[4]*3
    hyper = 3*g1**2/5
    gauge_cost = [16*g3**2+4*hyper/3]*3+[16*g3**2-2*hyper/3]*3+[6*hyper]*3
    gauge_cost[2] -= 20*g3**2+2*hyper/3
    partners = p[3:6]+p[:3]+[arb(0)]*3
    eps = arb(".001")
    trace = sum(dim*u for dim,u in zip(n,p))
    # V_S=[(Σb_i S_i²)²+eps(ΣS_i²)²]/8，各块的 S_i² 为范数平方。
    required = [-u*(u+3*v+12*arb(".125")+cost)-eps*(trace+2*u)
                for u,v,cost in zip(p,partners,gauge_cost)]
    required[2] -= hyper**2/3
    def image(A):
        b = [2*r/(A+(A**2+8*u*r).sqrt()) for u,r in zip(p,required)]
        assert all(v.is_finite() for v in b)
        return sum(dim*u*v for dim,u,v in zip(n,p,b))
    lo, hi = arb(2), arb(5)
    assert all(lo**2+8*u*r>0 for u,r in zip(p,required))
    assert image(lo)>lo and image(hi)<hi
    for _ in range(16):
        mid = (lo+hi)/2
        value = image(mid)
        if value>mid:
            lo = mid
        elif value<mid:
            hi = mid
        else:
            break
    print("严格正四次势实现九层门户同时相切的 A 区间",lo.str(7),hi.str(7))


def sextet():
    # 取忠实商群允许的最小绝对超荷 |Y|=1/3，弱单态。
    indices = [arb(2)/5,arb(0),arb(5)/2]
    casimirs = [arb(1)/15,arb(0),arb(10)/3]
    b = [v/3 for v in indices]
    B = [[v*(4*casimirs[j]+(arb(2)*(0,2,3)[i]/3 if i==j else 0))
          for j in range(3)] for i,v in enumerate(indices)]
    assert (B[2][2]-arb(115)/3).contains(0)
    return arb("172.5"),b,B


def verify():
    chiral = gauge["chiral"]
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    color_choices()
    mismatch()
    all_portals_tangent()
    read = chiral["source"]["spectral_supplier"](model)
    for n,seed_L,seed_alpha in ((22027,"29.197",".08187"),(28325,"29.198",".08185")):
        rows = read(acb.zeta_zero(n).imag)
        layers = {q:chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}
        print("正零点",n)
        mu,_,alpha = gauge["crossing_enclosure"](model,initial,layers,False,724,seed_L,seed_alpha,
                                                extra_thresholds=(sextet(),))
        assert alpha < arb(".1")


if __name__=="__main__":
    verify()
