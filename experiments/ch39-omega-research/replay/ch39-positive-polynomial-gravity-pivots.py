"""两种首枢轴的非负系数截止多项式在普朗克上限下的全次数引力界。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def support(weights, ratios, tau, bound):
    # 二分端点用整数编号，一个网格单位是 1/256 个次数。
    unit, stop = 256, 1600*256

    @cache
    def powers(index):
        return [r**(arb(index)/unit) for r in ratios]

    @cache
    def moment(lo, hi):
        return sum(w*a.union(b) for w,a,b in zip(weights,powers(lo),powers(hi),strict=True))

    error = sum(abs(c)*r**1600 for c,r in zip(weights[:-1],ratios[:-1],strict=True))
    tail = weights[-1] + arb(0, error.upper())
    assert tail/2 > -bound and tail-tau*tail > -bound
    pending = [(0,stop)]
    while pending:
        a,b = pending.pop()
        m = moment(a,b)
        if m/2 > -bound and (m+tail)/2-tau*tail > -bound:
            continue
        assert b-a > 1
        mid = (a+b)//2
        pending.extend(((a,mid),(mid,b)))

    rectangles = [(0,stop,0,stop)]
    count = 0
    while rectangles:
        a,b,c,d = rectangles.pop()
        if (moment(a,b)+moment(c,d))/2 - tau*moment(a+c,b+d) > -bound:
            count += 1
            continue
        assert max(b-a,d-c) > 1
        if b-a >= d-c:
            mid = (a+b)//2
            rectangles.extend(((a,mid,c,d),(mid,b,c,d)))
        else:
            mid = (c+d)//2
            rectangles.extend(((a,b,c,mid),(a,b,mid,d)))
    return count


def verify(root):
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    _,heat,(rows,)=gravity["broken_fields"](root,[2_306_888_100_085])
    ctx.prec=256
    pi=arb.pi()
    alpha=-sum(a for g,a,z in rows)/(2*pi)
    beta=-sum(g for g,a,z in rows)/(8*pi)
    error_k=sum(abs(a)*(1-heat(2,z)) for g,a,z in rows)/(2*pi)
    error_w=sum(abs(g)*z for g,a,z in rows)/(4*pi)
    for q,levels,direct,_ in native["shell_weights"](root):
        if q not in (-4,6):
            continue
        weights=[z.real for z in direct[0]]
        ratios=[z/levels[-1] for z in levels[:-1]]+[arb(1)]
        degree,tau,bound=(108,arb(".72"),arb(".308")) if q==-4 else (62,arb(".76"),arb(".348"))
        count=support(weights,ratios,tau,bound)
        def moment(n):
            return sum(c*r**n for c,r in zip(weights,ratios,strict=True))
        tv=sum(abs(c) for c in weights)
        limit=alpha*bound+tv*(error_k+alpha*tau*error_w/beta)
        assert limit<1
        overrun = ((1-alpha*tau*tv*error_w/beta)/(alpha*bound+tv*error_k)).sqrt()
        assert overrun > 1
        print("荷、对偶矩下界、质量误差后的 K 上界",q,(-bound).str(9),limit.str(12),"连续指数矩形数",count,flush=True)
        print("两项引力同时匹配所需 Qmax/MP 下界",overrun.str(13),flush=True)
        assert moment(2*degree)<0<moment(2*degree+2)
        def zeroth(t):
            return sum(c*(r**degree*(1-t*(1-r)))**2 for c,r in zip(weights,ratios))
        lo,hi=arb(0),arb(1)
        for _ in range(44):
            mid=(lo+hi)/2
            if zeroth(mid)<0:
                lo=mid
            else:
                assert zeroth(mid)>0
                hi=mid
        t=lo.union(hi)
        kappa=-1/(alpha*sum(c*r**degree*(1-t*(1-r)) for c,r in zip(weights,ratios)))
        assert kappa>1
        print("无质量相邻幂变号、t、尺度",degree,moment(2*degree).str(12),
              moment(2*degree+2).str(12),t.str(12),kappa.str(12),flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
