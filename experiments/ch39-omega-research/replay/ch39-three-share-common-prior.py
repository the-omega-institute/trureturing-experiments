"""同一三线族采用公共伯努利参考概率时，能标的全域选择与有限计数门槛。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def logit(p):
    return (p/(1-p)).log()


def weighted_logit(p):
    if not p.contains(0):
        return p*logit(p)
    hi = p.upper()
    assert hi < arb(".1")
    if hi.is_zero():
        return arb(0)
    # x log(x) 在 [0,0.1] 上递减，端点取连续延拓。
    return (hi*hi.log()).union(0)-p*(1-p).log()


def bisect(lo,hi,increasing_residual):
    assert increasing_residual(lo)<0<increasing_residual(hi)
    for _ in range(50):
        mid = ((lo+hi)/2).mid()
        v = increasing_residual(mid)
        if v<0:
            lo = mid
        elif v>0:
            hi = mid
        else:
            break
    return lo.union(hi)


def verify():
    global_energy = run_path(str(ROOT/"ch39-three-share-global-energy.py"))
    global_energy["verify"]()
    point,mz = run_path(str(ROOT/"ch39-three-share-energy-interval.py"))["setup"]()

    def stationary_prior(s):
        _,p,_,dp,_ = point(s,derivatives=True)
        b = sum(dp)
        assert b>0 and all(0<v<1 for v in p)
        z = sum(d*logit(v) for d,v in zip(dp,p,strict=True))/b
        return 1/(1+(-z).exp())

    def monotone(s):
        result = point(s,share_curvature=True)
        if result is None:
            return False
        p,dp,ddp = result
        if p[1]<0:
            return True
        if not (0<p[0]<1 and 0<p[2]<1 and p[1]<1 and dp[1]>0):
            return False
        x = p[1].intersection(arb(0).union(1))
        b,db = sum(dp),sum(ddp)
        if not b>0 or (x.contains(0) and not x.upper()<arb(".1")):
            return False
        derivative = b*dp[1]**2/(1-x)+(b*ddp[1]-db*dp[1])*weighted_logit(x)
        for k in (0,2):
            derivative += x*(b*dp[k]**2/(p[k]*(1-p[k]))+(b*ddp[k]-db*dp[k])*logit(p[k]))
        return derivative>0

    boxes = global_energy["cover"](((arb("1.00060376910"),arb("1.00417102488")),),monotone)
    print("驻点参考概率严格递增的覆盖盒数",boxes)
    for energy,low,high in ((mz,".55366743",".55366747"),(arb(10)**13,".18748425",".18748428"),
                            (arb(10)**15,".06344454",".06344458")):
        def residual(s):
            return energy-point(s)[0]
        s = bisect(arb("1.00060376912"),arb("1.00418"),residual)
        q = stationary_prior(s)
        assert arb(low)<q<arb(high)
        print("能标门槛、参考概率",energy.str(9),q.str(12),flush=True)
    cases = ((1,6,"2.58e13","2.59e13"), (1,16,"1.0248e15","1.0249e15"),
             (1,19,"1.30954e15","1.30955e15"), (2,19,"2.83940e14","2.83942e14"),
             (3,19,"3.77328e13","3.77331e13"), (4,19,"3.24696e12","3.24699e12"))
    for k,d,lo,hi in cases:
        q = arb(k)/d
        s = bisect(arb("1.00060380"),arb("1.00417102"),lambda s: stationary_prior(s)-q)
        energy,_,floor = point(s)
        assert arb(lo)<energy<arb(hi) and floor>0
        print("有限计数及唯一能标",k,d,energy.str(12),flush=True)


if __name__ == "__main__":
    verify()
