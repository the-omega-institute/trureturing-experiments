"""由 1278 组有限时钟码约束未分辨的 CKM 读数所需的 11431 占比。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq

finite = run_path(str(Path(__file__).with_name("ch39-rational-local-flavor.py")))
anchor = tuple(fmpq(v,10**12) for v in (50610287000,1745347280,13900795,50550619000))
s12, s23, s13 = (fmpq(22501,10**5),fmpq(4183,10**5),fmpq(3732,10**6))
half = fmpq(55,10**6)
coeff = (fmpq(370538746, 100000), fmpq(2018931693, 100000), fmpq(78377747019, 100000))
tan_lower = fmpq(218162207,10**8)


def points():
    for n, codes in finite["records"]():
        w = finite["flavor_operation"](codes)
        def probability(i,j):
            return w[i,j]**2+w[i+3,j]**2
        _, j, _ = finite["measurements"](w)
        yield n, [probability(0,1),probability(1,2),probability(0,2),probability(1,0),j]


def linear_real(z):
    a,b,u,d = anchor
    v = (1+u)/(1-u)
    base = (1-u)*d-a-u*b+v*a*b
    gradient = [-1+v*b,-u+v*a,-d-b+2*a*b/(1-u)**2,1-u]
    return (base+sum(g*(x-y) for g,x,y in zip(gradient,z,anchor)))/2


def verify():
    ctx.prec = 160
    assert arb(tan_lower)<arb("1.141").tan()
    weights = {11431:470383,11831:123998,11032:98830,11432:306789}
    assert sum(weights.values())==10**6
    means = [fmpq(0) for _ in range(5)]
    minimum = None
    for n,z in points():
        a,b,u,d,j = z
        assert 0 < j < fmpq(57,10**6)
        slack = (int(n==11431) - coeff[0]*(a-(s12-half)**2*(1-u))
                 - coeff[1]*(b-(s23-half)**2*(1-u)) - coeff[2]*(j-tan_lower*linear_real(z)))
        minimum = slack if minimum is None else min(slack,minimum)
        if n in weights:
            for i,x in enumerate(z):
                means[i] += fmpq(weights[n],10**6)*x
    assert minimum > fmpq(467,1000)
    def squared_box(center):
        return arb((center-half)**2).union(arb((center+half)**2))
    u = squared_box(s13)
    a,b = [squared_box(s)*(1-u) for s in (s12,s23)]
    v = (1+u)/(1-u)
    d = (2*arb(0).union(arb("0.000057"))/arb(tan_lower).union(arb("2.26"))
         +a+u*b-v*a*b)/(1-u)
    da,db,du,dd = [abs(x-arb(y)).upper() for x,y in zip((a,b,u,d),anchor)]
    a0,b0,u0,_ = map(arb,anchor)
    remainder = (du*(dd+db) + (1+u0)/(1-u0)*da*db
                 + 2*du/(1-u.upper())**2*(abs(a0)*db+abs(b0)*da+da*db)
                 + 2*du**2/(1-u.upper())**3*a0*b0)/2
    lower = arb(minimum)-arb(coeff[2])*arb(tan_lower)*remainder
    assert lower > arb(".466")
    a,b,u,d,j = means
    ss,tt = a/(1-u), b/(1-u)
    assert all((s-half)**2 < val < (s+half)**2 for val,s in zip((ss,tt,u),(s12,s23,s13)))
    real = (1-u)*(d-ss*(1-tt)-(1-ss)*tt*u)/2
    phase = acb(arb(real),arb(j)).arg()
    assert arb("1.141")<phase and phase<arb("1.153")
    print("Minimum forced weight above",lower.lower().str(10),"feasible",weights[11431]/10**6,
          "witness phase",phase.str(12))


if __name__ == "__main__":
    verify()
