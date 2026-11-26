"""区间核验 #11661 可见谱的两个同侧三截止排除区间与曲率方向反例。

依赖 python-flint。沿用 ch39-cutoff-scale-stationarity.py 的固有时模型和谱输入。
"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

SOURCE = run_path(str(Path(__file__).with_name("ch39-cutoff-scale-stationarity.py")))
chiral = SOURCE["minimal"]["chiral"]
model, *_ = chiral["source"]["mean"]["inputs"]()
ctx.prec = 192
read_spectrum = chiral["source"]["spectral_supplier"](model)
points = read_spectrum(acb.zeta_zero(11661).imag)
layers = {
    charge: chiral["ordered"]([value for value, label in points if label == charge])
    for charge in (0, *model["charges"])
}
rows = SOURCE["minimal"]["fields"](layers)
ctx.prec = 192
log_term_bound = SOURCE["two_cutoff_obstruction"](rows)
assert len(rows) == 99
assert log_term_bound < arb("0.0001221")

target = arb("1e-122")

def expint(n, u):
    if u.is_zero():
        return arb(1)/(n-1)
    return u.lower().expint(n).union(u.upper().expint(n))

def response_derivatives(logx):
    """各行给出函数及其对 log x 的一、二阶导数。"""
    x = logx.exp()
    argument = [z/x for _, _, z in rows]
    g = [g for g, _, _ in rows]
    a = [a for _, a, _ in rows]
    uexp = [[(u**m)*(-u).exp() for u in argument] for m in range(3)]
    phi0, phi1, phi2 = [x * sum(c*v for c,v in zip(a,mom)) for mom in uexp]
    psi0, psi1, psi2 = [x*x * sum(c*v for c,v in zip(g,mom)) for mom in uexp]
    psi = [psi0, 2*psi0+psi1, 4*psi0+3*psi1+psi2]
    phi = [phi0, phi0+phi1, phi0+phi1+phi2]
    K = x*sum(c*expint(2,u) for c,u in zip(a,argument))/(2*arb.pi())
    newton = [K-1, phi0/(2*arb.pi()), (phi0+phi1)/(2*arb.pi())]
    G = [x*z*g/(8*arb.pi()) for g,_,z in rows]
    A = [target*a*z/(2*arb.pi()) for _,a,z in rows]
    residual = [
        sum(c*expint(2,u)+b*expint(1,u) for c,b,u in zip(G,A,argument) if not u.is_zero()),
        sum((c+b)*e for c,b,e in zip(G,A,uexp[0])),
        sum((c*(1+u)+b*u)*e for c,b,u,e in zip(G,A,argument,uexp[0]))
    ]
    return psi, newton, phi, residual

def curvature(logx):
    psi, newton, phi, residual = response_derivatives(logx)
    return [arb_mat([psi,newton,column]).det() for column in (phi,residual)]

lim = arb("0.0001221")
pi2 = 2*arb.pi()

def ball(mid, radius):
    return (mid-radius).union(mid+radius)

def newton_wronskian(s, p, log_term, dp, ddp):
    s0,s1,s2 = s
    p0,p1,p2 = p
    return arb_mat([
        [s0, 2*s0+s1, 4*s0+3*s1+s2],
        [-1-log_term, -dp/pi2, -ddp/pi2],
        [p0, p0+p1, p0+p1+p2],
    ]).det()

def interval(x):
    u = [z/x for _,_,z in rows]
    g = [g for g,_,_ in rows]
    a = [a for _,a,_ in rows]
    moments = [[(y**m)*(-y).exp() for y in u] for m in range(3)]
    s = [sum(c*v for c,v in zip(g,ms)) for ms in moments]
    p = [sum(c*v for c,v in zip(a,ms)) for ms in moments]
    L = sum(a*z*expint(1,y) for (_,a,z),y in zip(rows,u) if not z.is_zero())/pi2
    return newton_wronskian(s,p,L,x*p[1],x*p[2])

def boundary():
    lo = arb(-29)
    massless = [(g,a) for g,a,z in rows if z.is_zero()]
    massive = [(g,a,z) for g,a,z in rows if not z.is_zero()]
    zmin = min(z.lower() for _,_,z in massive)
    M = arb(32)
    assert lo.exp()*M < zmin
    def left_moments(coefficients, main):
        return [ball(main if i==0 else 0, coefficients*M**i*(-M).exp()) for i in range(3)]
    s = left_moments(sum(abs(g) for g,_,_ in massive),sum(g for g,_ in massless))
    p = left_moments(sum(abs(a) for _,a,_ in massive),sum(a for _,a in massless))
    left = newton_wronskian(s,p,ball(0,lim),ball(0,lo.exp()*abs(p[1]).upper()),
               ball(0,lo.exp()*abs(p[2]).upper()))
    assert left>0, left
    print("low tail", left.str(12), flush=True)


def confirm(left, right, expression, sign, max_boxes=12500):
    pending = [(arb(left),arb(right))]
    count = 0
    margin = None
    while pending:
        count+=1
        a,b = pending.pop()
        v = sign * expression(a.union(b))
        if v>0:
            v = v.lower()
            margin = v if margin is None else min(margin,v)
            continue
        if count>=max_boxes or b-a<arb("1e-7"):
            raise RuntimeError((left,right,count,a,b,v))
        m = ((a+b)/2).mid()
        pending.extend([(a,m),(m,b)])
    print(left,right,count,"minimum certified margin",margin.str(12),flush=True)

def verify():
    nodes = ("-26", "-21.5", "-14.5")
    vectors = [curvature(arb(x)) for x in nodes]
    for x, (u, v) in zip(nodes, vectors):
        print("曲率行列式", x, u.str(15), v.str(15))
    for (u, v), (a, b) in zip(vectors, vectors[1:] + vectors[:1]):
        orientation = u * b - v * a
        assert orientation > 0, orientation
        print("相邻方向行列式", orientation.str(15))

    boundary()
    confirm("-29", "-26.8", lambda logx: interval(logx.exp()), 1)
    confirm("-14.55", "0", lambda logx: curvature(logx)[1], 1)


if __name__ == "__main__":
    verify()
