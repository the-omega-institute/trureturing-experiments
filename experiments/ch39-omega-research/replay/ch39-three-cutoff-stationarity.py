"""区间核验 #11661 可见谱的同侧与混合三截止排除区间、曲率方向反例。

依赖 python-flint。沿用 ch39-cutoff-scale-stationarity.py 的固有时模型和谱输入。
"""
from fractions import Fraction
from functools import cache
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


def small_node_pair():
    """排除 0 < x₁ < x₂ ≤ exp(−30.6) 及 exp(−26.8) < x₃ ≤ 1。"""
    cutoff_log = arb("-30.6")
    small_cutoff = cutoff_log.exp()
    x = arb(0).union(small_cutoff)
    massless = [(g, a) for g, a, z in rows if z.is_zero()]
    massive = [(g, a, z) for g, a, z in rows if not z.is_zero()]
    min_argument = arb(32)
    assert all(z > min_argument * small_cutoff for _, _, z in massive)
    tail = (-min_argument).exp()
    g_tail = sum(abs(g) for g, _, _ in massive) * tail
    a_tail = sum(abs(a) for _, a, _ in massive) * tail
    g_zero = sum(g for g, _ in massless)
    a_zero = sum(a for _, a in massless)

    # E₂(u) ≤ exp(−u)，故 ΣαE₂(u) 与 Σαexp(−u) 共用包围区间。
    a_sum = ball(a_zero, a_tail)
    g_exp = ball(g_zero, g_tail)
    denominator = x * a_sum / pi2 - 1
    d_k = a_sum / pi2
    phi = x * a_sum
    # u exp(−u) 在 u ≥ 32 上由 32 exp(−32) 控制。
    d_phi = a_sum + ball(0, min_argument * a_tail)
    d_psi = x * (2 * g_exp + ball(0, min_argument * g_tail))
    d_a = (d_phi * denominator - phi * d_k) / denominator**2
    d_t = (d_psi * denominator - x**2 * g_exp * d_k) / denominator**2
    assert a_sum > 0 and denominator < 0
    assert d_a < -arb(1)/7
    assert d_t < 130 * small_cutoff

    # 由 A'(x) + slope*T'(x) < 0，左弦斜率小于 −slope。
    slope = arb("2e10")
    assert 1 / (910 * small_cutoff) > slope
    assert d_a + slope * d_t < 0

    # T=Ψ/(K−1) 为正且严格递增，A=Φ/(K−1) 在小节点处为负。
    # 右弦分子 A₃−A₂ > A₃ > response_floor*T₃，右弦分母 > (1−ratio_limit)*T₃。
    psi_s, k_s, _, _ = response_derivatives(cutoff_log)
    psi_h, k_h, _, _ = response_derivatives(arb("-26.8"))
    ratio_limit = arb("0.01")
    ratio = (psi_s[0] / k_s[0]) / (psi_h[0] / k_h[0])
    assert ratio < ratio_limit
    response_floor = arb("-1.9e10")
    assert response_floor / (1 - ratio_limit) > -slope

    def ratio_margin(logx):
        psi, _, phi, _ = response_derivatives(logx)
        return phi[0] / psi[0] - response_floor

    confirm("-26.8", "0", ratio_margin, 1, max_boxes=5000)
    print("small-node derivative margins", (-d_a).str(12), (d_t / small_cutoff).str(12))
    print("small-node T ratio", ratio.str(12))


def separated_small_node():
    """排除 x₁ ≤ exp(−31)、x₂ ≥ exp(−30.6)、x₃ ≥ exp(−26.8) 的有序三点解。"""
    cutoff = arb(-31).exp()
    small = arb(0).union(cutoff)
    massless = [(g, a) for g, a, z in rows if z.is_zero()]
    massive = [(g, a, z) for g, a, z in rows if not z.is_zero()]
    assert all(z > 32*cutoff for _, _, z in massive)
    tail = arb(-32).exp()
    alpha = ball(sum(a for _, a in massless), tail * sum(abs(a) for _, a, _ in massive))
    polarization = ball(sum(g for g, _ in massless), tail * sum(abs(g) for g, _, _ in massive))
    denominator = small * alpha / pi2 - 1
    assert denominator < 0
    # 对 u≥32，两项 E₁(u) 与 E₂(u) 均由 exp(−32) 控制。
    residual_bound = tail * (
        cutoff * sum(abs(g*z) for g, _, z in massive) / (4*pi2)
        + target * sum(abs(a*z) for _, a, z in massive) / pi2
    )
    initial_t = small**2 * polarization / denominator
    initial_responses = (small * alpha / denominator, ball(0, residual_bound) / denominator)

    @cache
    def chart(left, right):
        """T，两个 U=响应/Ψ，及其对 log x 的导数与 T'/T。"""
        lo = arb(left.numerator) / left.denominator
        hi = arb(right.numerator) / right.denominator
        psi, newton, phi, residual = response_derivatives(lo.union(hi))
        t = psi[0] / newton[0]
        if not (t > 0 and t-initial_t > 0):
            return None
        values = tuple(response[0] / psi[0] for response in (phi, residual))
        derivatives = tuple((response[1]-value*psi[1]) / psi[0]
                            for response, value in zip((phi, residual), values))
        log_t_derivative = psi[1]/psi[0] - newton[1]/newton[0]
        if left < right:
            midpoint = (left+right)/2
            center = chart(midpoint, midpoint)
            if center is not None:
                radius = (right-left)/2
                width = arb(radius.numerator) / radius.denominator
                values = tuple(value.intersection(mid + ball(0, width*abs(derivative)))
                               for value, mid, derivative in zip(values, center[1], derivatives))
        return t, values, derivatives, log_t_derivative

    def slope_derivatives(t, values, derivatives, log_t_derivative):
        # S=(U−A₁/T)/(1−T₁/T)，这里计算 S' 乘以严格正的 (1−T₁/T)²。
        ratio = initial_t/t
        return tuple(derivative*(1-ratio) + log_t_derivative*(initial/t-ratio*value)
                     for value, derivative, initial in zip(values, derivatives, initial_responses))

    lower = Fraction(-134, 5)
    middle = Fraction(-153, 5)
    upper = Fraction(0)
    directions = []
    for point in (lower, Fraction(-26), Fraction(-21)):
        entry = chart(point, point)
        assert entry is not None
        directions.append(slope_derivatives(*entry))
    for (a, b), (c, d) in zip(directions, directions[1:] + directions[:1]):
        assert a*d-b*c > 0

    def separated(left, left_end, right_start, right):
        if right-left < Fraction(1, 2):
            span = chart(left, right)
            if span is not None and any(v > 0 or v < 0 for v in slope_derivatives(*span)):
                return True
        first = chart(left, left_end)
        second = chart(right_start, right)
        if first is None or second is None:
            return False
        t0, values0, _, _ = first
        t1, values1, _, _ = second
        inverse_difference = 1/t0 - 1/t1
        # 约去共同的 A₁T₁ 项。让每个 U 乘另一节点的分母，避免区间重复计入 U₂T₁/T₂。
        for u0, u1, initial in zip(values0, values1, initial_responses):
            difference = u1*(1-initial_t/t0) - u0*(1-initial_t/t1) + initial*inverse_difference
            if difference > 0 or difference < 0:
                return True
        return False

    pending = [(lower, upper, lower, upper), (middle, lower, lower, upper)]
    checked = 0
    while pending:
        a, b, c, d = pending.pop()
        if a >= d:
            continue
        checked += 1
        if separated(a, b, c, d):
            continue
        if checked >= 180000 or max(b-a, d-c) < Fraction(1, 100000):
            raise RuntimeError(("unresolved ordered box", checked, a, b, c, d))
        if b-a >= d-c:
            midpoint = (a+b)/2
            pending.extend([(a, midpoint, c, d), (midpoint, b, c, d)])
        else:
            midpoint = (c+d)/2
            pending.extend([(a, b, c, midpoint), (a, b, midpoint, d)])
    print("separated-small-node boxes", checked, flush=True)


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
    def mixed_curvature(logx):
        newton, vacuum = curvature(logx)
        return arb("2e-8") * newton + vacuum

    confirm("-25.3", "0", mixed_curvature, 1)
    small_node_pair()
    separated_small_node()


if __name__ == "__main__":
    verify()
