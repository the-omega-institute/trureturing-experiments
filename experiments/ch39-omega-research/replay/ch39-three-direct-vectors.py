"""求固定带荷几何与带色质量时伪装两处完整交换的最小正矢量数。"""
from pathlib import Path
from runpy import run_path

from flint import arb_poly

contact = run_path(str(Path(__file__).with_name("ch39-contact-current-direction.py")))
poles, mat, scheme, q = (contact[key] for key in ("poles", "mat", "scheme", "q"))
projection = mat["hypercharge"]


def effective_masses(first, second, e, x):
    # 正定同时对角化；第三个基本对称多项式由直接行列式计算。
    pencil = first.inv()*second
    trace = lambda m: sum(m[i,i] for i in range(3))
    s1 = trace(pencil)
    s2 = (s1**2-trace(pencil*pencil))/2
    roots = arb_poly([-pencil.det(), s2, -s1, 1]).complex_roots(tol=1e-13)
    masses = []
    assert len(roots) == 3
    for root in roots:
        assert root.imag.contains(0) and 1 < root.real < x
        masses.append((e*x*(root.real-1)/(x-root.real)).sqrt())
    return sorted(masses, key=lambda a:float(a.mid()))


def classify(logs, time, target):
    (length,z,rho,raw_d,raw_b,masses), _, _ = poles["retune"]["predictions"](logs,time)
    e, x = scheme["anchor"]/length, q("5/3")
    k, b, d = 1/rho, raw_b/rho, raw_d/rho
    kk = 1/(poles["offset"]+target/(15*length))
    assert kk < k
    h = mat["flow"]["running"]["hypercharge_window"](time+scheme["anchor"].log()/2).mid()**2
    hs, alphas, p = [], [], []
    for energy in (e, x*e):
        side, _, _, alpha, _ = contact["axes"](energy,z)
        hs.append(side)
        alphas.append(alpha)
        p.append(k*energy/(b+d*energy+4/side))
    # 两处分别使超荷、族系数相等，检查其余两方向严格为正。
    hh = h+(k-kk)*(alphas[0]-q("1/24"))
    dd = (kk*x*e/p[1]-4/hs[1]-b)/(x*e)
    assert 0 < hh < 1 and dd > 0 and alphas[1] > alphas[0]
    crit1 = (kk*e/p[0]-4/hs[0]-b)/e
    assert crit1 < dd
    first = mat["massive_exchange"](e,z,b,d)
    second = mat["massive_exchange"](x*e,z,b,d)
    assert mat["positive"](first) and mat["positive"](second-first) and mat["positive"](first*x-second)
    ratios = [m/masses[0] for m in effective_masses(first, second, e, x)]
    for ratio, bounds in zip(ratios, ((".636", ".638"), ("4.61", "4.63"), ("195.9", "196.3")), strict=True):
        assert q(bounds[0]) < ratio < q(bounds[1])
    print("three vector masses / triplet", *(r.str(11) for r in ratios), flush=True)
    for i,energy in enumerate((e,x*e)):
        side,g,v,alpha,beta = contact["axes"](energy,z)
        h_diff = (k-kk)*(alpha-q("1/24"))+h-hh
        p_diff = p[i]-kk*energy/(b+dd*energy+4/side)
        rank_two = g*g.transpose()*((k-kk)*beta)
        if i == 0:
            assert p_diff > 0 and h_diff.contains(0)
            rank_two += v*v.transpose()*p_diff
        else:
            assert h_diff > 0 and p_diff.contains(0)
            rank_two += projection*h_diff
        original = mat["massive_exchange"](energy,z,b,d)*k+projection*h
        other = mat["massive_exchange"](energy,z,b,dd)*kk+projection*hh
        contact["same_matrix"](original, other+rank_two)


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11430,time)
    for logs in poles["retune"]["log_sets"]:
        classify(logs, time, target)


if __name__ == "__main__":
    verify()
