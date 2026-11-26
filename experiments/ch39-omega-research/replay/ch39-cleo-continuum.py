"""两种轻色谱的两圈 QCD 跑动与 CLEO 七点连续强子截面。"""
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


light = run_path(str(Path(__file__).with_name("ch39-broken-electroweak-vacuum.py")))["light_spectrum"]
# CLEO 0706.2813，表 VII。E、R、统计误差、共用系统误差、独立系统误差。
cleo = """
10.538 3.591 .003 .067 .049
10.330 3.491 .006 .058 .055
9.996  3.497 .004 .064 .043
9.432  3.510 .005 .066 .037
8.380  3.576 .024 .058 .025
7.380  3.550 .019 .058 .020
6.964  3.597 .033 .057 .020
"""
quark_charges = {"u": 2, "d": -1, "s": -1, "c": 2}


def color_coefficients():
    # 复标量计一次，实标量计半次。列为有效 Dynkin 指数、色 Casimir。
    triplet = (Fraction(1, 2), Fraction(4, 3))
    sextet = (Fraction(5, 2), Fraction(10, 3))
    real_adjoint = (Fraction(3, 2), 3)
    for reps in ([triplet, sextet] + 3*[real_adjoint], 5*[real_adjoint]):
        db = sum(t/3 for t, c in reps)
        dB = sum(t*(4*c + 2) for t, c in reps)
        assert (db, dB) == (Fraction(5, 2), 105)
    return arb(23)/3 - arb(5)/2, arb(116)/3 - 105


def strong_at(alpha_z, energy, coeffs):
    b0, b1 = coeffs
    c = b1/(4*arb.pi()*b0)
    def primitive(u):
        return u - c*(u+c).log()
    target = primitive(1/alpha_z) + b0/(2*arb.pi())*(energy/arb("91.1876")).log()
    lo, hi = arb(4), arb(10)
    assert lo+c > 0 and primitive(lo) < target < primitive(hi)
    def root(endpoint):
        left, right = lo, hi
        for _ in range(60):
            mid = (left + right)/2
            if primitive(mid) < endpoint:
                left = mid
            elif primitive(mid) > endpoint:
                right = mid
            else:
                break
        return left.union(right)
    inverse = root(target.lower()).union(root(target.upper()))
    return 1/inverse


def corrections(beta):
    # Hikasa 与 Hisano，式 (3.16)–(3.18)，总截面含虚修正与实胶子发射。
    b2 = beta**2
    z = (1-beta)/(1+beta)
    log_ratio = -z.log()
    common = (1+b2)/beta*(4*z.polylog(2) + 2*(-z).polylog(2)
              - (3*(2/(1+beta)).log() + 2*beta.log())*log_ratio)
    common -= 3*(4/(1-b2)).log() + 4*beta.log()
    scalar = common + (-3+10*b2+5*b2**2)*log_ratio/(4*beta**3) + 3*(1+b2)/(2*b2)
    vector = common + (33+22*b2-7*b2**2)*log_ratio/(8*beta*(3-b2)) + 3*(5-3*b2)/(4*(3-b2))
    return scalar, vector


def rates(energy, alpha, masses):
    quarks = arb(0)
    fundamental, sextet = {}, {}
    for name, mass in masses.items():
        if 2*mass >= energy:
            fundamental[name] = sextet[name] = arb(0)
            continue
        beta = (1 - (2*mass/energy)**2).sqrt()
        sc, vec = corrections(beta)
        # R 对点状无质量 μ 对截面归一。色标量电荷绝对值均为 1/3。
        fundamental[name] = beta**3/12*(1 + arb(4)/3*alpha/arb.pi()*sc)
        sextet[name] = beta**3/6*(1 + arb(10)/3*alpha/arb.pi()*sc)
        if name in quark_charges:
            quarks += (arb(quark_charges[name]**2)/3 * beta*(3-beta**2)/2
                       * (1 + arb(4)/3*alpha/arb.pi()*vec))
    return quarks, fundamental, sextet


def chisq(measured, predicted):
    # 对一个单位方差的共同系统偏移解析极小化。
    w = [1/(stat**2 + unc**2) for _, _, stat, _, unc in measured]
    r = [p-row[1] for p, row in zip(predicted, measured)]
    h = [row[3] for row in measured]
    pull = sum(wi*ri*hi for wi, ri, hi in zip(w, r, h))/(1+sum(wi*hi**2 for wi, hi in zip(w, h)))
    return pull**2 + sum(wi*(ri-hi*pull)**2 for wi, ri, hi in zip(w, r, h))


def compare(number, mz_alpha, data, masses, coeffs, scale):
    rates_at = [rates(row[0], strong_at(mz_alpha, scale*row[0], coeffs), masses) for row in data]
    neutral = [q for q, _, _ in rates_at]
    charged = {}
    crossings = []
    null = chisq(data, neutral)
    for a, b in permutations(masses, 2):
        signal = [fund[a] + sext[b] for _, fund, sext in rates_at]
        full = chisq(data, [q+s for q, s in zip(neutral, signal)])
        half = chisq(data, [q+s/2 for q, s in zip(neutral, signal)])
        # χ²(f)−χ²(0)=A f²+B f，f 为标量相对普通强子事例的共同接收效率。
        quadratic = 2*(full-2*half+null)
        linear = full-null-quadratic
        assert quadratic > 0
        crossings.append((-linear/quadratic).upper())
        charged[a, b] = full
    best = min(charged, key=lambda k: charged[k].lower())
    assert all(charged[best] < value for key, value in charged.items() if key != best)
    high = max(v.upper() for v in charged.values())
    print(number, len(data), "μ/E", scale.str(4), "零直接项 χ²、最小最大带电 χ²、最小的(基本,六维)",
          null.str(7), charged[best].str(7), high.str(7), best)
    crossing = max(crossings)
    print("超过此相对效率，全部带电分配的 χ² 高于电中性分配", crossing.str(7))
    assert null < charged[best] and crossing < 1
    return null, charged[best], crossing


def verify():
    ctx.prec = 224
    coeffs = color_coefficients()
    scalar, vector = corrections(1-arb("1e-16"))
    assert abs(scalar-3) < arb("1e-12") and abs(vector-arb(3)/4) < arb("1e-12")
    b = arb("1e-6")
    scalar, vector = corrections(b)
    assert abs(scalar-arb.pi()**2/(2*b)+2) < arb(".001")
    assert abs(vector-arb.pi()**2/(2*b)+4) < arb(".001")
    data = [tuple(map(arb, row.split())) for row in cleo.splitlines() if row.strip()]
    masses = {name: z.sqrt() for name, _, _, z in light(arb(1)) if name in (*quark_charges, "b")}
    for number, alpha in [(22027, ".1188549474"), (28325, ".1187872682")]:
        alpha_z = arb(alpha, "3e-10")
        print(number, "两圈纯 QCD 在 10 GeV 的强耦合", strong_at(alpha_z, arb(10), coeffs).str(9))
        null, best, crossing = compare(number, alpha_z, data, masses, coeffs, arb(1))
        assert 6 < null < arb("6.1") and arb("13.5") < best < arb("13.6")
        assert arb(".338") < crossing < arb(".340")
        off_threshold = [row for row in data if not row[0].overlaps(arb("8.380"))]
        assert len(off_threshold) == 6
        compare(number, alpha_z, off_threshold, masses, coeffs, arb(1))
        for scale in (arb(3)/4, arb(3)/2):
            assert all(scale*row[0] > masses["b"] for row in data)
            compare(number, alpha_z, data, masses, coeffs, scale)


if __name__ == "__main__":
    verify()
