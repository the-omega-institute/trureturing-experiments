"""相同重阈值下，电荷重分配对规范交点的一次两圈插入之差。"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


color = run_path(str(Path(__file__).with_name("ch39-all-light-color-assignments.py")))
broken = run_path(str(Path(__file__).with_name("ch39-broken-electroweak-vacuum.py")))


def gauge_matrix():
    # 每行给出实/复系数、色维数、色指标、色 Casimir 和超荷平方。
    fund = (1, 3, F(1, 2), F(4, 3), F(1, 9))
    six = (1, 6, F(5, 2), F(10, 3), F(1, 9))
    adj = (F(1, 2), 8, 3, 3, 0)
    singlet = (1, 1, 0, 0, 1)
    def terms(fields):
        B = [[F(0), F(0)], [F(0), F(0)]]
        for n, d, t, c, y2 in fields:
            c1 = F(3, 5)*y2
            B[0][0] += n*4*d*c1**2
            B[0][1] += n*4*d*c1*c
            B[1][0] += n*4*t*c1
            B[1][1] += n*t*(4*c+2)
        return B
    old = terms([fund, six] + 3*[adj])
    new = terms([singlet] + 5*[adj])
    diff = [[a-b for a, b in zip(x,y)] for x,y in zip(new,old)]
    assert diff == [[F(32,25), -F(32,5)], [-F(4,5), 0]]
    return [[arb(x.numerator)/x.denominator for x in row] for row in diff]


def neutral_layers(bands):
    heavy_mass = 32*arb("2.435e18")/724**2*bands[2][15]
    mu = arb("246.22")**2/heavy_mass
    levels = [(name, mass2.sqrt(), g, -4 if mass2.is_zero() else 0)
              for name,g,_,mass2 in broken["light_spectrum"](arb(1))]
    levels.extend((f"差分{j}", arb(j).sqrt()*mu, 0, w)
                  for j,w in ((1,12), (2,-12), (3,4)))
    mz = arb("91.1876")
    assert arb(3).sqrt()*mu < arb("1e-6")
    candidates = []
    for name,mass,g,target in levels:
        if mass > mz:
            continue
        scalars, dirac = broken["minimal_partners"](g,target)
        if 2*mass < mz:
            # 电中性弱单态色表示中心平凡。实伴随需 8 个实标量或 8 个 Dirac。
            assert dirac < 8 and scalars <= 12
            if scalars >= 8:
                candidates.append(name)
        else:
            assert name in ("W","Z") and scalars < 6 and dirac < 3
            assert (scalars, dirac) == ((2,2) if name=="W" else (1,1))
    assert set(candidates) == {"u","d","s","c","b","差分1"}
    choices = list(combinations(candidates,5))
    assert len(choices) == 6
    print("电中性轻色谱的带标签选择数、四层首个非零质量/GeV", len(choices), mu.str(9))


def displacement(model, initial, bands, matrix, guess):
    two_pi = 2*arb.pi()
    beta = [b+d for b,d in zip(model["beta"], (arb(1)/5, 0, arb(5)/2))]
    rows = [(t, [arb(v)/30 for v in weights]) for t, weights in color["shell_table"](model, bands)]
    rows.append(((arb("172.5")/model["mz"]).log(), [arb(2)/15, arb(0), arb(5)/6]))
    window = (arb(guess, "3e8")/model["mz"]).log()
    assert all(not window.overlaps(t) and t>0 for t,_ in rows)
    open_rows = [(t,w) for t,w in rows if t < window]
    slopes = [b+sum(w[i] for t,w in open_rows) for i,b in enumerate(beta)]
    moments = [sum(t*w[i] for t,w in open_rows) for i in range(3)]
    root = (two_pi*(initial[0]-initial[1])+moments[0]-moments[1])/(slopes[0]-slopes[1])
    assert window.contains_interior(root)
    unified = initial[1] - (slopes[1]*root-moments[1])/two_pi
    inverse3 = unified + (slopes[2]*root-moments[2])/two_pi
    starting = [initial[0], initial[1], inverse3]
    # 固定中心折点的分段积分。真实折点平移至多 rad(t)，于是逆耦合误差
    # 在整段受 Σ |w|rad(t)/(2π) 控制。用同一界扩张每段的两个端点。
    centered = sorted([(t.mid(), w) for t,w in open_rows],key=lambda row: row[0])
    errors = [sum(t.rad()*abs(w[i]) for t,w in open_rows)/two_pi for i in range(3)]
    last = arb(0)
    inv = list(starting)
    slope = list(beta)
    integrals = [arb(0), arb(0), arb(0)]
    for end, weight in centered + [(root, [0,0,0])]:
        length = end-last
        assert length>=0
        for i in (0,2):
            nxt = inv[i] - slope[i]*length/two_pi
            # 上下平移后的两条仿射线包住真实逆耦合，倒数积分端点有序。
            assert inv[i]-errors[i]>0 and nxt-errors[i]>0
            def integral(offset):
                a, b = inv[i]+offset, nxt+offset
                if slope[i].is_zero():
                    return length/a
                return two_pi/slope[i]*(a/b).log()
            integrals[i] += integral(errors[i]).union(integral(-errors[i]))
            inv[i] = nxt
        slope = [s+w for s,w in zip(slope,weight)]
        last = end
    assert inv[0].overlaps(unified) and inv[2].overlaps(unified)
    i1,i3 = integrals[0],integrals[2]
    d1,d3 = [-(row[0]*i1+row[1]*i3)/(8*arb.pi()**2) for row in matrix]
    dt = two_pi*d1/(slopes[0]-slopes[1])
    dinv3 = -d3 + (slopes[2]-slopes[1])*dt/two_pi
    dalpha3 = -dinv3/inverse3**2
    print("I1、I3、δlog M、δα_s(M_Z)", *(v.str(10) for v in (i1,i3,dt,dalpha3)))
    return dt, dalpha3


def verify():
    chiral = color["gauge"]["chiral"]
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 224
    source = chiral["source"]["spectral_supplier"](model)
    matrix = gauge_matrix()
    for number, guess in [(22027,"8.02037076e15"), (28325,"8.02121470e15")]:
        values = source(acb.zeta_zero(number).imag)
        bands = {q:chiral["ordered"]([h for h,z in values if z==q]) for q in (0,*model["charges"])}
        print("正零点", number)
        neutral_layers(bands)
        dt, da = displacement(model,initial,bands,matrix,guess)
        assert arb(".01227") < dt < arb(".01229") and arb(".00156") < da < arb(".00158")


if __name__ == "__main__":
    verify()
