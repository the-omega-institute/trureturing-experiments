"""用谱前缀证明全双钟相位的共同单调质量障碍。"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
source = run_path(str(ROOT / "ch39-chiral-all-phase-thresholds.py"))
torus = source["torus"]
phase = torus["phase"]
ordered = source["ordered"]


def shell_bands(model):
    wide = torus["dimensionless"](model)
    coarse = {q: ordered([v for v,z,k in wide if z==q for _ in range(k)])
              for q in (0,*model["charges"])}
    # 既有谱带把单位边的五种荷限定在 (0,24)，三辐边与 c 的各自误差 <1e-7。
    pad = arb("1e-7")
    for q in model["charges"]:
        lower = [arb(0)]*19
        upper = [arb(24)]*19
        for i in range(1,96):
            cut = (arb(i) + arb(17)/37)/4
            core, terms = phase["schur_coefficients"](q,cut)
            npos, nneg = core, 16-core
            for term in terms:
                try:
                    sign = phase["scalar_range"](*term)
                except ValueError:
                    continue
                npos += sign
                nneg += 1-sign
            for j in range(nneg):
                upper[j] = upper[j].min(cut)
            for j in range(19-npos,19):
                lower[j] = lower[j].max(cut)
        values = []
        for lo,hi in zip(lower,upper):
            a,b = lo-pad,hi+pad
            assert a < b
            pts = [a,b] + [x for x in (arb("7.5"),arb(15)) if a < x < b]
            h = a*abs(a-15)
            for x in pts:
                h = h.union(x*abs(x-15))
            values.append(h + arb(0,b*pad))
        coarse[q] = [a.intersection(b) for a,b in zip(coarse[q],ordered(values))]
        assert all(x > 0 for x in coarse[q])
    return coarse


def verify():
    model, _, _, initial, _ = torus["local"]["inputs"]()
    ctx.prec = 224
    bands = shell_bands(model)
    initial[2] = 1 / arb(".117").union(arb(".119"))
    lmin,lmax = [(arb(n)/model["mz"]).log() for n in ("1e13","1e17")]
    # 用有理数累计 beta 投影与区间端点事件。
    d = (-11,17,-6)
    def proj(w):
        return sum(a*v for a,v in zip(d,w))
    def exact(weights):
        # 各 beta 系数是分母整除 30 的有理数。
        return [F(int((v*30).unique_fmpz()),30) for v in weights]
    fermions = {q:2*proj(exact(w)) for q,w in zip(model["charges"],model["weyl"])}
    scalar = proj(exact(model["higgs"]))
    vec = proj([0,-7,-F(21,2)])
    assert scalar>0 and sum(d)==0
    for full in (True,False):
        events = defaultdict(F)
        def add(weight, rows):
            for h in rows:
                left = int((10000*h.lower()).floor().unique_fmpz())
                right = int((10000*h.upper()).ceil().unique_fmpz())
                events[left if weight>0 else right] += weight

        for q,w in fermions.items():
            add(w,bands[q][3:])
        add(vec,bands[0])
        add(scalar,bands[3] if full else bands[3][3:])
        prefix = cap = F(0)
        for cutoff in sorted(events):
            prefix += events[cutoff]
            cap = max(cap,prefix)
        assert prefix == 16*sum(fermions.values()) + 18*vec + (19 if full else 16)*scalar
        def ball(v):
            return arb(v.numerator)/v.denominator
        sharp = (2*arb.pi()*proj(initial)
                 -(proj(model["beta"])+ball(cap))*lmin.union(lmax))/(arb.pi()*sum(map(abs,d)))
        # E_F−E_V 在 [0,5/42]，E_S−E_V = −16(E_F−E_V)/5。
        # 十六个荷三态同时带一个标量；另三个标量的额外误差只会降低投影。
        correction = 16*sum(max(F(0),w-(F(16,5)*scalar if q==3 else 0))
                            for q,w in fermions.items())*F(5,42)
        smooth = sharp - ball(correction)/(arb.pi()*sum(map(abs,d)))
        print("十九重标量" if full else "十六重标量", "前缀上界",cap,
              "尖阈值和有限减法下界",sharp.lower().str(12),smooth.lower().str(12))
        assert sharp > arb(".7048" if full else "1.9424")
        assert smooth > arb(".2426" if full else "1.4802")


if __name__=="__main__":
    verify()
