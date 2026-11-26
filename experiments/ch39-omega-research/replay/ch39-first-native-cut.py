"""在真实零点前缀内以荷二曲率带及最小割分界逐枚核验三窗。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root):
    window = run_path(str(root / "ch39-curvature-window-density-cap.py"))
    inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))
    actual = run_path(str(root / "ch39-first-three-mode-observables.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    geometry = shell["geo"]
    model, initial = actual["model_inputs"]()
    physical = actual["joint"]
    ctx.prec = 128
    low, high = arb("18.94085"), arb("18.944146")
    left, right = window["eigen_band"]()
    assert low < left < right < high
    # 五荷各恰有三个轻模时，第十六质量超过最小割，故也超过 9。
    assert shell["soft"]["point"]["study"]["offset"] > 0
    assert 9+(1-arb(23)/25*arb(2)**-20/14**2)**2 > 9
    separators = [inertia["schur_coefficients"](2, cut) for cut in (low, high)]

    def below(terms, u, v):
        fixed, scalars = terms
        phases = [c+2*(a*u+b*v+d*u.conjugate()*v).real for c,a,b,d in scalars]
        assert all(s > 0 or s < 0 for s in phases)
        return 19-fixed-sum(s>0 for s in phases)

    ctx.prec = 76
    zeros = acb.zeta_zeros(1, 8359)
    assert len(zeros) == 8359 and zeros[0].imag > 14
    ctx.prec = 128
    survivors = []
    for n, zero in enumerate(zeros, 1):
        phases = geometry["phases"](zero.imag)
        u, v = (phases[d]**2 for d in (724,804))
        if below(separators[0], u, v) < 16 <= below(separators[1], u, v):
            survivors.append(n)
    assert survivors == [2465, 5882, 6514, 8035, 8236, 8359]
    print("曲率必要带的候选编号",survivors, flush=True)
    coeff = (9*arb(3).log()-12*arb(2).log())/2
    scale = arb("246.22")**2/(32*arb("2.435e18")/724**2*arb("1.22089e19"))
    offset = shell["soft"]["point"]["study"]["offset"]
    accepted = []
    profiles = []
    for n in survivors:
        ctx.prec = 192
        g = acb.zeta_zero(n).imag
        delta = arb(23)/25*arb(2)**-20/(g*g+offset)
        threshold = 9+(1-delta)**2
        raw = model["shell"]["shell_spectra"]((0,*model["charges"]), zero_index=n)
        center = sum(raw[0])/18
        assert center.overlaps(15-delta/3)
        bands = {q: physical["chiral"]["ordered"]([v*abs(v-center) for v in raw[q]])
                 for q in raw}
        counts = []
        for q in (1, 2, 3, 4, 6):
            assert all(x < threshold or x > threshold for x in bands[q])
            counts.append(sum(x < threshold for x in bands[q]))
        profiles.append(counts)
        print(n, "各荷的轻模数", counts, flush=True)
        if counts != [3]*5:
            continue
        a = (scale/bands[2][15])**2
        vac = physical["vacuum"]
        fields = vac["paired_spectrum"](vac["minimal"]["fields"](bands), min_gap="1e-18")
        rows = fields + vac["low_sector"](4, a)
        assert vac["newton"](rows, arb(".03")) < 1 < vac["newton"](rows, arb(".033"))
        ratio = a**2*coeff/(arb.pi()*arb("1e-122"))*(1-arb(0).union(a/(arb(".03")*coeff)))
        assert arb(".998") < ratio < arb("1.002") or ratio < arb(".998") or ratio > arb("1.002")
        print(n, "四层曲率比", ratio.str(15), flush=True)
        if ratio < arb(".998") or ratio > arb("1.002"):
            continue
        with redirect_stdout(StringIO()):
            _,_,alpha = physical["gauge"]["crossing_enclosure"](
                model, initial, bands, False, 724, "29.04", ".1184", recenter=True)
        print(n, "两圈强耦合下上界", alpha.lower().str(12), alpha.upper().str(12), flush=True)
        assert arb(".117") < alpha < arb(".119") or alpha < arb(".117") or alpha > arb(".119")
        if arb(".117") < alpha < arb(".119"):
            accepted.append(n)
            if n == 8236:
                assert bands[6][1] < 9 < bands[6][2] < threshold < bands[6][3]
                assert all(bands[q][2] < 9 < bands[q][3] for q in (1, 2, 3, 4))
                assert arb("9.20607325462") < max(bands[q][2] for q in model["charges"]) < arb("9.20607325464")
                assert arb("10.58546452370") < min(bands[q][3] for q in model["charges"]) < arb("10.58546452372")
                assert arb(".11830581072") < alpha < arb(".11831197167")
                assert arb("1.00093825580") < ratio < arb("1.00093825581")
            print(n, "高度、最小割、最大第三质量、最小第四质量",
                  g.str(18), threshold.str(18), max(bands[q][2] for q in model["charges"]).str(16),
                  min(bands[q][3] for q in model["charges"]).str(16), flush=True)
    assert profiles == [[3,3,4,3,2], *[[3,3,4,3,3]]*3, *[[3]*5]*2]
    assert accepted == [8236, 8359]
    print("前8359枚实际零点的三窗命中", accepted, flush=True)


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
