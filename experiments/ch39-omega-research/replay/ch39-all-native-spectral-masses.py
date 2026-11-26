"""全部合格参考味共同真空都强迫第 105 圈的份额超过原生正谱上界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


RECORDS = (
    (815000, 873125, (180637030, 2626240, 87306859, -4621681), (
        (-2800000, 4300000, 438994173), (-1600000, -14900000, 2888306385),
        (2000000, -6500000, 799212886), (2000000, -5300000, 875025786),
        (2000000, -4100000, 1288565300),
    )),
    (873125, 931250, (173424493, 3208660, 83354978, -10662457), (
        (-2800000, 4300000, 8225997), (-1600000, -14900000, 2665705225),
        (2000000, -6500000, 883982579), (2000000, -5300000, 961092176),
        (2000000, -4100000, 1577412627),
    )),
    (931250, 960312, (146988381, -2025248, 75765634, -9188087), (
        (-2800000, 4300000, 691552474), (-1785714, -14285714, 913094727),
        (-1600000, -14900000, 690647891), (2000000, -6500000, 1121228310),
        (2000000, -4100000, 1278916827), (2000000, -2900000, 159799611),
    )),
    (960312, 989375, (148172287, 297972, 75561895, -7764093), (
        (-2800000, 5500000, 186158899), (-1600000, -14900000, 1565345613),
        (2000000, -6500000, 1019360621), (2000000, -4100000, 1985578454),
        (2000000, -2900000, 116880255),
    )),
    (989375, 1018437, (151450969, 3602577, 75399939, -15509901), (
        (-2800000, 5500000, 26290146), (-1600000, -14900000, 1801762722),
        (2000000, -6500000, 1251135876), (2000000, -5300000, 3524615),
        (2000000, -4100000, 2075234256),
    )),
    (1018437, 1032968, (137697010, -278449, 71466242, -7181477), (
        (-2800000, 5500000, 97580238), (-1600000, -14900000, 1171957941),
        (-1600000, 3100000, 305473331), (2000000, -6500000, 1015945544),
        (2000000, -4100000, 1841585407),
    )),
    (1032968, 1047500, (137801957, 11332720, 71169358, 14939841), (
        (-1785714, 3571429, 1352772990), (-1600000, -14900000, 1256546042),
        (-1600000, 3100000, 609257993), (2000000, -6500000, 418702707),
        (2000000, -5300000, 586122946),
    )),
    (1047500, 1062031, (120819337, 3568386, 67035946, 29338400), (
        (-1600000, -14900000, 397710944), (-1600000, 3100000, 1572838857),
        (2000000, -5300000, 753914812), (2000000, -4100000, 533395437),
    )),
    (1062031, 1076562, (113688347, 7814674, 66989145, 39225804), (
        (-2800000, 4300000, 64123474), (2000000, -5300000, 220442374),
        (2000000, -4100000, 906713327), (2000000, -2900000, 1494550534),
    )),
    (1076562, 1105625, (113857911, 9607881, 67520350, 41353174), (
        (2000000, -6500000, 13391880), (2000000, -5300000, 138937551),
        (2000000, -4100000, 1017841457), (2000000, -2900000, 1495165402),
    )),
    (1105625, 1163750, (113061476, 13292882, 68975855, 47351603), (
        (-1785714, 3571429, 74771116), (2000000, -4100000, 1074966733),
        (2000000, -2900000, 1468217528),
    )),
    (1163750, 1280000, (115759296, 21271655, 74024885, 55778038), (
        (-2800000, 3100000, 950887260), (-1785714, 3571429, 1368068713),
        (2000000, -2900000, 231656961),
    )),
)


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    load("common-all-phase-global-order")["verify"](root)
    exact = load("mincut-real-triangle-tomography")
    native = load("intrinsic-convolution-weights")
    angles = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    ctx.prec = 224
    labels = sorted({z for k, z in exact["enumerate_faces"]().values() if k != (0, 0)})
    assert len(labels) == 32 and {exact["conjugate"](z) for z in labels} == set(labels)
    native["native_bound"](exact, labels)
    tri = lib["terms"](shell)
    d = -angles["reference_angles"](flavor)[0]
    keys = ((1, 0), (0, 1), (1, -1))
    modes = []
    for q in range(1, 106):
        bare = lib["coeff"](tri, q, arb(0))
        soft = [sum((z**q for z, m in tri[k] if m), acb(0)) for k in keys]
        modes.append((bare, soft))
    assert RECORDS[0][0] == 815000 and RECORDS[-1][1] == 1280000
    assert all(v[1] == w[0] for v, w in zip(RECORDS, RECORDS[1:]))
    worst = arb(1)
    total = 0
    threshold = arb(".072")
    for left, right, directions, comparisons in RECORDS:
        assert left < right
        lo, hi = [arb(q)/10**8 for q in (left, right)]
        gx, gy, hx, hy = [arb(q)/10**10 for q in directions]
        cuts = [(arb(a)/10**8, arb(b)/10**8, arb(v)/10**10) for a, b, v in comparisons]
        assert all(lam > 0 for _, _, lam in cuts)
        mass = sum(v for _, _, v in cuts)
        for q, (bare, soft) in enumerate(modes, 1):
            constant = arb(int(q == 105)) + sum(lam*sum((v*acb(0,q*t).exp()).real
                       for v,t in zip(bare,(x,y,x-y))) for x,y,lam in cuts)
            united = [(-mass-acb(0,q*g))*v-acb(0,q*h)*r for v,r,g,h
                      in zip(bare,soft,(gx,gy,gx-gy),(hx,hy,hx-hy))]
            shift = acb(0, q*d).exp()
            constant += (united[2]*shift.conjugate()).real
            harmonic = united[0]+united[1]*shift
            # 仅余 Re[h exp(-iqt)]。用中点斜率与统一二阶导数界包住闭段。
            pending = [(lo, hi)]
            while pending:
                l, r = pending.pop()
                mid = (l+r)/2
                radius = (r-l)/2
                z = harmonic*acb(0,-q*mid).exp()
                lower = constant+z.real-radius*q*abs(z.imag)-(radius*q)**2*abs(harmonic)/2
                if lower > threshold:
                    total += 1
                    worst = min(worst, lower.lower())
                else:
                    assert r-l > arb("1e-16"), (q, l, r, lower)
                    pending.extend(((l,mid),(mid,r)))
    print("有理证书数、单阶覆盖末段数、全部目标低界",len(RECORDS),total,worst.str(17))


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else here
    verify(root)
