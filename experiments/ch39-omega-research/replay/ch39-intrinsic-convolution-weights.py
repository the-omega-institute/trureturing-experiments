"""有向裸面谱的105阶份额、一般酉一阶差分与三面闭钟相位字。"""
from collections import Counter
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


DIRECTION_WITNESSES = (
    ((0, 1003479), (1, 1003479), (6, 1110976), (7, 1110976), (10, 368264), (11, 368264), (20, 413588), (21, 413588)),
    ((0, 868069), (1, 1070503), (6, 1060518), (7, 1085723), (10, 52748), (11, 658684), (20, 629167), (21, 169824)),
    ((0, 675687), (1, 1066716), (7, 997805), (11, 864859), (20, 786262), (28, 1054579)),
    ((0, 368996), (1, 953711), (2, 230165), (7, 864697), (8, 450510), (11, 1006289), (18, 189240), (20, 967486), (28, 630043)),
    ((0, 8484), (1, 773374), (2, 561186), (5, 174068), (7, 688314), (8, 1104432), (11, 1092892), (15, 8097), (17, 7652), (18, 457640), (20, 1094427)),
    ((1, 620747), (2, 747538), (5, 406103), (7, 347672), (11, 1025337), (12, 1057894), (15, 252630), (17, 315226), (18, 691393), (20, 1138336)),
    ((1, 429429), (2, 879548), (5, 609429), (11, 827999), (12, 6282), (15, 486963), (17, 614874), (18, 849251), (20, 1036628), (26, 984998), (31, 131991)),
    ((1, 208511), (2, 955266), (5, 776367), (15, 663423), (17, 759515), (18, 998623), (20, 957957), (22, 816503), (26, 151928), (31, 914930)),
    ((2, 963600), (5, 894006), (15, 789667), (17, 836479), (18, 1087768), (20, 825917), (22, 518658), (24, 463688), (31, 840309)),
    ((2, 872089), (4, 79947), (5, 966119), (15, 862822), (17, 858062), (18, 1105707), (20, 632026), (24, 1014574), (31, 719670)),
    ((2, 704926), (4, 329744), (5, 968595), (15, 909194), (16, 252878), (17, 891164), (18, 1013576), (20, 352017), (24, 1007831), (25, 284500), (31, 438007)),
    ((2, 492975), (3, 16394), (4, 556800), (5, 904607), (14, 178041), (15, 912195), (16, 554031), (17, 897349), (18, 708412), (20, 38353), (24, 948567), (25, 645737), (31, 102661)),
    ((2, 263207), (3, 263207), (4, 755940), (5, 755940), (14, 749526), (15, 749526), (16, 797989), (17, 797989), (18, 122820), (19, 122820), (24, 867819), (25, 867819)),
)


ADJACENT_WITNESSES = (
    ((2, 272529), (3, 272529), (4, 700315), (5, 700315), (18, 896327), (19, 896328), (20, 1020261), (21, 1020261), (26, 981754), (27, 981754)),
    ((2, 25284), (3, 501399), (4, 862952), (5, 490214), (13, 991961), (18, 1013290), (19, 736244), (20, 925861), (21, 1098592), (22, 475104), (26, 459816)),
    ((3, 712450), (4, 944453), (5, 252339), (9, 1005882), (18, 1078054), (19, 497654), (20, 783845), (21, 1038484), (22, 728460), (24, 189290)),
    ((3, 875690), (4, 958470), (9, 553755), (18, 1053284), (19, 234872), (20, 564964), (21, 930762), (24, 922289), (29, 459595)),
    ((3, 952608), (4, 899241), (16, 254733), (18, 950068), (20, 280517), (21, 767235), (24, 950011), (29, 1002365)),
    ((0, 186782), (3, 943329), (4, 797209), (7, 1037948), (11, 79858), (14, 275941), (16, 555033), (18, 571849), (21, 637847), (24, 930191)),
    ((0, 407206), (3, 863990), (4, 644612), (7, 1088236), (11, 337766), (14, 829189), (16, 778967), (21, 437595), (24, 830336)),
    ((0, 602191), (3, 725053), (4, 446888), (7, 1051244), (11, 588509), (14, 898541), (16, 867299), (21, 196540), (24, 584473), (30, 111801)),
    ((0, 756925), (1, 9509), (3, 534770), (4, 220093), (7, 1005907), (11, 871319), (14, 882898), (16, 851900), (24, 192691), (30, 485986)),
    ((0, 933682), (1, 366075), (3, 203825), (7, 904387), (11, 1055288), (14, 836705), (16, 839291), (30, 754714)),
    ((0, 1046793), (1, 665853), (7, 710130), (11, 1113528), (14, 755732), (16, 819435), (30, 902275)),
    ((0, 1054026), (1, 857724), (7, 388947), (10, 78550), (11, 1053970), (14, 608039), (15, 150069), (16, 720991), (17, 216242), (30, 876022)),
    ((0, 989328), (1, 989328), (6, 28017), (7, 28017), (10, 1018646), (11, 1018646), (14, 395266), (15, 395266), (16, 523540), (17, 523540)),
)


def angular_certificate(wave, witnesses, corrections):
    maximum, slack = arb(0), arb(1)
    for j, row in enumerate(witnesses):
        direction = acb(0, -arb.pi()*j/12).exp()
        dual = [sum(arb(a)*wave[i][t]/10**8 for i, a in row) for t in range(105)]
        for t, correction in corrections.items():
            dual[t] += arb(correction)
        for vector in wave:
            gap = sum((a.conjugate()*b).real for a, b in zip(dual, vector))-(direction*vector[-1]).real
            slack = min(slack, gap.lower())
        norm_bound = sum(abs(b)**2 for b in dual)/(arb.pi()/24).cos()**2
        maximum = max(maximum, norm_bound.upper())
    return maximum, slack


def native_bound(exact, labels):
    # 算在实 Hilbert 空间 C^105 上，点积为复内积的实部。
    units = [exact["unit"](z) for z in labels]
    wave = [[u**q for q in range(1, 106)] for u in units]
    assert all(v[1].real > arb(".95") for v in wave)
    maximum, slack = angular_certificate(wave, DIRECTION_WITNESSES, {1: ".00001"})
    assert maximum < arb(".067563") and slack > arb("8.6e-6")
    print("原生32点的105阶份额上界、逐点方向余量", maximum.str(17), slack.str(15))
    # 正因子 |z-1| 吸收进任意非负面质量，只留下列方向的变化。
    adjacent = [[v*(u-1)/abs(u-1) for v in row] for u, row in zip(units, wave)]
    maximum, slack = angular_certificate(adjacent, ADJACENT_WITNESSES, {0: ".0001", 1: "-.0001"})
    assert maximum < arb(".064915") and slack > arb("3.8e-7")
    print("相邻圈差分的105阶份额上界、逐点方向余量", maximum.str(17), slack.str(15))


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    exact = load("mincut-real-triangle-tomography")
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    fixed = load("height-independent-flavor-vacuum")
    reference = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    angles = load("monotone-winding-action")
    ctx.prec = 224
    table = Counter((k, z) for k, z in exact["enumerate_faces"]().values() if k != (0, 0))
    labels = sorted({z for _, z in table})
    assert len(table) == len(labels) == 32 and set(table.values()) == {3}
    assert set(labels) == {exact["conjugate"](z) for z in labels}
    keys = ((1, 0), (0, 1), (1, -1))
    groups = [[exact["unit"](z) for key, z in table if key == k] for k in keys]
    assert list(map(len, groups)) == [14, 14, 4]
    native_bound(exact, labels)
    tri = lib["terms"](shell)
    x = arb("-.01046")
    y = x-reference["reference_angles"](flavor)[0]
    fixed["required_high_weights"](lib, tri, x, y)

    # 对 |theta|≥2*pi/105，由实部几何和控制分母，近零区用 sinc 递减。
    order = 105
    inv_peak = arb((order+1)*(2*order+1))/(6*order)
    far_floor = (order-1/(arb.pi()/order).sin())/2
    assert (104-1/(arb.pi()/104).sin())/2 > arb(105*209)/(6*104)
    assert far_floor > inv_peak and 1/inv_peak < arb(".028168")
    print("酉一阶差分的105阶锐上确界", (1/inv_peak).str(17))

    # 三类独立按面等权选取，x-y-(x-y)=0，故该字无外加钟角。
    raw = []
    for q in range(1, 106):
        a, b, c = [sum(z**q for z in group)/len(group) for group in groups]
        raw.append(1-(a*b.conjugate()*c.conjugate()).real)
    assert min(raw) > 0
    total = sum(raw)
    weights = [z/total for z in raw]
    assert arb(".00982") < weights[-2] < arb(".00983")
    assert arb(".00981") < weights[-1] < arb(".00982")
    print("无参三面闭钟字的104、105阶份额",*(w.str(17) for w in weights[-2:]))

    bands = [(arb(l), arb(r)) for l, r in angles["BANDS"]]
    d1, d2 = angles["verify_phases"](root, shell, bands)
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    modes = [(q, [w*c for c in lib["coeff"](tri, q, delta)])
             for q, w in enumerate(weights, 1)]
    opponent = sum(lib["energy"](q, c, arb("-.044"), arb("-.05")) for q, c in modes)
    lower, cells = arb(100), 0
    for sign, d in ((1, d1), (-1, d2)):
        for left, right in bands:
            for i in range(16):
                t = (left+(right-left)*i/16).union(left+(right-left)*(i+1)/16)
                for lift in (arb(0), arb.pi()):
                    a, b = sign*t+lift, sign*(t-d)+lift
                    gap = sum(lib["energy"](q, c, a, b) for q, c in modes)-opponent
                    lower = min(lower, gap.lower())
                    cells += 1
    assert lower > arb("8.66327")
    print("两支全部合窗参考味角的格数及作用减竞争点的下界", cells, lower.str(16))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
