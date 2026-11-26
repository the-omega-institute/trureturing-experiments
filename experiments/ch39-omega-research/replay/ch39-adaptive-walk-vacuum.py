"""公平游走在任意有界停时选权时的全参考味角竞争证书。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


THRESHOLD = ".05"
TARGET = (".03", ".06")

# 方差与归一模平方的线性约束乘子，整数均按 10^7 缩放。
DUALS = (
    (207045, (
        0, 0, 0, 0, 0, 0, 438258083, 0, 0, 338256182,
        316743999, 0, 518591813, 34854267, 244904920, 497813498, 0, 203423189, 376133216, 80230034,
        233283229, 11660677, 159186698, 216548271, 0, 0, 71390362, 75934588, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 30349233,
        109477811, 0, 0, 1572748, 360350502, 0, 0, 73262566, 446509503, 0,
        218693870, 55109668, 245270742, 356379480, 0, 183706672, 220796921, 132475906, 119893769, 177562938,
        0, 166012257, 85352039, 0, 36511686, 121801021, 0, 0, 0, 95241012,
        58405745, 0, 10249276, 188766714, 0, 121158494, 156523722, 0, 288868647, 0,
        142590696, 274098637, 0, 163648362, 141268846, 169244087, 49269402, 288072052, 0, 37031763,
        378539708, 2047209, 128742327, 145400045, 0, 396344521, 0, 0, 162963179, 186682840,
        0, 0, 0, 0, 0,
    )),
    (3470664, (
        0, 0, 0, 0, 0, 1185808343, 1291284244, 0, 0, 0,
        1187693066, 668269688, 0, 0, 0, 1077743699, 165846123, 0, 0, 73237737,
        691598949, 0, 0, 0, 66491526, 361716569, 0, 0, 0, 0,
        205435873, 37806071, 0, 0, 0, 54897959, 163891197, 0, 0, 0,
        0, 292925965, 5836833, 0, 0, 26683281, 400699435, 0, 0, 0,
        141689127, 428497973, 0, 0, 0, 289297110, 370293433, 0, 0, 0,
        424443343, 256961444, 0, 0, 0, 514954237, 84114228, 0, 0, 0,
        459495921, 9140270, 0, 0, 0, 403834473, 66325609, 0, 0, 0,
        550427054, 36188211, 0, 0, 177416811, 575744452, 0, 0, 0, 527032544,
        273848575, 0, 0, 148406353, 510707057, 0, 0, 0, 366311967, 0,
        0, 0, 0, 0, 0,
    )),
    (269733, (
        0, 0, 0, 0, 0, 0, 542937206, 0, 0, 22398182,
        718533470, 0, 14390892, 396768062, 0, 554644057, 19856603, 30964579, 202084643, 298400236,
        0, 54358912, 0, 137244115, 34649036, 0, 0, 0, 94340594, 0,
        0, 0, 43515469, 1353882, 0, 0, 0, 0, 35790232, 0,
        10010733, 92770545, 0, 84988404, 113204428, 0, 228691483, 90921927, 107995232, 156241933,
        240067531, 247801612, 20903634, 296635600, 258391296, 164487179, 266113342, 61413462, 239963905, 237102814,
        227823440, 8256003, 0, 469118238, 21816102, 0, 127018424, 63770317, 316724945, 0,
        0, 0, 394675207, 0, 0, 0, 206061928, 229917460, 0, 0,
        112509564, 401530299, 0, 50902363, 222051563, 68857492, 429075532, 0, 0, 421676922,
        224774991, 0, 374190301, 0, 242691943, 268628660, 0, 0, 249510184, 0,
        0, 0, 0, 0, 0,
    )),
    (3449913, (
        0, 0, 0, 0, 0, 1260116808, 1246317413, 0, 0, 0,
        1401768523, 440600371, 0, 0, 216563808, 940536611, 0, 0, 0, 339521787,
        228860774, 0, 0, 0, 141552031, 12381944, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 57028576, 50972225,
        0, 0, 0, 264329497, 13661269, 0, 0, 166538488, 297342068, 0,
        0, 0, 369737630, 139925592, 0, 0, 32019226, 357695816, 0, 0,
        0, 40181004, 195865136, 0, 0, 0, 0, 67912898, 89347051, 0,
        0, 0, 0, 223302380, 52670135, 0, 0, 0, 425616109, 44614931,
        0, 0, 106897650, 554604054, 0, 0, 0, 344525074, 469756118, 0,
        0, 0, 590655212, 163538646, 0, 0, 111780172, 406780685, 0, 0,
        0, 119515461, 43107757, 0, 0,
    )),
)



def lower_energy(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    exact = load("mincut-real-triangle-tomography")
    shell = load("three-ray-chiral-boundary")
    lib = load("mincut-global-winding")
    angles = load("monotone-winding-action")
    ctx.prec = 224
    theta = exact["unit"]((108469,109)).arg()
    bands = [(arb(l),arb(r)) for l,r in angles["BANDS"]]
    dplus, dminus = angles["verify_phases"](root,shell,bands)
    d = dplus.union(dminus)
    tri = lib["terms"](shell)
    cap = arb(23)/25*arb(2)**-20/14**2
    assert sum(m for key in ((1,0),(0,1),(1,-1)) for _,m in tri[key]) == 6
    model = load("reference-flavor-global-winding")
    model["conjugation_symmetry"](root)
    ctx.prec = 224
    modes = [lib["coeff"](tri,q,arb(0)) for q in range(1,106)]
    lower = []
    for lo, hi in bands:
        for lift in (0, 1):
            row = []
            for q,(x,y,diagonal) in enumerate(modes,1):
                # 裸面按共轭配对，可同时反转两钟。两个等味距离都在 d 内。
                assert x.imag.contains(0) and y.imag.contains(0) and diagonal.imag.contains(0)
                a, b = map(arb,TARGET)
                rival = sum(v.real*(q*t).cos() for v,t in zip((x,y,diagonal),(a,b,a-b)))
                coefficient = (x.real+y.real*acb(0,-q*d).exp())*(-1)**(q*lift)
                constant = rival-diagonal.real*(q*d).cos()-12*cap-arb(THRESHOLD)
                uppers = []
                radius = (hi-lo)/128
                for i in range(64):
                    center = lo+(2*i+1)*radius
                    z = coefficient*acb(0,q*center).exp()
                    uppers.append((z.real+radius*q*abs(z.imag)+(radius*q)**2*abs(coefficient)/2).upper())
                # 网格只在单圈能差的一个谐波上使用；其余步骤是有限的整数索引。
                row.append((constant-max(uppers)).lower())
            lower.append(row)
    return theta, lower


def verify(root):
    theta, bounds = lower_energy(root)
    costs = [max(arb(0),1-53*(q*theta)**2)**2 for q in range(1,106)]
    characters = arb_mat([[(q*j*theta).cos() for q in range(1,106)] for j in range(213)])
    margins = []
    for row, (variance, integers) in zip(bounds, DUALS, strict=True):
        beta = arb(variance)/10**7
        nu = [arb(v)/10**7 for v in integers]
        assert beta >= 0 and len(nu) == 105 and min(nu) >= 0
        moments = characters*arb_mat([[a-b] for a,b in zip(row,nu)])
        baseline = sum(c*v for c,v in zip(costs,nu))
        margin = min((moments[j,0]+beta*(j*j-212)+baseline).lower() for j in range(213))
        assert margin > arb(".26056")
        margins.append(margin)
        print("逐点对偶余量",margin.str(16))
    assert costs[0] > 0
    gap = arb(THRESHOLD)+min(margins)/105
    assert gap > arb(".05248")
    print("全部公平停时的归一能差下界",gap.str(17))


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else here)
