"""固定重复 W 的第六壳顶底共同可达三角形、锐下界与窗顶最小放宽量。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def blocks(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    number = 2_306_888_100_091
    reference = lib("common-embedding-central-right-ckm")
    source = 2_306_888_100_086
    _, train, _ = reference["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                            ((6, 7, 10, 13), (6, 7, 10, 12)))
    projectors = lib("free-four-site-native-positive-spectrum")["native"](root, number)
    coords = lib("independent-unitary-quark-right-holdout")["coordinates"]
    ctx.prec = 256
    table = [line.split() for line in (root/"ch39-sixth-quark-distance-cells.tsv").read_text().splitlines()]
    scale = 10**14
    assert len(table) == 82 and table[:4] == [["scale", str(scale)], ["shell", str(number)],
            ["sites", "14", "8", "7", "6"], ["dw_upper", "281000", str(scale)]]
    char = []
    trivial_blocks = []
    for iq, (charge, power) in enumerate(((-4, 2), (2, 1))):
        data = table[4+39*iq:4+39*(iq+1)]
        weights = train[source, charge]
        assert data[0] == ["train", str(charge), *[str(int((w*scale).lower().floor().unique_fmpz()))
                                                   for w in weights]]
        cap = arb(1) if iq == 0 else arb(".3")
        assert all(0 < w < cap for w in weights)
        # 表中的 p 是相干中心重叠的上界。对 cap-w 再用一次它，给出同一中心读数的下界。
        error = 2*power*arb("2.81e-9")*(max(w.upper() for w in weights)-min(w.lower() for w in weights))
        amplitudes = []
        for i, row in enumerate(data[1:20]):
            assert row[:3] == ["prob", str(charge), str(i)] and len(row) == 22
            prob = [arb(n)/scale for n in row[3:]]
            lower = cap-sum(p*(cap-w) for p, w in zip(prob, weights, strict=True))
            upper = sum(p*w for p, w in zip(prob, weights, strict=True))
            amplitudes.append(lower.lower().union(upper.upper())+arb(0, error.upper()))
        four = sum((a*p for a, p in zip(amplitudes, projectors[iq], strict=True)), acb_mat(4, 4))
        values = coords(four)[4:]
        trivial = acb_mat([[four[i, j] for j in range(2)] for i in range(2)])
        trivial_blocks.append(trivial)
        top_trivial = max(v.real.upper() for v in trivial.eig(algorithm="rump"))
        assert top_trivial < min(v.lower() for v in values)
        print("荷、平凡块最高值上界、两字符区间", charge, top_trivial.str(12),
              [(v.lower().str(13), v.upper().str(13)) for v in values], flush=True)
        char.append(values)

    return trivial_blocks, char


def verify(root):
    trivial_blocks, char = blocks(root)
    order = sorted(range(2), key=lambda i: char[0][i].mid())
    ua, ub = [char[0][i] for i in order]
    da, db = [char[1][i] for i in order]
    assert ua < ub and db < da
    tx0, tx1 = arb(".8573"), arb(".8659")
    by0, by1 = arb(".01392"), arb(".01410")
    assert ua < tx0 < tx1 < ub and db < by0 < by1 < da

    def frontier(x):
        prob = (ub-x)/(ub-ua)
        return prob, db+prob*(da-db)

    for x in (tx0, tx1):
        p, y = frontier(x)
        print("顶窗端点、字符 a 概率、最小底耦合", x.str(7),
              [(v.lower().str(14), v.upper().str(14)) for v in (p, y)], flush=True)
    bottom_min = frontier(tx1)[1]
    widening = bottom_min-by1
    assert arb(".01652462565") < bottom_min < arb(".01652464838")
    assert arb(".00242462565") < widening < arb(".00242464838")
    first = frontier(tx0)[1]
    assert arb(".01844922621") < first < arb(".01844924967")
    # 取一个三角形内部点，在同一法向量上重放两个压缩特征值方程。
    x = (ua+2*ub)/3
    y = (db+2*da)/3
    alpha, gamma = (ub-x)/(x-ua), (da-y)/(y-db)
    identity = acb_mat([[1, 0], [0, 1]])
    ru, rd = [(identity*v-h).inv()[0, 0].real for v, h in zip((x, y), trivial_blocks, strict=True)]
    beta, delta = (ub-x)*ru, (da-y)*rd
    na = (gamma*beta+delta)/(1-alpha*gamma)
    nb = alpha*na+beta
    assert na > 0 and nb > 0
    assert (na/(x-ua)+nb/(x-ub)+ru).contains(0)
    assert (na/(y-da)+nb/(y-db)+rd).contains(0)
    # 两条重特征线的夹角若固定为非零，Rayleigh 商缺口还满足一个统一正下界。
    ctb = arb(".00377").cos()*arb(".04275").cos()
    denominator = (ub-ua)/((tx0-ua)*(da-db))+1/(bottom_min-db)
    mixing_gap = (1-ctb**2)/denominator
    print("目标顶底夹角另迫使底耦合高于斜边至少", mixing_gap.lower().str(13), flush=True)
    assert mixing_gap > arb(".0000221267")
    print("顶窗固定时底窗上端最小放宽量", widening.lower().str(14), widening.upper().str(14), flush=True)
    return bottom_min, widening


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
