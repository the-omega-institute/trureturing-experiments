"""第六壳字符相交的矩阵铅笔与 CKM 顶底夹角给出共同三代平面的位移障碍。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


RADIUS = ".06709"
HEAVY_RADIUS = ".0169"


def ray_bounds(probs, spectrum, factor, shifts):
    inverse = arb_mat([[1/(v+s) for s in shifts] for v in spectrum])
    denominator = probs*inverse
    return [max([arb(0)]+[(factor/denominator[i, j]-s).lower() for j, s in enumerate(shifts)])
            for i in range(19)]


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    ctx.prec = 256
    number = 2_306_888_100_091
    table = [line.split() for line in (root/"ch39-sixth-quark-distance-cells.tsv").read_text().splitlines()]
    scale = 10**14
    assert len(table) == 82 and table[:4] == [["scale", str(scale)], ["shell", str(number)],
            ["sites", "14", "8", "7", "6"], ["dw_upper", "281000", str(scale)]]
    # 只使用既定训练谱和第六壳。原复算表给出所有重叠概率的向上有理舍入。
    source = 2_306_888_100_086
    _, train, _ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root, ((source, ".0294"),), ".1", (), ((6, 7, 10, 13), (6, 7, 10, 12)))
    rows = lib("free-four-site-native-positive-spectrum")["native"](root, number)
    ctx.prec = 256
    radius = arb(RADIUS)
    training_radius = arb("2.81e-9")
    shifts = [arb(2)**n*k/64 for n in range(-36, 7) for k in range(64, 128)]
    coords = lib("independent-unitary-quark-right-holdout")["coordinates"]
    character = []
    heavy_bounds = []
    for iq, (q, power) in enumerate(((-4, 2), (2, 1))):
        data = table[4+39*iq:4+39*(iq+1)]
        source_values = train[source, q]
        assert data[0] == ["train", str(q), *[str(int((v*scale).lower().floor().unique_fmpz()))
                                              for v in source_values]]
        assert all(a[:3] == ["prob", str(q), str(i)] and len(a) == 22 for i, a in enumerate(data[1:20]))
        probs = arb_mat([[arb(x)/scale for x in line[3:]] for line in data[1:20]])
        delta = power*(2*training_radius+radius)
        factor = (1-delta**2/2)**2
        cap = arb(1) if iq == 0 else arb(".3")
        assert all(0 < v < cap for v in source_values) and delta < 1
        lower = ray_bounds(probs, source_values, factor, shifts)
        reverse = ray_bounds(probs, [cap-v for v in source_values], factor, shifts)
        upper = [cap-v for v in reverse]
        chars = [coords(row)[4:] for row in rows[iq]]
        character.append([sum(w*col[j] for w, col in zip(lower, chars, strict=True)).lower().union(
                          sum(w*col[j] for w, col in zip(upper, chars, strict=True)).upper())
                          for j in range(2)])
        # 两个重窗也可单独检验。较小闭球中，分别选平面时每个重窗仍全窗可达。
        small_delta = power*(2*training_radius+arb(HEAVY_RADIUS))
        small_factor = (1-small_delta**2/2)**2
        low_heavy = ray_bounds(probs, source_values, small_factor, shifts)
        up_heavy = [cap-v for v in ray_bounds(probs, [cap-v for v in source_values], small_factor, shifts)]
        heavy_bounds.append([sum(w*col[j] for w, col in zip(low_heavy, chars, strict=True)).lower()
                             for j in range(2)])
        heavy_lo, heavy_hi = map(arb, ((".8573", ".8659"), (".01392", ".01410"))[iq])
        upper_matrix = sum((w*p for w, p in zip(up_heavy, rows[iq], strict=True)), acb_mat(4, 4))
        third = sorted((z.real for z in upper_matrix.eig(algorithm="rump")), key=lambda z: z.mid())[2]
        fourth_lower = max(heavy_bounds[iq])
        assert third < heavy_lo and fourth_lower > heavy_hi
        print("荷、可覆盖重窗的第三值上界与第四值下界", q, third.upper().str(12), fourth_lower.str(12))
        if iq == 0:
            trivial = [acb_mat([[p[i, j] for j in range(2)] for i in range(2)]) for p in rows[iq]]
            light_lower = sum((v*h for v, h in zip(lower, trivial, strict=True)), acb_mat(2, 2))

    # 两个字符按上型小到大编号，同一编号贯穿两个荷。
    order = sorted(range(2), key=lambda i: character[0][i].mid())
    ua, ub = [character[0][i] for i in order]
    da, db = [character[1][i] for i in order]
    assert ua < ub and db < da
    lua, lub = [heavy_bounds[0][i] for i in order]
    lda, ldb = [heavy_bounds[1][i] for i in order]
    assert lua < arb(".8659") < lub and ldb < arb(".01410") < lda
    # 交线中同一个单位向量在字符 a 上的概率 p 须同时满足两个 Rayleigh 上界。
    p_lower = (lub-arb(".8659"))/(lub-lua)
    p_upper = (arb(".01410")-ldb)/(lda-ldb)
    print("共平面重窗所需的同一字符概率下界、上界", p_lower.lower().str(14), p_upper.upper().str(14))
    assert p_lower > arb(".13250193") and p_upper < arb(".13241985")
    print("两枚字符的上、下型区间", [(v.lower().str(10), v.upper().str(10)) for v in (ua, ub, da, db)])
    t = (ub-ua)/(da-db)
    # 两个表达式给同一个字符交点。正的 t 也给出交点 >= ub。
    crossing_lower = ub.lower()
    xmax, charm_lo, charm_hi, bottom_lo, bottom_hi = map(arb, (".00000629", ".00306", ".00316", ".01392", ".01410"))
    assert crossing_lower > charm_hi+t.upper()*bottom_hi
    theta13, theta23 = arb(".00377"), arb(".04275")
    ctb = theta13.cos()*theta23.cos()
    sin_tb = (1-ctb**2).sqrt()
    sin_et = t.upper()*bottom_hi*sin_tb/(crossing_lower-charm_hi-t.upper()*bottom_hi)
    sin_eb = sin_et+sin_tb
    mean_lower = bottom_lo*(1-sin_eb**2)
    assert db < mean_lower and bottom_hi < da
    leakage = (bottom_hi*sin_eb)**2/((mean_lower-db.upper())*(da.lower()-bottom_hi))
    # 上窗存在则上型平凡二阶块的最低本征值不超过它的窗顶。
    vals = sorted((v.real for v in light_lower.eig(algorithm="rump")), key=lambda z: z.mid())
    lo, hi = vals
    assert 0 < lo < xmax < hi and ua.lower() > xmax
    ident = acb_mat([[1, 0], [0, 1]])
    ground = (ident*hi-light_lower)/(hi-lo)
    sin2 = ((xmax-lo)/(hi-lo)).upper()
    assert 0 < sin2 < 1
    ranks = []
    tiny = arb(0)
    for h in trivial:
        tr = h.trace().real
        if tr < arb("1e-45"):
            tiny += abs(tr).upper()
            continue
        # 每个非零二阶项是秩一的实际四站谱投影压缩。
        share = (ground*h).trace().real/tr
        assert 0 < share < 1
        gap = ((1-sin2)*share).sqrt()-(sin2*(1-share)).sqrt()
        assert gap > 0
        ranks.append((gap**2).lower())
    assert len(ranks) == 7
    trace_cap = (lo+hi+(xmax-lo)/min(ranks)+tiny).upper()
    total = trace_cap+ub.upper()*leakage
    print("铅笔系数上界、字符泄漏平方、平凡块迹上界、粲质量上界",
          *[v.upper().str(13) for v in (t, leakage, trace_cap, total)], flush=True)
    assert t < arb("6.404453") and leakage < arb(".000649316")
    assert trace_cap < arb(".002017735") and total < arb(".002628008") < charm_lo
    return total


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
