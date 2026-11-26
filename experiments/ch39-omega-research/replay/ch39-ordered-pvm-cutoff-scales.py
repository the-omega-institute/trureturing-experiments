"""中心实际零点的五荷三枢轴，在三种保序标尺下同时匹配牛顿与弱差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


REFERENCE = ("-.074256577194663365546", "-.074256577194663365544")


def symmetric(radius):
    return arb(0, radius.upper())


def weak_radius(x):
    return arb("2.398e-21") + arb("7.468e-36")/x


def normalized_positions(levels, family):
    scale = [x.log() for x in levels] if family == "幂律" else levels
    return [arb(0), *((z-scale[0])/(scale[-1]-scale[0]) for z in scale[1:-1]), arb(1)]


def exponential_moment(weights, distances, low, high=None):
    # b 在 [low,high] 内，None 代表上端无限。末点的物理位置总是 v。
    total = weights[-1]
    for c, r in zip(weights[:-1], distances[:-1], strict=True):
        left = arb(0) if high is None else (-r*high).exp()
        total += c * left.union((-r*low).exp())
    return total


def verify(root):
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    _, fields = run_path(str(root / "ch39-three-index-free-atoms.py"))["spectra"](root)
    heat = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    weak_column = run_path(str(root / "ch39-full-neutral-free-cutoff.py"))["weak_column"]
    rows = fields[0]
    ctx.prec = 256
    p0, p1 = [(arb(m)/arb("1.22089e19"))**2 for m in ("91.1876", "1.176397e14")]
    ref = arb(REFERENCE[0]).union(arb(REFERENCE[1]))
    weak_series = native["weak_series"]()
    v_low, v_high = arb("1e-5"), arb(".99")
    u_high = 1000*p0
    u_small = p0/200
    lip = sum(abs(a) for _, a, _ in rows)/(2*arb.pi())
    slope = sum(a for _, a, _ in rows)/(2*arb.pi())
    error = sum(abs(a)*(1-heat(2,z)) for _, a, z in rows)/(2*arb.pi())
    assert lip < 33 and slope < 0 and error < arb(".0002"), (lip,slope,error)
    k_high = v_high*sum(a*heat(2,z/v_high) for _, a, z in rows)/(2*arb.pi())
    print("牛顿项绝对斜率界、一次近似斜率、统一误差、K(0.99)",
          lip.str(12), slope.str(12), error.str(12), k_high.str(15))
    counts = {family: [0, 0] for family in ("自由", "幂律", "指数")}
    for q, levels, direct, _ in native["shell_weights"](root):
        # 过原点的非负凹标尺，端点比至多为本征值比；也检验小幂次。
        ratio = levels[-1]/levels[0]
        assert ratio < 306 and 61/ratio.log() > arb("10.65")
        for pivot, projected in enumerate(direct, 1):
            weights = [v.real for v in projected]
            tv = sum(abs(c) for c in weights)
            assert sum(weights).contains(1) and tv < 712
            # 如果 u<p0，直接用 p0。否则对所有十九项使用大截止余量。
            u_limit = max(p0, tv*arb("7.468e-36")/(-ref-tv*arb("2.398e-21")))
            span_min = (1/(lip*tv*u_limit)).log()
            assert span_min > 61
            assert tv*lip*v_low < 1

            tails = [sum(weights[j:]) for j in range(1, len(weights))]
            cut = min(range(len(tails)), key=lambda j: float(tails[j]))
            tail = tails[cut]
            assert all(tail < s for j, s in enumerate(tails) if j != cut)
            assert tail < arb("-.91")
            floor_ok = pivot > 1 or q == 6
            a = 1-tail
            small_d = weak_column(u_small, p0, p1)
            large_v_d = ref+symmetric(weak_radius(v_low))
            assert a*small_d + tail*large_v_d > 0
            assert a*weak_series(arb(1000)) + tail*large_v_d < 0
            assert a*symmetric(lip*u_high) + tail*k_high > 1
            d_at_mz = weak_series(arb(1))
            if floor_ok:
                assert a*d_at_mz + tail*large_v_d > 0
            else:
                # d=-D 严格递增。后缀分部求和给所有有序点上的上界。
                assert a*d_at_mz + tail*weak_column(arb(1),p0,p1) < 0
            counts["自由"][0] += 1
            counts["自由"][1] += int(floor_ok)
            print("荷、枢轴、最负后缀首层和深度",q,pivot,cut+2,(-tail).str(9),end="")

            for family in ("幂律", "指数"):
                positions = normalized_positions(levels, family)
                distances = [1-t for t in positions]
                no_solution = (pivot == 3 and q in (1,-3)) or (
                    family == "指数" and (q,pivot) == (6,1))
                if no_solution:
                    # 所有可能的对数带宽 b>61。对长尾保持最后谱权的常数项。
                    bounds = []
                    for k in range(2*61,2*384):
                        moment = exponential_moment(weights, distances, arb(k)/2, arb(k+1)/2)
                        bounds.append((slope*moment).upper())
                    moment = exponential_moment(weights, distances, arb(384))
                    bounds.append((slope*moment).upper())
                    bound = max(arb(0), max(bounds)) + tv*error
                    assert bound < 1, (q,pivot,family,bound)
                    print("，",family,"排除，K 上界",bound.upper().str(6),end="")
                    continue

                lower = p0 if pivot > 1 else u_small
                first = weights[0]
                rest = sum(abs(c) for c in weights[1:])
                assert first > 1
                # 第二层已远离 p0，全部高层的误差共用这个下界。
                for u, d, sign in ((lower,d_at_mz if pivot>1 else small_d,1),
                                   (u_high,weak_series(arb(1000)),-1)):
                    x2 = ((1-positions[1])*u.log()+positions[1]*v_low.log()).exp()
                    assert x2 > p0
                    assert sign*(first*d + (1-first)*ref + symmetric(rest*weak_radius(x2))) > 0

                b0, b1 = (v_high/u_high).log(), (v_high/lower).log()
                bounds = []
                for k in range(128):
                    left, right = [b0+(b1-b0)*j/128 for j in (k,k+1)]
                    moment = exponential_moment(weights, distances, left, right)
                    bounds.append((v_high*slope*moment-tv*error).lower())
                lower_k = min(bounds)
                assert lower_k > 1, (q,pivot,family,lower_k)
                # u>=p0 时所有单点贡献 D(u)<=D(p0)。结合 b>61 夹住其余十八层。
                if pivot == 1:
                    next_floor = p0*(61*positions[1]).exp()
                    assert first*d_at_mz + (1-first)*ref + rest*weak_radius(next_floor) < 0
                if (q,pivot,family) == (1,2,"幂律"):
                    c_lo, c_hi = arb("4.1463246"), arb("4.1463248")
                    v_lo, v_hi = arb(".126"), arb(".129")
                    for c, sign in ((c_lo,1),(c_hi,-1)):
                        x2 = ((1-positions[1])*(c*p0).log()+positions[1]*v_lo.log()).exp()
                        assert x2 > p0
                        weak = first*weak_series(c)+(1-first)*ref+symmetric(rest*weak_radius(x2))
                        assert sign*weak > 0
                    for v, sign in ((v_lo,-1),(v_hi,1)):
                        low = (v/(c_hi*p0)).log()
                        high = (v/(c_lo*p0)).log()
                        moment = exponential_moment(weights, distances, low, high)
                        assert sign*(v*slope*moment+symmetric(tv*error)-1) > 0
                    gap = (levels[-1]/levels[0]).log()
                    nu = ((v_lo/(c_hi*p0)).log()/gap).union((v_hi/(c_lo*p0)).log()/gap)
                    assert arb("25.385") < nu < arb("25.394")
                    qmin = arb("91.1876")*c_lo.union(c_hi).sqrt()
                    assert arb("185.68098") < qmin < arb("185.68100")
                counts[family][0] += 1
                counts[family][1] += int(pivot>1)
                print("，",family,"存在，牛顿上端面下界",lower_k.str(6),end="")
            print(flush=True)
    assert counts == {"自由":[15,11], "幂律":[13,8], "指数":[12,8]}
    print("允许任意正位置、再要求 x1≥p0 的可行项数", counts)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
