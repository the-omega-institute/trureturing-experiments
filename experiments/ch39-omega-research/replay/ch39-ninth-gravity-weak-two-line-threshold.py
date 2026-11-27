"""以真空项压低第九壳的两线弱界，并在同时施加牛顿项时认证整 GeV 锐度。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    squared, weights, finite, large, bound_from = lib("ninth-weak-two-line-threshold")["verify"](root)
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [2_306_888_100_094])
    ctx.prec = 256
    planck, mz = arb("1.22089e19"), arb("91.1876")
    p0, p1 = (mz/planck)**2, (arb("1.176397e14")/planck)**2
    floor, heavy, top = (squared[i].sqrt() for i in (0, 17, 18))
    trace = 17*floor + heavy + top
    positive = [i for i,w in enumerate(weights) if w > 0]
    negative = [i for i,w in enumerate(weights) if w < 0]
    assert positive == [0,9,12,15] and len(negative) == 15
    reference = arb("-.074256577194663365545")

    # 将正重数的质量向左、负重数的质量向右推，所得有符号分布函数始终给出上界。
    events = sorted(((arb(m).lower() if g > 0 else arb(m).upper(), g) for g,_,m in field),
                    key=lambda item: (item[0], -item[1]))
    assert events[0][0] == 0 and events[0][1] < 0
    cumulative = 0
    for _,g in events:
        cumulative += g
        assert cumulative <= -20
    # E3(m/a²) 随 m 下降至零，E3(0)=1/2。按质量分部求和得单线 W(a)>=20 a⁴/(8π)。
    small = arb(".00001")
    big = arb(".05")
    assert trace > 4*small+15*big
    vacuum_upper = (sum(weights[i] for i in positive)*sum(abs(g) for g,_,_ in field)*small**4
                    + max(weights[i] for i in negative)*20*big**4)/(8*arb.pi())
    print("四条正权线均小于指定振幅时的真空上界", vacuum_upper.str(12), flush=True)
    assert vacuum_upper < arb("-1e-9")

    upper = arb(".00003850")
    count = bound_from((arb(1888)/mz)**2, upper)
    far = (large(small**2,p0,p1)-reference).upper()
    totals = []
    for high in positive:
        bounds = {i: (upper if i in (0,15) else arb(".020483945")) for i in positive}
        bounds[high] = far
        totals.append(reference + sum(weights[i]*v for i,v in bounds.items())
                      + sum(weights[i] for i in negative)*arb("-1.4e-11"))
    print("若第一与十六条均至少 1888 GeV，四种高正权线的弱差上界与覆盖段数",
          *[z.upper().str(13) for z in totals], count, flush=True)
    assert max(totals) < arb("-9e-6")

    slope, column = lib("fifth-seventh-two-pair-thermal-rotation")["thermal_columns"](field, moments, heat)
    w_a, w_b = weights[0], weights[17]
    w_floor = weights[9]+weights[12]
    w_low = weights[15]
    w_bulk = sum(weights[i] for i in negative if i != 17)
    bulk_count = len(negative)-1

    def gravity(point, q):
        a, b = point
        low = q/planck
        bulk = (trace-2*floor-low-a-b)/bulk_count
        terms = [(w_a,a), (w_b,b), (w_floor,floor), (w_low,low), (w_bulk,bulk)]
        readings = [column(v) for _,v in terms]
        residual = arb_mat([[sum(term[0]*row[k] for term,row in zip(terms,readings,strict=True))-target]
                            for k,target in enumerate((1,arb("1e-122")))])
        da, db, dc = slope(a), slope(b), slope(bulk)
        jac = arb_mat([[w_a*da[k]-w_bulk*dc[k]/bulk_count, w_b*db[k]-w_bulk*dc[k]/bulk_count]
                       for k in range(2)])
        return residual, jac, bulk

    center = [arb(".14376178"), arb(".18300802")]
    radius = arb("1e-6")
    box = [x+arb(0,radius) for x in center]
    q_range = arb(1887).union(arb(1888))
    residual, center_jac, _ = gravity(center, q_range)
    _, jac, bulk = gravity(box, q_range)
    pre = center_jac.mid().inv().mid()
    assert not pre.det().contains(0)
    defect = arb_mat([[1,0],[0,1]]) - pre*jac
    ratios = [(v.abs_upper()/radius+sum(x.abs_upper() for x in row))
              for v,row in zip((pre*residual).entries(),defect.tolist(),strict=True)]
    print("全能标段引力二式闭盒自映比", [z.str(12) for z in ratios],flush=True)
    assert max(ratios) < arb(".5")
    assert floor < q_range/planck < bulk < box[0] < box[1] < heavy < top
    # 总迹精确，且所有振幅在源谱上下界内。平移后的源谱秩二，故满足 Schur–Horn 主序。
    # 对每个 q，统一压缩映射给出连续的唯一引力根。弱差的面界对整个引力根盒有效。
    def face(q):
        low_correction = finite((arb(q)/mz)**2)-reference
        return (reference + w_floor*(finite(squared[0]/p0)-reference) + w_low*low_correction
                + sum(w*(large(a.lower()**2,p0,p1)-reference)
                      for w,a in ((w_a,box[0]),(w_b,box[1]),(w_bulk,bulk))))

    left, right = face(1887), face(1888)
    print("连续引力精确解在 1887、1888 GeV 的弱差", left.str(16), right.str(16),flush=True)
    assert left > arb("4e-6") and right < arb("-9e-6")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
