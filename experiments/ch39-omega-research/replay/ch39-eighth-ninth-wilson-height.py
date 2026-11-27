"""同一最短弧读出上更换基本闭路相位，及以真实连续高度重算色支撑指标。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def read_weights(shell, gram, h, phases, *, with_gap=False):
    matrix = shell["laplacian"](6, h, phases=phases)
    values, vecs = matrix.eig(right=True, algorithm="rump")
    order = sorted(range(19), key=lambda i: values[i].real.mid())
    assert all(v.imag.contains(0) for v in values)
    assert all(values[i].real < values[j].real for i, j in zip(order, order[1:]))
    near = shell["boundary_sites"]()
    pairs = {i:j for i,j,length,_ in shell["geo"]["edges"] if length == 12}
    pole = shell["soft"]["pole"]
    weights = []
    for j in order:
        norm = sum(abs(vecs[i,j])**2 for i in range(19))
        read = [vecs[pole,j]]
        for i in near:
            b = pairs[i]
            link = -matrix[i,b]
            read.append((7*vecs[i,j]+24*link*vecs[b,j])/25)
        source = sum(x*y for x,y in zip(gram,read,strict=True))
        weights.append((source*read[2].conjugate()/(gram[2]*norm)).real)
    assert sum(weights).contains(1)
    if with_gap:
        return weights, min(values[j].real-values[i].real for i,j in zip(order,order[1:]))
    return weights


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    shell = lib("three-ray-chiral-boundary")
    geo = shell["geo"]
    native = lib("native-rh-ray-pvm")
    arc = lib("eighth-shortest-arc-pivot-support")
    support = lib("cross-shell-pivot-support")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    ctx.prec = 320
    base = 2_306_888_100_093
    h0, h1 = (acb.zeta_zero(base+i).imag for i in (0,1))
    span = h1-h0
    clocks = [geo["phases"](h) for h in (h0,h1)]
    gram = native["gram_cholesky"](root)[2]
    pole = shell["soft"]["pole"]
    for d in (724,804):
        edge = next((i,j) for i,j,dd,_ in geo["edges"] if dd == d)
        i,j = edge
        loops = []
        for h in (h0,h1):
            matrix = shell["laplacian"](6,h)
            phi = -matrix[pole,i]*matrix[i,j]*matrix[j,pole]
            loops.append(phi/abs(phi))
        ratio = loops[1]*loops[0].conjugate()
        expected = (clocks[1][d]*clocks[0][d].conjugate())**6
        assert ratio.overlaps(expected)
        angle = ratio.arg()
        assert 0 < angle < arb.pi()
        print("平方弦长、基本闭路、两高度间的通量增角", d, (pole,i,j), angle.str(15), flush=True)
    def score(table, weights):
        box = arc["score"](table, weights)
        assert box is not None
        return tuple(Fraction(int(x),arc["SCALE"]) for x in box)
    tables = [support["projected_tables"](core,root,base+i,None) for i in (0,1)]
    for a,b,sgn in ((0,0,1),(1,0,1),(0,1,-1),(1,1,-1)):
        phases = {724:clocks[a][724], 804:clocks[b][804]}
        w = read_weights(shell, gram, h0, phases)
        bounds = [score(table,w) for table in tables]
        assert all((l>0 if sgn>0 else r<0) for l,r in bounds)
        print("724、804 时钟（0/1 为原/新）、两张端点表的指标",a,b,*bounds,flush=True)
    # 固定端点表时，用谱隙与投影扰动界证明该窄窗内 S3 连续。
    mid = h0 + arb(".358235") * span
    halfwidth = arb(".000005") * span
    weights, gap = read_weights(shell, gram, mid, geo["phases"](mid), with_gap=True)
    speeds = [arb(0)] * 19
    for i,j,d,_ in geo["edges"]:
        if d in (724,804):
            norm = 6 * geo["a" if d==724 else "b"]
            speeds[i] += norm
            speeds[j] += norm
    offset = shell["soft"]["point"]["study"]["offset"]
    assert offset >= 0 and h0 > arb("1e11")
    # 三条软辐边的权导数各 <10^-35；含对角的绝对行和最多六倍。
    motion = (max(speeds)+arb("6e-35")) * halfwidth
    probe_norm = sum(c*c for c in gram).sqrt()/gram[2]
    error = probe_norm * motion/(gap-motion)
    major = (0,7,16)
    separation = min(abs(weights[i]).lower() for i in major) - max(
        abs(w).upper() for i,w in enumerate(weights) if i not in major
    )
    assert gap > 2*motion and separation > 2*error
    print("窄窗中心相邻谱隙、算子位移、单权误差、前三截断间隙",
          *[x.str(12) for x in (gap,motion,error,separation)], flush=True)
    for s, expected in (("0.35823",1),("0.35824",-1)):
        h = h0+arb(s)*span
        w = read_weights(shell, gram, h, geo["phases"](h))
        table = support["projected_tables"](core,root,base,None,height=h)
        l,r = score(table,w)
        for fixed in tables:
            fl,fr = score(fixed,w)
            assert fl>0 if expected>0 else fr<0
        assert (l>0 if expected>0 else r<0)
        print("高度分数、实际高度和本高度全量响应的 S3 包围",s,h.str(29),l,r,flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
