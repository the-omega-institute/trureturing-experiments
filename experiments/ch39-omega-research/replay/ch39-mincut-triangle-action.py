"""单位荷三角作用在最小割三窗内的双钟变分。"""
from collections import Counter, defaultdict
from itertools import combinations
from math import prod
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, arb, ctx


def triangles(shell):
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    edges = {}
    for i, j, d, z in shell["geo"]["edges"]:
        spoke = pole in (i, j)
        soft = int(spoke and (i in tuned or j in tuned))
        phase = acb(1) if spoke else z
        a, b = int(d==724 and not spoke), int(d==804 and not spoke)
        edges[i,j] = (phase, a, b, soft)
        edges[j,i] = (phase.conjugate(), -a, -b, soft)
    counts = Counter()
    coeff = defaultdict(acb)
    for i,j,k in combinations(range(19),3):
        cycle = [(i,j),(j,k),(k,i)]
        if not all(e in edges for e in cycle):
            continue
        z, a, b, softened = zip(*(edges[e] for e in cycle), strict=True)
        c, a, b = prod(z), sum(a), sum(b)
        if a<0 or (a==0 and b<0):
            a,b,c = -a,-b,c.conjugate()
        key=(a,b,sum(softened))
        counts[key]+=1
        coeff[key]+=c
    return counts, coeff


def curvature_hole(model):
    radius = arb(".10124")
    cutoff = arb("18.944146")
    _, upper = model["sign"]["density"]["eigen_band"]()
    assert upper < cutoff
    count = model["cell_count"](2, cutoff)
    pending = [(1, 0, 0)]
    included, outside = 0, 0
    while pending:
        n, i, j = pending.pop()
        a, b = [max(0, abs(2*k+1-n)-1) for k in (i,j)]
        if (a*a+b*b)*(arb.pi()/(16*n))**2 > radius**2:
            outside += 1
            continue
        # 只用四层曲率必要带的上端：若最多 15 根在 cutoff 以下，第 16 根已经过高。
        if count(n,i,j)[1] <= 15:
            included += 1
            continue
        assert n < 2**19
        pending.extend((2*n, 2*i+p, 2*j+q) for p in range(2) for q in range(2))
    print("圆盘半径、惯性排除格、圆外格", radius.str(7), included, outside)
    return radius


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    model = run_path(str(root/"ch39-mincut-three-window-angle.py"))
    ctx.prec = 192
    counts, coefficients = triangles(shell)
    assert counts == {(0,0,0):412, (0,0,1):21, (0,1,0):39, (0,1,1):3,
                      (1,0,0):39, (1,0,1):3, (1,-1,0):12}
    assert sum(counts.values()) == 529
    # 在固定高度的两钟变分中，三辐边亏额保持常数。用整段亏额统一认证。
    assert shell["soft"]["point"]["study"]["offset"] > 0 and acb.zeta_zero(1).imag > 14
    delta_cap = arb(23)/25*arb(2)**-20/14**2
    assert 0 < delta_cap < arb("4.48e-9")
    delta = arb(0).union(delta_cap)
    def coeff(a,b):
        c0 = coefficients[a,b,0]
        c1 = coefficients.get((a,b,1), acb(0))
        return (c0+c1)-delta*c1
    a,b,d = coeff(1,0), coeff(0,1), coeff(1,-1)
    print("三个工作谐波的实部和虚部", *[c.str(14) for c in (a,b,d)])
    assert a.real > arb("29.88") and b.real > arb("29.90") and d.real > arb(".0399")
    assert abs(a.imag)<arb("7.75e-10") and abs(b.imag)<arb("6.75e-10") and abs(d.imag)<arb("1e-48")
    gradient0 = ((abs(a.imag).upper()+abs(d.imag).upper())**2
                 +(abs(b.imag).upper()+abs(d.imag).upper())**2).sqrt()
    assert gradient0 < arb("1.03e-9")
    half = arb.pi()/16
    error = abs(a.imag).upper()+abs(b.imag).upper()+abs(d.imag).upper()
    minimum = arb("29.88")
    # 方框外的作用高于原点；方框内的 Hessian 是两个对角项加半正定的钟差项。
    assert minimum*(1-half.cos())-error > arb(".57")
    assert d.real*(2*half).cos()-abs(d.imag) > 0
    hessian = minimum*half.cos()-error
    assert hessian > arb("29.3058") and gradient0/hessian < arb("3.52e-11")
    radius = curvature_hole(model)
    # 逐项使用 x sin(x)≥x² sinc(half)、1-cos(x)≥x²(1-cos(half))/half²。
    # 混合项的实部贡献非负；虚部梯度误差范数至多 gradient0。
    energy = minimum*radius**2*(1-half.cos())/half**2-error
    force = minimum*radius*half.sin()/half-gradient0
    assert energy > arb(".15263") and force > arb("3.00565")
    # 相位无关的逐面相对误差 eps，使每只钟的一阶导误差至多 (42+12)eps。
    assert sum(v for (i,k,j),v in counts.items() if i) == 54
    assert sum(v for (i,k,j),v in counts.items() if k) == 54
    epsilon = arb(".03935")
    torque = 54*arb(2).sqrt()*epsilon
    assert epsilon < 1 and torque < force
    print("逐面相对误差3.935%引起的梯度误差与剩余力下界", torque.str(12), (force-torque).str(12))
    print("Hessian下界、唯一极小点至原点的上界、合窗点力和作用差下界",
          hessian.str(12), (gradient0/hessian).str(12), force.str(12), energy.str(12))

    # 与原磁拉普拉斯的三阶谱迹交叉检验，三角形每一方向均有三个循环起点。
    for n in (8236, 8359, 148341):
        height = acb.zeta_zero(n).imag
        phases = shell["geo"]["phases"](height)
        u,v = phases[724], phases[804]
        soft = 1-arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
        action = sum(soft**j*(number-(coefficients[i,k,j]*u**i*v**k).real)
                     for (i,k,j),number in counts.items())
        neutral = shell["laplacian"](0,height)
        magnetic = shell["laplacian"](1,height)
        trace = ((magnetic**3-neutral**3).trace()/6).real
        assert trace.overlaps(action)
        print("零点、三角作用", n, action.str(17))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
