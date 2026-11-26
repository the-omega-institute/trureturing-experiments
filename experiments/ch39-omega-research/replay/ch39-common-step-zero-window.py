"""固定两档谱规则的零点窗口。钟边只跨三个内点，带阈惯性由三阶 Schur 补确定。"""
from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).parent
audit = run_path(str(ROOT/"ch39-monotone-heavy-mass-ratio.py"))
geometry = audit["SHELL"]
charges = (1,4,2,3,6)
clock_edges = [(a,b,d,u) for a,b,d,u in geometry["edges"] if d in (724,804)]
core = [i for i in range(19) if i not in (6,7,8)]
inner = [6,7,8]


def sub(matrix, rows, cols):
    return acb_mat([[matrix[i,j] for j in cols] for i in rows])


def supplier(q, level):
    matrix = audit["cut_matrix"](q)
    for i in range(19):
        matrix[i,i] -= level
    fixed = sub(matrix,core,core)
    fixed_positive = audit["count_above"](fixed,arb(0))
    inverse = fixed.inv()
    cross = sub(matrix,core,inner)
    small = sub(matrix,inner,inner)
    inserts = [(core.index(a), inner.index(b), d, -u**q) for a,b,d,u in clock_edges]
    def count(z):
        off = acb_mat(cross)
        for i,j,d,v in inserts:
            off[i,j] = v*z[d]**q
        schur = small - off.conjugate().transpose()*inverse*off
        # 总数 19，剔除最低的一根，所以 18 减去高于分界的根数。
        return 18-fixed_positive-audit["count_above"](schur,arb(0))
    return count


def verify():
    ctx.prec = 160
    target = (2,1,0,0,1)
    neutral = audit["cut_matrix"](0)
    for a,b,_,_ in clock_edges:
        neutral[a,b] = neutral[b,a] = -1
    eps = arb("1e-7")
    levels = [arb("7.7")-eps,arb("7.7")+eps,arb("7.66")-eps,arb("7.72")+eps]
    for t in levels:
        assert audit["count_above"](neutral,t) == 18
    counters = [[supplier(q,t) for q in charges] for t in levels]
    fitting, ramp = [], []
    types = Counter()
    zeros = acb.zeta_zeros(10855,1238)
    for index,zero in enumerate(zeros,10855):
        assert zero.imag > 10000
        z = geometry["phases"](zero.imag)
        below,above = [tuple(count(z) for count in row) for row in counters[:2]]
        assert below == above
        types[below] += 1
        if below == target:
            fitting.append(index)
            low, high = [tuple(count(z) for count in row) for row in counters[2:]]
            if low == high == target:
                ramp.append(index)
    assert len(zeros) == 1238 and sum(types.values()) == 1238
    assert (len(types),len(fitting),len(ramp)) == (30,17,4) and 11431 in ramp
    print("两档边界 7.70 所得计数组合数、目标命中数、整个斜坡间隙的保证数",
          len(types),len(fitting),len(ramp))
    print("硬分界与完整斜坡间隙命中编号", fitting, ramp)
    # 任意两只相位各偏离目标不足 1/1000 时，六种荷的扰动都小于 0.02。
    assert 2*arb(2).sqrt()*6/1000 + 2*eps < arb(".02")
    print("全谱命中的面积密度下界", (1/(1000*arb.pi())**2).str(14))


if __name__ == "__main__":
    verify()
