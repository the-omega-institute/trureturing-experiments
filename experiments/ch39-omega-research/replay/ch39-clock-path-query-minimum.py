"""固定路径叠加在十一枚上实现最优单次通过率的最少双钟查询数。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

root = Path(__file__).parent
selection = run_path(str(root / "ch39-single-pass-arc-selection.py"))
flag_modes = ((0,-363), (10,-353), (36,327), (55,308), (74,289), (94,269),
              (102,261), (151,212), (159,204), (167,196), (170,193))
separator = (53102423,-518738,-32119013,1978455,3953026,-2014649,3057618,1433243,-384238,-1438596)
label_modes = (
    (-1139, -932), (-995, -1016), (-802, -1318), (-757, -767), (-587, -1474), (-581, -492),
    (0, 1610), (0, 1706), (0, 1778), (0, 1841), (0, 1990), (0, 2126),
    (31, -2119), (38, 2112), (66, -2084), (153, 1997), (166, 1594), (181, 1969),
    (210, 1940), (257, 1893), (279, 1871), (356, 1794), (364, 1786), (385, 1593),
    (399, 1751), (406, 1744), (426, -1724), (454, 1696), (462, 1688), (517, 1633),
    (547, 1603), (558, 1592), (568, 710), (596, 1554), (604, 1546), (622, 1356),
    (629, 1521), (647, -1503), (707, 1443), (717, 1433), (798, 1352), (810, 1340),
    (844, -1306), (858, 1292), (874, 1276), (960, 1190), (1011, 879), (1035, 921),
    (1035, 1115), (1119, 1031), (1142, 1008), (1278, -385), (1864, -286), (2023, -127),
    (2119, -31), (2126, -24),
)


def phase_columns(angles, pairs, modes):
    return arb_mat([[1]*len(modes), *[
        [(u*(angles[i][0]-angles[j][0])+v*(angles[i][1]-angles[j][1])).cos()
         for u, v in modes] for i, j in pairs]])


def probability_weights(matrix):
    target = arb_mat([[int(i == 0)] for i in range(matrix.nrows())])
    p = matrix.solve(target).entries()
    assert all(v > 0 for v in p)
    return p


def positive_effect(angles, rates, others):
    p = probability_weights(phase_columns(angles, [(1,i) for i in others], flag_modes))
    pairs = list(combinations(others, 2))
    entries = (phase_columns(angles, pairs, flag_modes)*arb_mat([[v] for v in p])).entries()[1:]
    gram = dict(zip(pairs, entries, strict=True))
    # 检验 G 减去目标通过概率对角阵正定，故存在有界二元测量效应。
    remainders = [[1-rates[i] if i == j else gram[min(i,j),max(i,j)] for j in others] for i in others]
    for k in range(1, len(others)+1):
        assert arb_mat([r[:k] for r in remainders[:k]]).det() > 0
    assert max(abs(u)+abs(v) for u,v in flag_modes) == 363


def exclude_shorter(angles, rates, others):
    budget = 362
    grid = 2**44
    def table(angle):
        pairs = [(grid*(k*angle).cos(), grid*(k*angle).sin()) for k in range(-budget,budget+1)]
        return [(int(c.floor().unique_fmpz()), int(s.floor().unique_fmpz())) for c,s in pairs]

    # 先将一维三角函数向下读为整数，二维逐个频率只需整数乘法。
    x_tables, y_tables = [], []
    for i in others:
        x_tables.append(table(angles[i][0]-angles[1][0]))
        y_tables.append(table(angles[i][1]-angles[1][1]))
    lowest = None
    for u in range(budget+1):
        row = []
        for tab, factor in zip(x_tables, separator, strict=True):
            c, s = tab[budget+u]
            row.append((c*factor, s*factor))
        for v in range(u-budget,budget-u+1):
            score = sum(c*tab[budget+v][0]-s*tab[budget+v][1]
                        for (c,s),tab in zip(row,y_tables,strict=True))
            if lowest is None or score < lowest:
                lowest = score
    total = sum(abs(k) for k in separator)
    minimum = arb(lowest-total*(4*grid+2))/grid**2
    allowed = sum(abs(k)*rates[i].sqrt() for i,k in zip(others,separator,strict=True))
    margin = (minimum-allowed)/10**8
    assert margin > arb(".00035")
    print("362 步频域整数分离证书的净裕量",margin.str(14))


def verify():
    joint = selection["entropy"]["joint"]
    rows, times, *_ = joint["columns_at"](selection["scale"], selection["ranges"]["indices"])
    weights = joint["solve"](rows[:8], times)+[arb(0)]*3
    geometry = run_path(str(root / "ch39-arc-relative-determinant.py"))
    ctx.prec = 384
    xs = [row[1] for row in rows]
    clocks = [geometry["phases"](selection["scale"]/x.sqrt()) for x in xs]
    angles = [[z[d].arg() for d in (724,804)] for z in clocks]
    rates = [p*xs[1]/(x*weights[1]) for p,x in zip(weights,xs,strict=True)]
    others = [i for i in range(len(rows)) if i != 1]
    positive_effect(angles,rates,others)
    exclude_shorter(angles,rates,others)
    label_weights = probability_weights(phase_columns(angles,list(combinations(range(11),2)),label_modes))
    assert len(label_weights) == 56 and max(abs(u)+abs(v) for u,v in label_modes) == 2150
    print("最优放行的固定路径查询数、完整标签辨认的构造查询数",363,2150)


if __name__ == "__main__":
    verify()
