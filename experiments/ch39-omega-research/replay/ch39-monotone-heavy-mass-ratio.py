"""两荷十八重层在一切共同单调质量读数下的锐界，以及交点模型的独立质量比。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


ROOT = Path(__file__).parent
SHELL = run_path(str(ROOT / "ch39-arc-relative-determinant.py"))


def cut_matrix(charge):
    # 保留十二条钟边对顶点度的贡献，仅删除其非对角元。
    matrix = acb_mat(19,19)
    for i,j,d,u in SHELL["edges"]:
        matrix[i,i] += 1
        matrix[j,j] += 1
        if d not in (724,804):
            matrix[i,j] = -u**charge
            matrix[j,i] = matrix[i,j].conjugate()
        else:
            assert i in range(6) and j in (6,7,8)
    return matrix


def compress(matrix, removed):
    keep = [i for i in range(19) if i not in removed]
    return acb_mat([[matrix[i,j] for j in keep] for i in keep])


def count_above(matrix, level):
    """在有理分隔值处用厄米消元和惯性定律数出严格较大的特征值。"""
    data = matrix.tolist()
    for i in range(len(data)):
        data[i][i] -= level
    count = 0
    for k in range(len(data)):
        pivot = data[k][k]
        assert pivot.imag.contains(0)
        assert pivot.real > 0 or pivot.real < 0
        count += int(pivot.real > 0)
        for i in range(k+1,len(data)):
            for j in range(k+1,i+1):
                data[i][j] -= data[i][k]*data[j][k].conjugate()/pivot.real
    return count


def uniform_ordering():
    odd, even = (cut_matrix(q) for q in (1,2))
    # 钟边用 2|xy| <= 3|x|²+|y|²/3 作统一的算子上界。
    majorant = acb_mat(even)
    clock_edges = [(a,b) for a,b,d,_ in SHELL["edges"] if d in (724,804)]
    assert len(clock_edges) == 12
    for a,b in clock_edges:
        majorant[a,a] += 3
        majorant[b,b] += arb(1)/3
    even_core = compress(even,{6,7,8})
    odd_core = compress(odd,{6,7,8})
    odd_inner = compress(odd,set(range(6)))
    guard = arb(1)/30
    # 递减排序，每组的最大偶荷特征值低于奇荷的第 j 项。
    for j,cut in enumerate(("21.5","19.045","18","16.7","14","13"),1):
        k = 3*j-2
        if j == 1:
            upper, lower, shift = majorant, odd_inner, 0
        else:
            upper, lower, shift = even_core, odd_core, 3
        n_upper = count_above(upper,arb(cut)-guard)+shift
        n_lower = count_above(lower,arb(cut)+guard)
        assert n_upper < k and n_lower >= j
        print("偶荷序号、奇荷序号、分隔点与两侧惯性计数",k,j,cut,n_upper,n_lower)


def verify():
    ctx.prec = 192
    uniform_ordering()
    center = run_path(str(ROOT / "ch39-native-center-parity-crossing.py"))
    levels, _ = center["setup"](11431)
    even, odd = levels[2], levels[1]
    assert len(even) == len(odd) == 18
    for row,number in [(even,3),(odd,1)]:
        assert sum(int(x > arb("19.26")) for x in row) == number
        assert sum(int(x < arb("19.24")) for x in row) == 18-number
    # 第 11431 枚上一条从 19.24 至 19.26 的饱和线性坡精确达到 3。
    u, _, _ = center["verify"]()
    baseline = sum(even)/sum(odd)
    ratio = u*u*baseline
    assert arb("8.0580925") < ratio < arb("8.0580927")
    assert arb("8.04419052") < 8*baseline < arb("8.04419054")
    floor = arb("2.42")**2*18/19
    assert floor > arb("5.54816")
    print("未分裂、拟合交点及刚度比八的十八层平方质量和比",
          baseline.str(15),ratio.str(13),(8*baseline).str(13))
    print("此前允许相位族内任一拟合交点的质量比下界",floor.str(14))


if __name__ == "__main__":
    verify()
