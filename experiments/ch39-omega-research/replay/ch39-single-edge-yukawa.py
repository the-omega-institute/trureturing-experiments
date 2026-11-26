"""三个轻 Higgs 下，同站点项与三种单腿壳边输运的共同复叠加。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


filters = run_path(str(Path(__file__).with_name("ch39-flavor-shell-filters.py")))
shell = filters["local"]["local"]["shell"]


def mode_columns(rows):
    ordered = sorted(rows, key=lambda row: row[1])
    return acb_mat(list(zip(*(vector for _, _, vector in ordered))))


def positive_definite(matrix):
    # 二阶厄米矩阵的 Sylvester 判据。
    return matrix[0, 0].real > 0 and matrix.det().real > 0


def largest_ratio(up, other):
    adjoint = filters["local"]["adjoint"]
    gram_up, gram_other = [adjoint(c) * c for c in (up, other)]
    assert positive_definite(gram_other)
    low, high = arb(0), arb(256)
    assert positive_definite(high**2 * gram_other - gram_up)
    for _ in range(32):
        mid = (low + high) / 2
        diff = mid**2 * gram_other - gram_up
        if positive_definite(diff):
            high = mid
        else:
            assert diff[0, 0].real < 0 or diff.det().real < 0
            low = mid
    return low.union(high)


def channels(fields, conjugate_higgs):
    # 每条腿同时提供原模与由相邻顶点搬来的模。
    rows = [[] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                coefficients = []
                for moving in (None, 0, 1, 2):
                    a, b, h = [transported if n == moving else original for n, (original, transported) in enumerate(fields)]
                    coefficients.append(sum(a[x,i] * b[x,j] * h[x,k] for x in range(19)))
                # k 始终指正荷 Higgs 的模编号，取共轭时轮换荷变号。
                if (i + j + (-k if conjugate_higgs else k)) % 3 == 0:
                    rows[k].append(coefficients)
                else:
                    assert all(c.contains(0) for c in coefficients)
    return [acb_mat(entries) for entries in rows]


def select_columns(matrix, columns):
    return acb_mat([[matrix[i,j] for j in columns] for i in range(matrix.nrows())])


def joint_ratio(up, down, leptons):
    adjoint = filters["local"]["adjoint"]
    gu, gd, ge = [adjoint(c) * c for c in (up, down, leptons)]
    # 向下取整窗口下端，给予候选更多余量。
    denominator = (arb("61.36519")**2 * gd + arb("84.98147")**2 * ge) / 2
    for matrix in (denominator, arb(".02363")**2 * denominator - gu):
        for size in range(1, 5):
            minor = acb_mat([[matrix[i,j] for j in range(size)] for i in range(size)])
            assert minor.det().real > 0
    roots = (denominator.inv() * gu).eig(algorithm="rump", multiple=True)
    assert all(r.imag.contains(0) for r in roots)
    return max(r.real.upper() for r in roots).sqrt()


def verify():
    ctx.prec = 448
    for number in (22027, 28325):
        height = acb.zeta_zero(number).imag
        neutral = shell["laplacian"](0, height)
        center = sum(neutral[i,i].real for i in range(19)) / 18
        fields = {}
        for q in (1, -4, 2, -3, 6, 3):
            lap = shell["laplacian"](q, height)
            adjacency = acb_mat([[-lap[i,j] if i != j else 0 for j in range(19)] for i in range(19)])
            modes = mode_columns(filters["spectra"](q, height, center))
            fields[q] = (modes, adjacency * modes)
        # 轻子双重态和共轭 Higgs 同荷，保留各自的特征标编号和相位约定。
        conjugate_higgs = tuple(block.conjugate() for block in fields[3])
        up = channels((fields[1], fields[-4], fields[3]), False)
        down = channels((fields[1], fields[2], conjugate_higgs), True)
        leptons = channels((fields[-3], fields[6], conjugate_higgs), True)
        for column, name in enumerate(("左费米子", "右费米子", "Higgs"), start=1):
            maxima = []
            for sector in (down, leptons):
                ratios = [largest_ratio(select_columns(u, (0, column)), select_columns(d, (0, column)))
                          for u, d in zip(up, sector)]
                largest = max(ratios, key=lambda r: r.upper())
                assert all(r is largest or r.upper() < largest.lower() for r in ratios)
                maxima.append(largest)
            assert maxima[0] < arb("1.229") and maxima[1] < arb("1.152")
            print(number, name, "同站点加单边输运，上/下与上/带电轻子范数比上确界",
                  *(value.str(10) for value in maxima))

        joint = [joint_ratio(u, d, e) for u, d, e in zip(up, down, leptons)]
        print(number, "同站点与三条单边项相干叠加时的联合参考比上限",
              *(value.str(10) for value in joint))


if __name__ == "__main__":
    verify()
