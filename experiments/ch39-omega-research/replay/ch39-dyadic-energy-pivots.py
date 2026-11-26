"""用整数矩阵和有理逆核验冻结倍频壳能量的舒尔枢轴。"""
from fractions import Fraction

from flint import fmpq_mat, fmpz_mat


def inverse_row(n):
    """计算从各壳和恢复第 n 个输入所需的系数，省去恒为零的第零壳。"""
    row = {}
    while True:
        row[n] = row.get(n, 0) + 1
        if n > 1:
            row[n - 1] = row.get(n - 1, 0) - 1
        if n % 2:
            return row
        n //= 2


def squared_length(n):
    odd = n
    order = 0
    while odd % 2 == 0:
        odd //= 2
        order += 1
    if odd > 1:
        return 2 * (order + 1)
    return max(1, 2 * order - 1)


def verify():
    for cutoff in range(2, 35):
        indices = range(1, cutoff)
        inc = fmpz_mat([[int(n <= j < 2*n) for n in indices] for j in indices])
        rows = [inverse_row(n) for n in indices]
        inverse = fmpz_mat([[row.get(j, 0) for j in indices] for row in rows])
        identity = fmpz_mat([[int(i == j) for j in indices] for i in indices])
        assert inverse * inc == identity
        band = fmpz_mat([[max(0, min(cutoff, 2*min(m, n)) - max(m, n))
                          for n in indices] for m in indices])
        assert band == inc.transpose() * inc and band.det() == 1
        band_inverse = fmpq_mat(band).inv()
        assert all(band_inverse[n-1, n-1] == squared_length(n) for n in indices)

        # 按降序消元时，每个尾主子式都是一。
        for first in indices:
            tail = range(first-1, cutoff-1)
            assert fmpz_mat([[band[i, j] for j in tail] for i in tail]).det() == 1
    assert squared_length(6) == 4 and squared_length(8) == 5
    # 取一个有限算术截止，求所有坐标与消元次序的锐下界。
    cutoff = 12093
    norms = [squared_length(n) for n in range(1, cutoff)]
    maximum = max(norms)
    minimizers = [i for i, norm in enumerate(norms, 1) if norm == maximum]
    assert maximum == 25 and minimizers == [8192]
    print("算术截止、最小舒尔枢轴、取等坐标", cutoff, Fraction(1, maximum), minimizers)


if __name__ == "__main__":
    verify()
