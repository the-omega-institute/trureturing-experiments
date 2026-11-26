"""在所有硬截止低于 216 的情形达成交替外点读数的逐项最小改动。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

previous = run_path(str(Path(__file__).with_name("ch39-three-clock-alternating-outer.py")))


def verify():
    points, _, outside, inside, signs, coupling = previous["verify"]()
    _, _, g0, delta, _, _ = previous["plane"]["shell_data"]()
    submatrix = previous["plane"]["pairs"]["submatrix"]
    size = len(points)
    a, b = fmpq_mat(size, 1), fmpq_mat(size, 1)
    for k in range(6):
        a[outside[k], 0] = b[inside[k], 0] = signs[k, 0]
    alpha = (delta*a)[outside[0], 0]
    assert delta*a == alpha*a + coupling*b
    endpoint = delta - (coupling/6)*(a*b.transpose() + b*a.transpose())
    for i in range(size):
        endpoint[i, i] += -alpha if i in outside else 700
    assert endpoint*a == fmpq_mat(size, 1)
    # 删除核向量中一个非零坐标，剩余主块正定便可推出整个端点正半定。
    keep = [i for i in range(size) if i != outside[0]]
    for k in range(1, len(keep) + 1):
        assert submatrix(endpoint, keep[:k], keep[:k]).det() > 0
    assert endpoint.rank() == 18
    differences = [abs(delta[i, j]-endpoint[i, j]) for i, j in combinations(range(size), 2)]
    assert max(differences) == coupling/6 and differences.count(coupling/6) == 36
    bound = fmpq(944, 10**9)
    inv_norm = max(sum(abs(x) for x in row) for row in g0.inv().tolist())
    for matrix in (delta, endpoint):
        assert max(sum(abs(x) for x in row) for row in matrix.tolist())/10**10 < bound
    assert inv_norm*bound < (1-fmpq(98, 1000))**2
    print("All cutoffs below 216, optimal entrywise bound",
          float(coupling/(6*10**10)), "attained with at most 37 auxiliary states")


if __name__ == "__main__":
    verify()
