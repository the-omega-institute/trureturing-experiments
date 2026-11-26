"""在截止 200 上求交替读数与截止外单项改动幅度的精确关系。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq_mat, fmpq

previous = run_path(str(Path(__file__).with_name("ch39-three-clock-alternating-outer.py")))


def verify():
    points, distance, outside, inside, signs, c = previous["verify"]()
    identity = fmpq_mat([[int(i == j) for j in range(6)] for i in range(6)])
    sign_projection = signs*signs.transpose()/6
    complement = identity - sign_projection
    assert complement*complement == complement and complement.rank() == 5
    insertion = fmpq_mat(len(points), 6)
    for j, (i, k) in enumerate(zip(outside, inside)):
        insertion[i, j] = insertion[k, j] = 1
    add_full = c*insertion*insertion.transpose()
    add_null = c*insertion*complement*insertion.transpose()
    short = {tuple(sorted(edge)) for edge in zip(outside, inside)}
    assert short == {(i, j) for i, j in combinations(range(19), 2)
                     if 200 < distance[i][j] <= 216}
    # 216 构造在六条新边上缺 c，在其余所有截止外边上已精确匹配。
    extreme = []
    for i, j in combinations(range(19), 2):
        if distance[i][j] <= 200:
            continue
        baseline_error = -c if (i, j) in short else fmpq(0)
        assert baseline_error + add_full[i, j] == 0
        extreme.append(abs(baseline_error + add_null[i, j]))
    assert max(extreme) == c/6
    for correction, rank, coefficient in [(add_full, 6, c), (add_null, 5, 0)]:
        assert correction.rank() == rank
        cross = previous["plane"]["pairs"]["submatrix"](correction, inside, outside)
        assert cross*signs == coefficient*signs
    _, _, g0, _, _, _ = previous["plane"]["shell_data"]()
    inv_norm = max(sum(abs(x) for x in row) for row in g0.inv().tolist())
    max_trace = fmpq(2, 10**6) + sum(add_full[i, i] for i in range(19))/10**10
    assert inv_norm*max_trace < (1-fmpq(98, 1000))**2
    print("Optimal off-cutoff entry bound for a null outer probe:",
          float(c/(6*10**10)), "at cutoff 200")


if __name__ == "__main__":
    verify()
