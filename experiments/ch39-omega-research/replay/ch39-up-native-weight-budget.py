"""给出右外壳可入轻上型窗的原生层振幅及谱转移概率的必要条件。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


RAYS = (1, 6, 7, 10, 12, 16, 18)
SELECTED = (1, 16, 18)
BUDGET = (".02873", ".20725", ".23702")


def trivial_rays(target):
    matrices = {}
    for j in range(19):
        ray = acb_mat([[target[14, j]],
                       [sum(target[i, j] for i in (8, 7, 6)) / arb(3).sqrt()]])
        matrix = ray * ray.conjugate().transpose()
        if matrix.trace().real > 0:
            matrices[j + 1] = matrix
        else:
            assert matrix.contains(acb_mat(2, 2))
    assert tuple(matrices) == RAYS
    assert sum(matrices.values(), acb_mat(2, 2)).contains(acb_mat([[1, 0], [0, 1]]))
    return matrices


def verify(root):
    character = run_path(str(root / "ch39-up-unrestricted-character-floor.py"))
    weights, target = character["shell_data"](root)
    matrices = trivial_rays(target)
    source = sorted(w.lower() for w in weights)
    limit = arb(".00000629")
    identity = acb_mat([[1, 0], [0, 1]])
    guaranteed = {j: arb(1) for j in SELECTED}
    eligible = []
    for smallest in combinations(RAYS, 3):
        for fourth in (j for j in RAYS if j not in smallest):
            # 用每个等级的严格下端构造向下包围的二阶矩阵。
            coefficients = [source[0] if j in smallest else source[3] if j == fourth else source[4]
                            for j in RAYS]
            matrix = sum((a * w for a, w in zip(matrices.values(), coefficients, strict=True)),
                         acb_mat(2, 2))
            gap = ((matrix[0, 0].real - matrix[1, 1].real)**2 + 4*abs(matrix[0, 1])**2).sqrt()
            low = (matrix.trace().real - gap) / 2
            if low > limit:
                continue
            assert low < limit < low + gap
            eligible.append((smallest, fourth))
            cap = (limit - low) / gap
            low_projector = ((low + gap) * identity - matrix) * (1 / gap)
            for j in SELECTED:
                near = (low_projector * matrices[j]).trace().real
                far = matrices[j].trace().real - near
                amplitude = ((1-cap)*near).sqrt() - (cap*far).sqrt()
                assert amplitude > 0
                guaranteed[j] = min(guaranteed[j], (amplitude**2).lower())

    assert eligible == [(SELECTED, j) for j in RAYS if j not in SELECTED]
    print("保留的四种重排（从一开始的目的层）", eligible)
    for j, coefficient in zip(SELECTED, BUDGET, strict=True):
        assert guaranteed[j] > arb(coefficient)
        print("目的层、读出方向投影概率下界、目的层振幅上界",
              j, guaranteed[j].str(14), (limit/guaranteed[j]).str(14))
    # 两条目的线在三个最低源层外的谱转移，以及在六个最高源层内的谱转移。
    pair_weight = arb(BUDGET[1])
    assert source[3] > arb(".00025507312") and source[-6] > arb(".947554")
    complement = limit / (pair_weight * source[3])
    highest = limit / (pair_weight * source[-6])
    assert 2 - complement > arb("1.8810") and highest < arb(".00003204")
    print("第十六、十八层回到最低三源层的概率和下界、回到最高六源层的概率和上界",
          (2-complement).str(14), highest.str(14))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
