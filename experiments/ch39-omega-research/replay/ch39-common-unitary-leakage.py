"""同一站点酉矩阵跨实际零点保持轻空间时的味混合限制与最小跃迁概率。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


mix = run_path(str(Path(__file__).with_name("ch39-two-boundary-site-mixing.py")))
adjoint = mix["base"]["adjoint"]


def diagonal(values):
    return acb_mat([[v if i == j else 0 for j in range(3)] for i,v in enumerate(values)])


def square_norm(matrix):
    return sum(abs(z)**2 for z in matrix.entries())


def leakage(first, other, target):
    overlap = adjoint(first)*other
    assert all(overlap[i,j].contains(0) for i in range(3) for j in range(3) if i != j)
    c = [abs(overlap[i,i]) for i in range(3)]
    s = [(1-v**2).sqrt() for v in c]
    assert all(0 < v < 1 for v in s)
    assert all(not s[i].overlaps(s[j]) for i in range(3) for j in range(i))
    d = [a*b for a,b in zip(c,s)]
    cosine, sine, cs = map(diagonal, (c,s,d))

    # 每个出射行的长度给重定相不变的算子范数下界。
    lower = []
    for rows in permutations(range(3)):
        for cols in permutations(range(3)):
            lengths = [(sum(abs(target[rows[i],cols[j]])**2*d[j]**2 for j in range(3))).sqrt()
                       for i in range(3)]
            lower.append(max(abs(lengths[i]-d[i]).lower() for i in range(3)))
    print("任意行列标号的最坏态跃迁振幅下界", min(lower).str(10))
    assert min(lower) > arb("1e-6")

    # 最优化包含任意十六维重酉。取三维压缩的极因子作可行见证。
    b = cs*target*cs
    squared = adjoint(b)*b
    eig = [z.real for z in squared.eig(multiple=True)]
    assert all(v > 0 for v in eig)
    assert all(not v.overlaps(w) for i,v in enumerate(eig) for w in eig[:i])
    unit = diagonal([1,1,1])
    polar = acb_mat(3,3)
    for i,v in enumerate(eig):
        projector = unit
        for j,w in enumerate(eig):
            if j != i:
                projector = projector*(squared-w*unit)/(v-w)
        polar += b*projector/v.sqrt()
    assert max(abs(z) for z in (adjoint(polar)*polar-unit).entries()) < arb("1e-30")
    least = (3-square_norm(cosine*target*cosine)
             - 2*sum(v.sqrt() for v in eig) - sum(v**4 for v in s)) / 3
    example = square_norm(cosine*polar*sine-sine*target*cosine) / 3
    assert 0 < least < example
    assert example-least < arb("1e-12")*example + arb("1e-12")
    print("另一背景均匀三轻态泄漏概率的最优值夹逼",least.str(13),example.str(13))


if __name__ == "__main__":
    ctx.prec = 448
    first = mix["light_modes"](28325)
    _, target = mix["calibration"](mix["boundary_frames"](first))
    for number in (22027, 48067, 58214):
        print("第 28325 枚实现参考 CKM，在第",number,"枚读取跃迁")
        leakage(first, mix["light_modes"](number), target)
