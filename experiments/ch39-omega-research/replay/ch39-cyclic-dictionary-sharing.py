"""105圈权的三个等分相位字典锐门槛，以及原壳两面字的正质量表示。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


N = 105


def orbit(size, q, difference):
    paired = -q-int(difference)
    return min(q % size, paired % size)


def cyclic_thresholds():
    # 实质量的傅里叶矩按 q ↔ -q 共轭。边 q→q+1 的反射索引为 -q-1。
    for difference, minimum in ((False, 210), (True, 211)):
        for size in range(1, minimum):
            assert len({orbit(size, q, difference) for q in range(1, N+1)}) < N
        assert len({orbit(minimum, q, difference) for q in range(1, N+1)}) == N
    assert (-105) % 211 == 106
    # 每个目标振幅作为独立形式变量，以下用基向量精确检验傅里叶赋值。
    for size, diff, symmetric in ((210, False, True), (211, True, False), (212, True, True)):
        for col in range(1, N+1):
            # 记录未缩放非零频矩的两倍，全部实虚分量都在整数内。
            moments = [0j]*size
            if not diff:
                moments[col] = 2
            else:
                if size == 211:
                    moments[N] = -1j*int(col == N)
                    last = N-1
                else:
                    last = N
                for q in range(last, 0, -1):
                    moments[q] = moments[q+1]-2*int(q == col)
            for q in range(1, N+1):
                moments[-q] = complex(moments[q]).conjugate()
            assert all(moments[q] == moments[-q].conjugate() for q in range(size))
            if symmetric:
                assert all(v.imag == 0 for v in moments)
            for q in range(1, N+1):
                value = (moments[q+1]-moments[q]) if diff else moments[q]
                desired = (2j if diff and size == 211 and q == N else 2)*int(q == col)
                assert value == desired
    # 各非零傅里叶矩的模不超过 ε√105。反演后的质量全都为正。
    assert 211*arb(N).sqrt()/4096 < arb(".53")
    print("等分单位根字典的最小点数", 210, 211, 212)


def shell_dictionary(root):
    geometry = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 224
    faces = geometry["enumerate_faces"]()
    first, second = (7514, -3), (29189, -41)
    native = {z for k, z in faces.values() if k == (1, 0)}
    assert first in native and second in native
    label = geometry["product"](first, geometry["conjugate"](second))
    assert label == (108469, 109)
    unit = geometry["unit"](label)
    # 两个相反取向的同类面净钟频率为零。余弦块使用一对 ±j 原子的总质量。
    n = N+1
    angle = unit.arg()
    assert 0 < angle and 2*n*angle < 2*arb.pi()+arb(".1")
    assert all(not (unit**j-1).contains(0) for j in range(1, 2*n+1))
    matrix = arb_mat([[(q*j*angle).cos() for j in range(n+1)] for q in range(n+1)])
    inverse = matrix.inv()
    h = min(inverse[j, 0].lower() for j in range(n+1))
    norm = max(sum(abs(inverse[j, q]) for q in range(n+1)).upper() for j in range(n+1))
    assert h > arb(".00470733") and norm < arb("1.289122")
    positive = h-norm*arb(N).sqrt()/4096
    assert positive > arb(".00148234")
    print("原壳213点的零矩基准最小组质量、逆矩阵行和、扰动后最小组质量下界",
          *(v.str(16) for v in (h, norm, positive)))


if __name__ == "__main__":
    cyclic_thresholds()
    shell_dictionary(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
