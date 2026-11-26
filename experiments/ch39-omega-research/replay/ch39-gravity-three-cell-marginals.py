"""三胞中心化标记的正概率证书与真实零点的三位联合读数。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
import numpy as np
import sympy as sp


ROOT = Path(__file__).parent


def three_cell_coefficients():
    model = run_path(str(ROOT / "ch39-hinge-curvature-inverse.py"))
    matrix, sines, _ = model["readouts"]()
    scaled = [6 * op for op in model["operator_coefficients"](matrix)]
    assert all(v.is_Integer for op in scaled for v in op)
    operators = np.array([op.tolist() for op in scaled], dtype=np.int64)
    assert np.abs(operators).max() == 12
    # 每个收缩至多含 216 个绝对值不超过 12³ 的乘积。
    tensor = np.einsum('iab,jbc,kca->ijk', operators, operators, operators)
    assert np.array_equal(tensor, tensor.transpose(1,0,2))
    assert np.array_equal(tensor, tensor.transpose(1,2,0))
    assert sum(map(abs, map(int, tensor.flat))) == 216 * sp.Rational(15400, 3)
    sine_sq = np.array([int(4*s**2) for s in sines], dtype=np.int64)
    score = tensor**2 * np.einsum('i,j,k->ijk', sine_sq, sine_sq, sine_sq)
    assert score.max() == (216 * 8)**2 * 64
    i = next(i for i, s in enumerate(sines) if s == 1 and abs(tensor[i,i,i]) == 8*216)
    hinges = sorted((a,b) for a,b in combinations(range(1,16),2) if a & b == 0)
    coeff = int(tensor[i,i,i]) // 216
    print("选定同一四价反向面对、三胞有序迹系数", hinges[i], coeff)
    assert hinges[i] == (1, 6) and coeff == 8
    return coeff


def check_marginals():
    # 三胞最短周期，三个块交错。每胞有三个相异槽，验证两个完整胞的全部配置。
    m = sp.Rational(2,5)
    alpha = sp.Rational(1,3)
    joint = {}
    reference = {}
    for bits in product((0,1), repeat=9):
        law = m**sum(bits) * (1-m)**(9-sum(bits))
        weight = sp.S.One
        for z in range(3):
            centered = [bits[3*((z+t)%3)+t] - m for t in range(3)]
            weight *= 1 + alpha * sp.prod(centered)
        assert weight > 0
        joint[bits], reference[bits] = law * weight, law
    assert sum(joint.values()) == 1
    for a,b in combinations(range(3),2):
        residual = {}
        for bits,w in joint.items():
            key = bits[3*a:3*a+3] + bits[3*b:3*b+3]
            residual[key] = residual.get(key, 0) + w - reference[bits]
        assert len(residual) == 64 and all(v == 0 for v in residual.values())
    triple = sum(w - reference[bits] for bits,w in joint.items() if bits[0]*bits[4]*bits[8])
    assert triple == alpha * (m*(1-m))**3


def verify():
    coeff = three_cell_coefficients()
    check_marginals()
    shell = run_path(str(ROOT / "ch39-arc-relative-determinant.py"))
    ctx.prec = 160
    lower, count = 10855, 1238
    zeros = acb.zeta_zeros(lower, count)
    assert len(zeros) == count and all(u.imag < v.imag for u,v in zip(zeros, zeros[1:]))
    m = 4 - 2 * arb(3).sqrt()
    eps = arb(3)/72155
    # 选择每个反向面对的一面，每条腿的平均各给 1/2。
    coefficient = eps * coeff * (m*(1-m))**3 / 8
    readings = [m**3 + coefficient * shell["phases"](zero.imag)[724].real for zero in zeros]
    delta = readings[11366-lower] - readings[11431-lower]
    order = sorted(range(count), key=lambda i: float(readings[i]))
    gaps = [(readings[j] - readings[i], lower+i, lower+j) for i,j in zip(order, order[1:])]
    assert all(gap > 0 for gap,_,_ in gaps)
    tight, first, second = min(gaps, key=lambda row: float(row[0]))
    assert all(gap > tight for gap,i,j in gaps if (i,j) != (first,second))
    assert arb("4.69018860646e-7") < delta < arb("4.69018860648e-7")
    assert (first, second) == (11870, 11943)
    assert arb("4.47664151643e-13") < tight < arb("4.47664151645e-13")
    print("11366 与 11431 的三位全一概率差", delta.str(16))
    print("全部 1238 个读数的最小间隔及编号", tight.str(16), first, second)


if __name__ == "__main__":
    verify()
