"""固定右内壳上型源谱时，至多八站的局部酉不能达到轻夸克参考窗。"""
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb, ctx


def rows(matrix, indices):
    return acb_mat([[matrix[i,j] for j in range(matrix.ncols())] for i in indices])


def numeric(matrix):
    return np.array([[complex(x) for x in line] for line in matrix.tolist()])


def multipliers(block, rhs, norm2):
    # 浮点运算只选择每个目的层的一个实数乘子，符号和最终下界均用区间核验。
    eigenvalues, basis = np.linalg.eigh(numeric(block))
    weights = np.abs(basis.conj().T@numeric(rhs))**2
    probability = np.array(list(map(float,norm2)))
    left = np.full(19, -eigenvalues[0])
    right = np.sqrt(weights.sum(axis=0)/np.maximum(probability,1e-28))+1e-10
    for _ in range(30):
        middle = (left+right)/2
        denominator = np.maximum(eigenvalues[:,None]+middle,1e-25)
        above = (weights/denominator**2).sum(axis=0)>probability
        left = np.where(above,middle,left)
        right = np.where(above,right,middle)
    return [arb(format(float(x),".16g")) for x in right]


def certify_set(support, h, target, four, positive, threshold, fixed_adjoint=None):
    inside = rows(target, support)
    outside = acb_mat(target)
    for j in range(19):
        for i in support:
            outside[i,j] = 0
    if fixed_adjoint is not None:
        fixed = [i for i in range(19) if i not in support]
        old = rows(outside, fixed)
        transformed = fixed_adjoint * old
        for k,i in enumerate(fixed):
            for j in range(19):
                outside[i,j] = transformed[k,j]
    block = acb_mat([[h[i,j] for j in support] for i in support])
    mixed = h*outside
    rhs = rows(mixed, support)
    norm2 = [sum(abs(inside[i,j])**2 for i in range(len(support))) for j in range(19)]
    gamma = multipliers(block,rhs,norm2)
    floors = []
    for j,g in enumerate(gamma):
        mat = acb_mat(block)
        for i in range(len(support)):
            mat[i,i] += g
        if not g >= 0:
            try:
                positive(mat)
            except AssertionError:
                g = arb(0)
                mat = block
        col = acb_mat([[rhs[i,j]] for i in range(len(support))])
        correction = (col.conjugate().transpose()*mat.solve(col))[0,0].real
        constant = sum(outside[i,j].conjugate()*mixed[i,j] for i in range(19)).real
        floors.append(max(arb(0),(constant-g*norm2[j]-correction).lower()))
    measured = acb_mat([[four[i,j]*floors[j] for j in range(19)] for i in range(4)])
    bound = measured*four.conjugate().transpose()
    for i in range(4):
        bound[i,i] -= threshold
    positive(bound)


def verify(root, size=8, filename="ch39-up-eight-site-cover.tsv", threshold=".0000073", exceptions=()):
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spec = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    positive = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["positive"]
    weighted = run_path(str(root / "ch39-common-unitary-quark-right-trace.py"))["weighted"]
    n = 2_306_888_100_086
    _, labels, _ = training(root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    _, source = spec(shell,-4,n)
    _, target = spec(shell,-4,n+1)
    assert min(labels[n,-4]) > 0
    h = weighted(source, labels[n,-4])
    four = rows(target,(14,8,7,6))
    covered = set()
    entries = (root/filename).read_text().splitlines()
    for line in entries:
        support = tuple(map(int,line.split()))
        assert tuple(sorted(set(support))) == support and size <= len(support) <= 19
        assert 0 <= min(support) and max(support) < 19
        certify_set(support, h, target, four, positive, arb(threshold))
        covered.update(combinations(support,size))
    assert not covered.intersection(exceptions)
    assert len(covered)+len(set(exceptions))==comb(19,size)
    print("活动站数、已复核超集、覆盖支撑数、最轻质量下界",size,len(entries),len(covered),threshold)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
