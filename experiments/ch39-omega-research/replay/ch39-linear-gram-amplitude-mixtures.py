"""任意多个复系数线性 Gram 振幅的模平方之和无法恢复两个锁定零点的六质量。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    inputs = run_path(str(Path(__file__).with_name("ch39-two-zero-source-correlations.py")))
    (au,ad),_,factor,adj = inputs["certified_inputs"]()
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    gu,gd = adj(au)*au, adj(ad)*ad
    # 按下型迹内积消去近共线分量，三角换基保留全部线性振幅。
    basis = []
    for word in (identity,gu,gd):
        element = word
        for previous in basis:
            q = (adj(ad*previous)*(ad*element)).trace().mid()
            n = (adj(ad*previous)*(ad*previous)).trace().real.mid()
            element -= previous*(q/n)
        norm = sum(abs(z)**2 for z in (ad*element).entries()).sqrt().mid()
        assert norm > 0
        basis.append(element/norm)
    pencil = arb("3.53")*gd - factor*gu
    moment = acb_mat([[(adj(left)*pencil*right).trace() for right in basis] for left in basis])
    for n in (1,2,3):
        minor = acb_mat([[moment[i,j] for j in range(n)] for i in range(n)]).det().real
        assert minor > 0
        print("线性振幅矩阵的正定主子式",n,minor.lower().str(12))


if __name__ == "__main__":
    verify()
