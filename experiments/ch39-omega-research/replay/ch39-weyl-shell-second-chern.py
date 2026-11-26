"""把原壳的简单 Weyl 点接入线性复跃迁链，核对四维球面上的有向陈数。"""
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    root = Path(__file__).resolve().parent
    run_path(str(root / "ch39-first-chern-edge.py"))["verify"](root)
    u, v = sp.symbols("u v", real=True)
    z = u + sp.I*v
    identity = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.diag(1, -1)
    gamma = [sp.kronecker_product(sz, s) for s in (sx, sy, sz)]
    gamma += [sp.kronecker_product(s, identity) for s in (sx, sy)]
    product = sp.eye(4)
    for i, left in enumerate(gamma):
        product *= left
        for j, right in enumerate(gamma):
            assert left*right + right*left == 2*int(i == j)*sp.eye(4)
    assert product == -sp.eye(4)
    for k in (1, 2, 3):
        path = sp.zeros(2*k)
        for a in range(2*k-1):
            forward = z.conjugate() if a % 2 == 0 else 1
            path[a, a+1] = forward
            path[a+1, a] = sp.conjugate(forward)
        ends = [0, 2*k-1]
        result = path.extract(ends, ends)
        if k > 1:
            inner = list(range(1, 2*k-1))
            dimers = path.extract(inner, inner)
            assert dimers.det() == (-1)**(k-1)
            links = path.extract(ends, inner)
            result -= links*dimers.inv()*links.adjoint()
        target = (-1)**(k-1) * z.conjugate()**k
        assert (result-sp.Matrix([[0, target], [sp.conjugate(target), 0]])).applyfunc(sp.expand) == sp.zeros(2)
        print("复跃迁次数、十九点层的总数、外向四维球面第二陈数", k, 2*k, k)


if __name__ == "__main__":
    verify()
