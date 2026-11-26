"""有限 Born 铰链概率、相位求和与仿射度量的一阶曲率。"""
from collections import defaultdict
from itertools import combinations_with_replacement, permutations
from pathlib import Path
from runpy import run_path

import sympy as sp


def affine_curvature(flat):
    edges = flat["EDGES"]
    tensor_entries = list(combinations_with_replacement(range(4),2))
    squared = sp.Matrix([sp.Rational((j-i)*(5-j+i),5) for i,j in edges])
    # 每个面类型内共同的 sqrt(5)/面积 因子可在核对零时约去。
    rows = [(face,list(squared.T*matrix)) for face,_,matrix in flat["canonical"]()]
    moments = defaultdict(int)
    for order in permutations(range(4)):
        vertices = [sum(1<<v for v in order[:k]) for k in range(5)]
        for (first,middle,last),coeffs in rows:
            origin = vertices[first]
            face = (vertices[middle]-origin,vertices[last]-vertices[middle])
            for (i,j),coef in zip(edges,coeffs):
                if not coef:
                    continue
                a,b = vertices[i],vertices[j]
                step = b-a
                double_midpoint = [((a>>v)&1)+((b>>v)&1)-2*((origin>>v)&1) for v in range(4)]
                for l,r in tensor_entries:
                    if (step>>l)&(step>>r)&1:
                        for axis,factor in enumerate([2,*double_midpoint]):
                            moments[face,l,r,axis] += coef*factor*(1+int(l!=r))
    assert len({idx[0] for idx in moments}) == 50 and not any(moments.values())


def verify():
    root = Path(__file__).parent
    flat = run_path(str(root/"ch39-regge-flat-cell.py"))
    affine_curvature(flat)
    x = sp.Symbol("x")
    assert sp.Poly(x**4-12,x).is_irreducible
    beta_squared = sp.sqrt(3)/2
    four = 1/(1+beta_squared)
    six = beta_squared/(1+beta_squared)
    assert sp.simplify(four+six-1) == 0
    assert sp.simplify(30*four-20*sp.sqrt(3)*six) == 0
    assert sp.simplify(four*sp.sqrt(5)-six/sp.sqrt(sp.Rational(3,20))) == 0
    u,p,dp,b,c = sp.symbols("u p dp b c",real=True)
    phase = b*u+c*u**2/2
    # 只需求看平直背景处的二阶喷射。
    exp_jet = 1+sp.I*phase-phase**2/2
    dz = (p+dp*u)*(exp_jet-1)
    action = -sp.I*(dz-dz**2/2)
    assert sp.expand(action).coeff(u,1) == p*b
    assert sp.simplify(2*sp.expand(action).coeff(u,2)-p*c-2*dp*b-sp.I*p*(1-p)*b**2) == 0
    print("两类非零负载概率、实幅度比的最低数域次数、消去仿射度量的一阶角亏面数",
          sp.radsimp(four),sp.radsimp(six),4,50)


if __name__ == "__main__":
    verify()
