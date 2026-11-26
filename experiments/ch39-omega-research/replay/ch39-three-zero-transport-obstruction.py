"""三枚规范候选的共同八系数输运顶角与任意逐枚轻 Higgs 方向。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


study = run_path(str(Path(__file__).with_name("ch39-three-leg-yukawa.py")))
positive = study["positive"]


def deficits(number, fraction, weights):
    blocks, _, _, _ = study["tensors"](number)
    quark = (fraction*arb("61.36519"))**2
    lepton = (fraction*arb("84.98147"))**2
    a,b = weights
    equal, weighted = [], []
    for block in blocks:
        u,d,e = [study["adjoint"](c)*c for c in block]
        delta_d, delta_e = quark*d-u, lepton*e-u
        equal.append((delta_d+delta_e)/2)
        weighted.append(a*delta_d+b*delta_e)
    return equal, weighted


def verify():
    ctx.prec = 448
    fraction = arb(".987")
    cap = arb(".001")
    cases = ((11661,(2,40)), (11934,(0,21)), (22027,(0,37)))
    forms = [deficits(n,fraction,w) for n,w in cases]
    for indices in product(range(3), repeat=3):
        if indices != (2,2,2):
            positive(sum((pair[0][k] for pair,k in zip(forms,indices)), acb_mat(8,8)))
    for equal,_ in forms:
        for k in (0,1):
            positive(equal[2]+cap*equal[k])
    for vertex in product((arb(0),cap),repeat=6):
        total = acb_mat(8,8)
        for i,(_,weighted) in enumerate(forms):
            total += weighted[2]+vertex[2*i]*weighted[0]+vertex[2*i+1]*weighted[1]
        positive(total)
    print("正定证书通过：26 个特征标组合，6 个纯度约束，64 个 Higgs 权重盒顶点")
    print("六个部门比至少一个低于下限的",fraction.str(5))


if __name__ == "__main__":
    verify()
