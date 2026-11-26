"""偶数子槽驻点的五个严格零动量方向与双偶偏振的长波方程。"""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-regge-fierz-pauli.py")))
    edges = model["EDGES"]
    metric = sp.Matrix(15,10,lambda a,b: (1+int(model["METRIC"][b][0]!=model["METRIC"][b][1]))
                       *int(all((a+1)&(1<<v) for v in model["METRIC"][b])))
    heavy = sp.Matrix.hstack(*metric.T.nullspace())
    even_difference = {(0,1):1, (2,3):-1}
    even_sum = {(0,1):-2, (2,3):-2, **{(i,j):1 for i in (0,1) for j in (2,3)}}
    cross = {(0,2):1, (0,3):-1, (1,2):-1, (1,3):1}
    modes = metric*sp.Matrix([[h.get(e,0) for h in (even_difference,even_sum,cross)]
                               for e in model["METRIC"]])
    parts = model["base"]["canonical"]()
    w,s,f,g = sp.symbols("w s f g",real=True)
    pencil = sp.zeros(3)
    for cosine, weight in [(0,f), (sp.Rational(1,2),g)]:
        local = sum((m for _,c,m in parts if c==cosine),sp.zeros(10))
        local = (local + local.T)/2
        coefficients = defaultdict(lambda: [0,0,0,0])
        for order in permutations(range(4)):
            vertices = [sum(1<<i for i in order[:j]) for j in range(5)]
            for (i,j), value in local.todok().items():
                a,b = [vertices[x] for x in edges[i]]
                c,d = [vertices[x] for x in edges[j]]
                delta = [((c>>k)&1)+((d>>k)&1)-((a>>k)&1)-((b>>k)&1) for k in range(4)]
                time = sum(delta)
                space = sum(delta[:2])-sum(delta[2:])
                row = coefficients[b-a-1,d-c-1]
                for pos, term in enumerate((1,time**2,2*time*space,space**2)):
                    row[pos] += value*term
        blocks = [sp.zeros(15) for _ in range(4)]
        for key, row in coefficients.items():
            for i, entry in enumerate(row):
                blocks[i][key] = entry
        mixed = metric.T*blocks[0]*heavy
        assert mixed.rank()==5 and metric.T*blocks[0]*metric==sp.zeros(10)
        null = metric * sp.Matrix.hstack(*mixed.T.nullspace())
        assert null[[0,1,3,7,14],:] == sp.zeros(5)
        assert blocks[0]*modes == sp.zeros(15,3)
        pencil -= weight*modes.T*(w*w*blocks[1]+w*s*blocks[2]+s*s*blocks[3])*modes/8
    assert sp.expand(pencil[2,2]+8*(f*s**2+(3*f+2*g)*w**2)/5)==0
    assert pencil[:2,2]==sp.zeros(2,1)
    det = sp.expand(pencil[:2,:2].det().subs(s,1))
    root = -f/(3*f+2*g)
    on_cone = sp.factor(det.subs(w**4,root**2).subs(w**2,root))
    assert sp.factor(on_cone-64*f*(4*f+g)*(6*f-g)/(25*(3*f+2*g)))==0
    assert sp.expand(det.subs(g,f) - 16*f**2*(75*w**4-90*w**2-1)/25)==0
    print("双奇零锥处双偶二维主块的行列式",on_cone)


if __name__=="__main__":
    verify()
