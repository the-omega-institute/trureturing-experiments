"""互补 Pell 支的交叉偏振、制备刚度与交替放缩的全槽位极限。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import arb, ctx


ROOT = Path(__file__).parent
flat = run_path(str(ROOT/"ch39-regge-flat-cell.py"))
stability = run_path(str(ROOT/"ch39-clock-hinge-stability.py"))
ground = run_path(str(ROOT/"ch39-clock-hinge-ground.py"))
CROSS = {5: 1, 9: -1, 6: -1, 10: 1}


def reference_constants():
    edges = flat["EDGES"]
    lengths = [sp.Rational((j-i)*(5-j+i),5) for i, j in edges]
    phase = -sum((part for _, _, part in flat["canonical"]()), sp.zeros(10))
    prep = sp.zeros(10)
    for i, j in edges:
        face = list(combinations([k for k in range(5) if k not in (i, j)], 2))
        total = sum(lengths[edges.index(e)] for e in face)
        row = sp.zeros(1,10)
        for edge in face:
            col = edges.index(edge)
            row[col] = (total/2 - lengths[col])/4
        prep += row.T*row / (6 if j-i in (1,4) else 4)
    assert prep.det() == sp.Rational(1,20873541058560000000000)
    alternate = sp.Matrix([q*((-1)**i+(-1)**j) for q, (i,j) in zip(lengths,edges)])
    assert 24*(alternate.T*phase*alternate)[0] == -sp.Rational(1824,125)
    assert 24*(alternate.T*prep*alternate)[0] == sp.Rational(732,125)
    cost = 0
    for order in permutations(range(4)):
        masks = [sum(1<<i for i in order[:j]) for j in range(5)]
        vector = sp.Matrix([CROSS.get(masks[j]-masks[i],0) for i,j in edges])
        cost += (vector.T*prep*vector)[0]/4
    assert cost == sp.Rational(17,200)


def static_symbol(local, space):
    value = arb(0)
    momentum = (1, 1, -1, -1)
    for order in permutations(range(4)):
        vertices = [sum(1<<i for i in order[:j]) for j in range(5)]
        selected = []
        for row, (i,j) in enumerate(flat["EDGES"]):
            weight = CROSS.get(vertices[j]-vertices[i],0)
            if weight:
                offset = sum(sign for d,sign in enumerate(momentum) if vertices[i] & (1<<d))
                selected.append((row,weight,offset))
        for i, u, a in selected:
            for j, v, b in selected:
                value += u*v*local[i,j]*((b-a)*space).cos()/4
    return value


def verify():
    ctx.prec = 256
    reference_constants()
    scale = (2*arb(3).sqrt()-3)/(arb(3)/20).sqrt()
    mass = 8*arb.pi()*scale**3/arb(5).sqrt()*17/200
    lower = (arb(1824)/732)*arb(5).sqrt()/(8*arb.pi()*scale**3)
    print("交叉偏振常数与全几何非负所需的强度比下界", mass.str(12), lower.str(12))
    w, k4 = 2, 1
    for step in range(1, 10):
        w, k4 = 7*w+12*k4, 4*w+7*k4
        if step not in (3, 6, 9):
            continue
        k6 = w//2
        n = k4+k6
        t, four, six, _ = ground["matched_geometry"](n,k4)
        assert max(ground["critical_strength"](n,k4,four,1),
                   ground["critical_strength"](n,k6,six,1)) < 1
        b = 4*arb(k6)/(n*(4*t-1).sqrt())
        _, phase, prep = stability["simplex_hessians"](t,b,n)
        symbol = (8*arb.pi()*scale/(n*arb(5).sqrt()))*(phase+prep)
        ns = arb(n)
        zero, crossover, window = [static_symbol(symbol,s) for s in
                                   (arb(0), 1/ns.sqrt(), 1/ns.sqrt().sqrt())]
        assert abs(n*zero-mass) < arb(4)/n
        assert abs(n*crossover-(mass+arb(1)/10)) < arb(5)/n
        assert abs(ns.sqrt()*window-arb(1)/10) < 2/ns.sqrt()
        print("槽位、零动量与过渡点乘 N、窗口内读数乘 √N", n, (n*zero).str(12),
              (n*crossover).str(12), (ns.sqrt()*window).str(12))


if __name__ == "__main__":
    verify()
