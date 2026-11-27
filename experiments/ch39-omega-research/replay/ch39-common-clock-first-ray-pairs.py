"""五荷恰有十对，逐对验证同一背景下三条原生射线的读出下界。"""
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).resolve().parent
TARGET = ".001"


def charge_reader(graph, pump, q):
    adjoint = pump["adjoint"]
    base, links = pump["clock_matrices"](graph, q)
    sectors = []
    for frame in pump["rotation_bases"](graph):
        orbit_index, = [j for j in range(frame.ncols()) if abs(frame[6,j]) > arb(".5")]
        blocks = [adjoint(frame)*link*frame for link in links]
        # 每一组钟边的首、末点不交。对任何向量 v，平方扰动残差可以按两端分别相加。
        squares = [adjoint(a)*a + a*adjoint(a) for a in blocks]
        sectors.append((adjoint(frame)*base*frame, blocks, squares, orbit_index))

    @lru_cache(maxsize=8192)
    def clock_cell(n, index):
        quarter, k = divmod(index,n)
        t = arb(2*k+1)/(2*n)
        def circle(x):
            return acb(1-x*x,2*x)/(1+x*x)
        center = circle(t)
        radius = max(abs(center-circle(arb(a)/n)).upper() for a in (k,k+1))
        z = (acb(0,1)**quarter * center)**q
        return z, min(arb(2), abs(q)*radius)

    def good(n, i, j):
        (u, ru), (v, rv) = clock_cell(n,i), clock_cell(n,j)
        shift = arb(2).sqrt()*(ru+rv) + 6*arb("1e-8")
        for fixed, (a,b), (sa,sb), idx in sectors:
            matrix = fixed+a*u+adjoint(a*u)+b*v+adjoint(b*v)
            vector, gap = pump["first_line"](matrix)
            if not gap > 2*shift:
                return False
            col = acb_mat([[z] for z in vector])
            error = 6*arb("1e-8")
            for square, radius in ((sa,ru),(sb,rv)):
                strength = (adjoint(col)*square*col)[0,0].real
                error += radius*strength.max(arb(0)).sqrt()
            sine = (error/(gap-shift)).upper()
            if not sine < 1:
                return False
            amplitude = abs(vector[idx]).lower()
            certified = amplitude*(1-sine**2).sqrt() - (1-amplitude**2).sqrt()*sine
            if not certified > arb(TARGET):
                return False
        return True
    return good


def certify_pair(pair):
    graph = run_path(str(ROOT/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(ROOT/"ch39-native-light-projector-pump.py"))
    assert sorted(graph["boundary_sites"]()) == [6,7,8]
    ctx.prec = 128
    readers = [charge_reader(graph,pump,q) for q in sorted(pair,key=abs)]
    pending = [(4,i,j) for i,j in product(range(16), repeat=2)]
    leaves = 0
    finest = 4
    while pending:
        n, i, j = pending.pop()
        if any(good(n,i,j) for good in readers):
            leaves += 1
            finest = max(finest,n)
        else:
            assert n < 4096, (pair,n,i,j)
            pending.extend((2*n,2*i+a,2*j+b) for a,b in product(range(2),repeat=2))
    return pair, leaves, finest


if __name__ == "__main__":
    pairs = list(combinations((1,-4,2,-3,6),2))
    with ProcessPoolExecutor(max_workers=len(pairs)) as pool:
        for pair, leaves, finest in pool.map(certify_pair,pairs):
            print("两荷、终格数、最细象限分数、至少一个荷的三站奇异值下界",
                  pair, leaves, finest, TARGET, flush=True)
