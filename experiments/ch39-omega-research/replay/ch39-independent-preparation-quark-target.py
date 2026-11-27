"""上下型各用一个制备背景时，以四个边界面的区间符号匹配两段夸克层级。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).resolve().parent


def lib(name):
    return run_path(str(ROOT/("ch39-"+name+".py")))


def verify():
    graph = lib("three-ray-chiral-boundary")
    pump = lib("native-light-projector-pump")
    ctx.prec = 160
    adj = pump["adjoint"]
    sources = lib("moving-preparation-quark-shape")["orthonormal_sources"](graph,adj)
    singular = lib("first-zero-relative-quark-spectrum")["singular_values"]
    ckm = lib("shell-clock-ckm-blocks")["standard_matrix"]
    v = acb_mat(ckm(*(arb(z).sin() for z in (".22704",".04275",".00377")),arb("1.139")))
    up = list(map(arb,(".00000615",".00311",".8616")))
    down = list(map(arb,(".0000135",".000268",".01401")))
    target = singular(acb_mat([[v[i,j]*up[i]/down[j] for j in range(3)] for i in range(3)]))
    goals = (target[1]/target[0],target[2]/target[1])
    settings = ((-4,"5.59475101",".87833265"),(2,".57241275","5.64213905"))
    width = arb(".0001")
    pieces = 512
    models = []
    for q,_,second in settings:
        base,links = pump["clock_matrices"](graph,q)
        sectors = []
        for frame in pump["rotation_bases"](graph):
            a,b = [adj(frame)*link*frame for link in links]
            term = b*acb(0,arb(second)).exp()
            sectors.append((adj(frame)*base*frame+term+adj(term),a,adj(frame)*sources))
        models.append(sectors)

    @lru_cache(maxsize=None)
    def interval(department,index,half_steps):
        angle = arb(settings[department][1])+width*index/pieces
        phase = acb(0,angle).exp()
        move = arb(2).sqrt()*width*half_steps/pieces + 6*arb("1e-8")
        rows = []
        minimum = None
        for fixed,a,inputs in models[department]:
            term = a*phase
            vec,gap = pump["first_line"](fixed+term+adj(term))
            assert gap > 2*move
            minimum = gap if minimum is None else min(minimum,gap)
            column = acb_mat([[z] for z in vec])
            rows.append((adj(column)*inputs).entries())
        matrix = acb_mat(rows)
        error = arb(2).sqrt()*move/(minimum-move) + arb("1.197e-8")
        relative_error = error/singular(matrix)[0]
        assert relative_error < 1
        return matrix, (1+relative_error)/(1-relative_error)

    for fixed in range(2):
        for side in (-1,1):
            lower_margin = arb(100)
            for position in range(-pieces+1,pieces,2):
                per_sector = [interval(k,side*pieces,0) if k==fixed else interval(k,position,1)
                              for k in range(2)]
                (mu,cu),(md,cd) = per_sector
                vals = singular(mu*md.inv())
                ratio = vals[fixed+1]/vals[fixed]
                distortion = cu*cd
                # 第一比沿上型角递增，第二比沿下型角递减，只须验证两端面的符号。
                above = (fixed==0 and side==1) or (fixed==1 and side==-1)
                margin = ratio/distortion-goals[fixed] if above else goals[fixed]-ratio*distortion
                assert margin > 0,(fixed,side,position,margin)
                lower_margin = min(lower_margin,margin.lower())
            print("固定角、边、相邻比相对中心目标的最小严格分离量",fixed,side,lower_margin.str(10))


if __name__ == "__main__":
    verify()
