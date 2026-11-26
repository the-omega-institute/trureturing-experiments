"""5432 槽驻点的交替形变、局部 Hessian 及任意周期盒上的正定界。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


EDGES = list(combinations(range(5), 2))


def gram(values):
    return arb_mat([[(values[0][i]+values[0][j]-values[i][j])/2
                     for j in range(1, 5)] for i in range(1, 5)])


def simplex_hessians(t, b, n):
    squared = [[arb(0) for _ in range(5)] for _ in range(5)]
    for i, j in EDGES:
        squared[i][j] = squared[j][i] = b*t if j-i in (1, 4) else b
    normals = arb_mat([[-1]*4] + [[int(i == j) for j in range(4)] for i in range(4)])
    frame = normals*gram(squared).inv()
    form = frame*normals.transpose()
    area_derivatives = arb_mat(10, 10)
    angle_derivatives = arb_mat(10, 10)
    shared_derivatives = arb_mat(10, 10)
    for col, (u, v) in enumerate(EDGES):
        dq = [[arb(0) for _ in range(5)] for _ in range(5)]
        dq[u][v] = dq[v][u] = arb(1)
        dform = -frame*gram(dq)*frame.transpose()
        for row, (i, j) in enumerate(EDGES):
            face = list(combinations([s for s in range(5) if s not in (i, j)], 2))
            six = j-i in (1, 4)
            area = b*((4*t-1) if six else t*(4-t)).sqrt()/4
            if (u, v) in face:
                area_derivatives[row,col] = (sum(squared[s][z] for s, z in face)/2
                                             - squared[u][v])/4
            denominator = (form[i,i]*form[j,j]).sqrt()
            cosine = -form[i,j]/denominator
            dcos = -dform[i,j]/denominator - cosine*(dform[i,i]/form[i,i]
                                                     + dform[j,j]/form[j,j])/2
            angle_derivatives[row,col] = -n*dcos/(16*arb.pi()*area*(1-cosine**2).sqrt())
            shared_derivatives[row,col] = area_derivatives[row,col]/(6 if six else 4)
    phase = area_derivatives.transpose()*angle_derivatives
    prep = area_derivatives.transpose()*shared_derivatives
    return squared, (phase + phase.transpose())/2, prep


def vertex_mode(squared, phase, prep, mask):
    phase_total, prep_total = arb(0), arb(0)
    for order in permutations(range(4)):
        vertices = [sum(1 << i for i in order[:j]) for j in range(5)]
        sign = [(-1)**((v & mask).bit_count()) for v in vertices]
        displacement = arb_mat([[squared[i][j]*(sign[i]+sign[j])] for i, j in EDGES])
        phase_total += (displacement.transpose()*phase*displacement)[0, 0]
        prep_total += (displacement.transpose()*prep*displacement)[0, 0]
    return phase_total, prep_total


def positive_lower_bound(matrix, lower):
    shifted = matrix - lower*(matrix**0)
    for size in range(1, matrix.nrows()+1):
        minor = arb_mat([[shifted[i,j] for j in range(size)] for i in range(size)])
        assert minor.det() > 0


def verify():
    ctx.prec = 192
    model = run_path(str(Path(__file__).with_name("ch39-clock-hinge-ground.py")))
    n = 5432
    t = model["matched_geometry"](n, 2911)[0]
    b = 4*arb(2521)/(n*(4*t-1).sqrt())
    squared, phase, prep = simplex_hessians(t, b, n)
    for mask in (0, 1, 15):
        angle, cost = vertex_mode(squared, phase, prep, mask)
        assert cost > 0
        if mask == 0:
            assert angle.contains(0)
        else:
            assert angle + cost < 0
        if mask == 15:
            assert arb("699.8682") < -angle/cost < arb("699.8684")
        print("交替掩码、角亏二阶项、四读二阶项、该方向的强度下界",
              mask, angle.str(14), cost.str(14), (-angle/cost).str(14))

    positive_lower_bound(phase + 927*prep, arb(1)/1000)
    print("强度 927 的每单纯形十维二阶型严格大于千分之一倍单位阵")
    zero_phase = arb_mat(15, 15)
    zero_prep = arb_mat(15, 15)
    for order in permutations(range(4)):
        vertices = [sum(1 << i for i in order[:j]) for j in range(5)]
        embed = arb_mat(10, 15)
        for row, (i, j) in enumerate(EDGES):
            embed[row, vertices[j]-vertices[i]-1] = 1
        zero_phase += embed.transpose()*phase*embed
        zero_prep += embed.transpose()*prep*embed
    positive_lower_bound(zero_phase + 699*zero_prep, arb(85))
    positive_lower_bound(zero_phase + zero_prep/2, arb(3)/50)
    print("强度 699 的零动量十五维二阶型严格大于 85 倍单位阵")
    print("强度 1/2 的零动量十五维二阶型严格大于 3/50 倍单位阵")


if __name__ == "__main__":
    verify()
