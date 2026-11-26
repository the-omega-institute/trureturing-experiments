"""在平方弦距截止 244 上用一个五维外点实现和一张二维平面补全三时钟响应。"""
from functools import cache
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

pairs = run_path(str(Path(__file__).with_name("ch39-three-clock-outer-pairs.py")))
previous = pairs["previous"]
outer, centre = pairs["outer"], pairs["centre"]


@cache
def shell_data():
    _, _, _, g0, g1, g2 = previous["MEM"]["coefficients"]()
    delta = 10**10*(g2-g1*g0.solve(g1))
    local = previous["LOCAL"]
    points = sorted(local["shell_grid"](1) | {(17,17,17)})
    distance = [[local["squared"](p,q) for q in points] for p in points]
    initial = pairs["projected_outer_gram"](delta,0)
    slope = pairs["projected_outer_gram"](delta,1)-initial
    return points,distance,g0,delta,initial,slope


def prepare(t):
    points,distance,g0,delta,initial,slope = shell_data()
    matrix = arb_mat if isinstance(t,arb) else fmpq_mat
    gram = matrix(initial)+t*matrix(slope)
    metric = matrix(7,7)
    for i in range(6):
        for j in range(6):
            metric[i,j] = gram[i,j]
    metric[6,6] = t
    base, direction = matrix(19,7), matrix(19,7)
    for k,i in enumerate(outer):
        base[i,k] = base[i,6] = 1
    base[centre,6] = delta[centre,outer[0]]/t
    for i in set(range(19))-set(outer)-{centre}:
        fixed = [k for k,j in enumerate(outer) if distance[i][j]>296]
        assert len(fixed)==4
        rhs = matrix([[delta[i,outer[k]],-t] for k in fixed])
        columns = matrix([[gram[k,l] for l in fixed] for k in fixed]).solve(rhs)
        for j,k in enumerate(fixed):
            base[i,k], direction[i,k] = columns[j,0], columns[j,1]
        direction[i,6] = 1
    aa = base*metric*base.transpose()
    ba = direction*metric*base.transpose()
    bb = direction*metric*direction.transpose()
    def entries(i,j):
        return aa[i,j]-delta[i,j],ba[i,j],ba[j,i],bb[i,j]
    return points,distance,g0,delta,entries


def verify():
    ctx.prec = 192
    points,distance,g0,delta,entries = prepare(fmpq(1))
    def residual(i,j,x,y):
        a,b,c,d = entries(i,j)
        return a+b*x+c*y+d*x*y

    # 296 处解除的六条近外点约束在本构造中也须精确满足。
    a,b,_,_ = entries(4,2)
    x = -a/b
    A,B = [residual(4,j,x,x) for j in (12,15)]
    v = (arb(2*A-B)/3).sqrt()
    u = (-v-arb(2*A+B).sqrt())/2
    assert v>0 and u<0
    def linear(i,j,s):
        a,b,c,d = entries(i,j)
        return a+b*s,c+d*s
    c0,c1,c2,c3 = entries(6,7)
    a,b = linear(4,16,x)
    # r=(a+b*y)/(3u) 和 3r²=c0+(c1+c2)y+c3*y²。
    A2 = b*b-3*u*u*c3
    A1 = 2*a*b-3*u*u*(c1+c2)
    A0 = a*a-3*u*u*c0
    disc = A1*A1-4*A2*A0
    assert disc>0 and 0 not in A2
    roots = [(-A1+sign*disc.sqrt())/(2*A2) for sign in (-1,1)]
    y, = [value for value in roots if -3<value<-2]
    r = (a+b*y)/(3*u)
    a4,b4 = linear(4,14,x)
    a6,b6 = linear(6,9,y)
    matrix = arb_mat([[b4,-3*u],[b6,-3*r]])
    assert 0 not in matrix.det()
    result = matrix.solve(arb_mat([[-a4],[-a6]]))
    w,s = result[0,0],result[1,0]

    permutations3 = list(permutations(range(3)))
    group = [[points.index(tuple(p[k] for k in perm)) for p in points] for perm in permutations3]
    # 所有指定内点方程均由以下六条在坐标置换下生成。
    representatives = [(4,12),(4,15),(6,7),(4,16),(4,14),(6,9)]
    inner = set(range(19))-set(outer)-{centre}
    edges = {tuple(sorted((g[i],g[j]))) for g in group for i,j in representatives}
    assert edges == {(i,j) for i,j in combinations(sorted(inner),2) if distance[i][j]>244}
    for g in group:
        assert all(entries(i,j)==entries(g[i],g[j]) for i,j in combinations(range(19),2))
    # 第三坐标由零和约束决定，三列只增加两维。用整数系数验证六类边的平面内积。
    templates = {4:((1,0,0,0),(0,1,0,0),(-1,-1,0,0)),
                 6:((0,0,1,0),(0,0,1,0),(0,0,-2,0)),
                 8:((0,0,0,1),(0,0,0,1),(0,0,0,-2))}
    linear_rows = {}
    parameters = [arb(0)]*19
    for i,template in templates.items():
        for perm,g in zip(permutations3,group):
            coeffs = fmpq_mat([template[k] for k in perm])
            if g[i] in linear_rows:
                assert coeffs==linear_rows[g[i]]
            linear_rows[g[i]] = coeffs
            parameters[g[i]] = {4:arb(x),6:y,8:w}[i]
    monomials = [(i,j) for i in range(4) for j in range(i,4)]
    expected = [{(0,0):-1,(0,1):-1,(1,1):-1}, {(1,1):1,(0,0):-2,(0,1):-2},
                {(2,2):-3}, {(0,2):-3}, {(0,3):-3}, {(2,3):-3}]
    for (i,j),target in zip(representatives,expected,strict=True):
        m = linear_rows[i].transpose()*linear_rows[j]
        assert all(m[k,l]+(m[l,k] if k!=l else 0)==target.get((k,l),0) for k,l in monomials)

    for i,j in combinations(range(19),2):
        if distance[i][j]<=244 or (i,j) in edges:
            continue
        if i in inner or j in inner:
            if j in inner:
                i,j = j,i
            a,b,c,d = entries(i,j)
            assert not c and not d
            assert (not a and not b) or (i in {g[4] for g in group} and a+b*x==0)
        else:
            assert entries(i,j)[0]==0
    vector = arb_mat([[u],[v],[r],[s]])
    plane = [arb_mat(linear_rows[i])*vector if i in inner else arb_mat(3,1) for i in range(19)]
    row_sums = []
    for i in range(19):
        row_sum = arb(0)
        for j in range(19):
            base = delta[i,j]+residual(i,j,parameters[i],parameters[j])
            row_sum += abs(base+(plane[i].transpose()*plane[j])[0,0])
        row_sums.append(row_sum.upper()/10**10)
    norm = max(row_sums)
    inv_norm = max(sum(abs(x) for x in row) for row in g0.inv().tolist())
    assert norm*inv_norm < (arb(1)-arb(".098"))**2
    print("Cutoff 244, auxiliary dimension at most 26, maximum Gram row sum",norm.str(12),
          "center coordinates",float(x),y.str(10),w.str(10),flush=True)


if __name__ == "__main__":
    verify()
