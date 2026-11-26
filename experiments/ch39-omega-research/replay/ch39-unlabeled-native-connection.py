"""以未编号的平方距离和一个手性选择恢复原磁拉普拉斯的边、钟、运输与软辐边。"""
from collections import Counter
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx, fmpq, fmpq_mat


def reconstruct(D, oriented_far):
    n = len(D)
    profiles = [tuple(sorted(row)) for row in D]
    pole, = [i for i, s in enumerate(profiles) if profiles.count(s) == 1]
    radii = sorted(set(D[pole])-{0})
    near = [i for i in range(n) if D[pole][i] == radii[0]]
    far = [i for i in range(n) if D[pole][i] == radii[-1]]
    assert len(near) == 3 and len(far) == 6 and oriented_far in far
    base = min(D[i][j] for i,j in combinations(far,2))
    refs, long = set(), set()
    iso_edges = set()
    for i,j in combinations(far,2):
        if D[i][j] != base:
            continue
        iso, = [k for k in near if D[i][k] == D[j][k]]
        refs.add(D[i][iso])
        iso_edges.update(tuple(sorted((k,iso))) for k in (i,j))
        longer = [tuple(sorted((D[i][k],D[j][k]))) for k in near if k != iso]
        assert len(set(longer)) == 1
        long.add(longer[0])
    ref, = refs
    (short, high), = long
    edges = set((i,j) for i,j in combinations(range(n),2) if D[i][j] <= high)
    assert iso_edges == {e for e in edges if D[e[0]][e[1]] == ref}
    third = Counter()
    mixed_edges = set()
    for tri in combinations(range(n),3):
        sides = list(combinations(tri,2))
        lengths = [D[i][j] for i,j in sides]
        if set(sides) <= edges and short in lengths and high in lengths:
            k, = [k for k, d in enumerate(lengths) if d not in (short,high)]
            third[lengths[k]] += 1
            mixed_edges.add(sides[k])
    assert len(third) == 2 and set(third.values()) == {6}
    negative = max(third)
    assert len([e for e in mixed_edges if D[e[0]][e[1]] == negative]) == 3
    ordered = {f: tuple(sorted(near,key=lambda p: D[f][p])) for f in far}
    assert all(tuple(D[f][i] for i in ns) == (ref, short, high) for f,ns in ordered.items())
    basis = ordered[oriented_far]
    def parity(ns):
        perm = [basis.index(x) for x in ns]
        return sum(perm[i] > perm[j] for i,j in combinations(range(3),2)) % 2
    soft = {f for f, ns in ordered.items() if not parity(ns)}
    assert len(soft) == 3
    # 以极点和三个有序近点作仿射基，只用距离计算 Gram 矩阵与球心。
    def gram(i,j):
        return fmpq(D[pole][i]+D[pole][j]-D[i][j],2)
    G = fmpq_mat([[gram(i,j) for j in basis] for i in basis])
    inv = G.inv()
    center = inv*fmpq_mat([[gram(i,i)/2] for i in basis])
    pos = [inv*fmpq_mat([[gram(i,j)] for i in basis])-center for j in range(n)]
    def dot(a,b): return (a.transpose()*G*b)[0,0]
    def det(a,b,c): return fmpq_mat([a.entries(),b.entries(),c.entries()]).det()
    R = dot(pos[pole],pos[pole])
    assert all(dot(v,v) == R for v in pos)
    native = {}
    for i,j in edges:
        a,b,c = pos[pole],pos[i],pos[j]
        re = R*(R+dot(a,b)+dot(b,c)+dot(c,a))
        im = -det(a,b,c)
        amplitude = acb(arb(re), arb(R*G.det()).sqrt()*arb(im))
        sign = -1 if D[i][j] == negative else 1
        distance = D[i][j]
        clock = 1 if distance == short else 2 if distance == high else 0
        sense = (1 if i in far else -1) if clock else 0
        native[i,j] = (sign*amplitude/abs(amplitude), clock, sense, pole in (i,j) and bool(soft & {i,j}))
    return R, (ref,short,high), pole, soft, negative, native


def clocks(geometry, height):
    R, sides = geometry[:2]
    angles = [(1-arb(d)/(2*arb(R))).acos() for d in sides]
    return [acb(1), *(acb(0,height*(a-angles[0])).exp() for a in angles[1:])]


def matrix(geometry, q, phi, delta):
    L = acb_mat(19,19)
    for (i,j),(z,clock,sense,soft) in geometry[-1].items():
        w = 1-delta if soft else arb(1)
        rotation = phi[clock] if sense >= 0 else phi[clock].conjugate()
        U = (z*rotation)**q
        L[i,i] += w
        L[j,j] += w
        L[i,j] = -w*U
        L[j,i] = L[i,j].conjugate()
    return L


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    graph = shell["geo"]
    ctx.prec = 224
    sites = graph["sites"]
    D = [[sum((a-b)**2 for a,b in zip(p,v)) for v in sites] for p in sites]
    pos, neg = (reconstruct(D, witness) for witness in (0,1))
    # 再用完全相反的编号，检验远到近的钟方向与恢复的取向随标签一同搬运。
    reverse = list(range(18,-1,-1))
    relabeled = reconstruct([[D[i][j] for j in reverse] for i in reverse],reverse.index(0))
    R, sides, pole, selected, minus, native = pos
    assert (R,sides,minus) == (867,(244,724,804),200)
    assert selected == set(shell["soft"]["star"]["orbit"]) == {0,3,4}
    assert neg[3] == {1,2,5} and pole == shell["soft"]["pole"]
    assert len(native) == 135 and sum(D[i][j] == minus for i,j in native) == 6
    print("参考、负号、两钟的平方弦长与各类边数",[(d,sum(D[i][j]==d for i,j in native)) for d in (sides[0],minus,*sides[1:])])
    print("两个可选的保向远点三元轨道",sorted(selected),sorted(neg[3]))
    for i,j,d,z in graph["edges"]:
        recovered, clock, sense, softened = native[i,j]
        assert (recovered-z).contains(0)
        assert clock == (1 if d == 724 else 2 if d == 804 else 0)
        assert sense == int(bool(clock))
        assert softened == (pole in (i,j) and (i in selected or j in selected))
    height = acb.zeta_zero(2_306_888_100_085).imag
    phi = clocks(pos,height)
    original_clocks = graph["phases"](height)
    assert all((v-original_clocks[d]).contains(0) for d, v in zip(sides[1:],phi[1:]))
    assert all(not v.imag.contains(0) for v in phi[1:])
    offset = shell["soft"]["point"]["study"]["offset"]
    delta = arb(23)/25*arb(2)**-20/(height**2+offset)
    reflection = [sites.index((y,x,z)) for x,y,z in sites]
    charges = (0,1,-1,2,-2,3,-3,4,-4,6,-6)
    for q in charges:
        L = matrix(pos,q,phi,delta)
        old = shell["laplacian"](q,height)
        renamed = matrix(relabeled,q,phi,delta)
        assert all((L[i,j]-renamed[18-i,18-j]).contains(0) for i in range(19) for j in range(19))
        mirrored = matrix(neg,q,phi,delta)
        assert all((L[i,j]-old[i,j]).contains(0) for i in range(19) for j in range(19))
        assert all((L[i,j]-mirrored[reflection[i],reflection[j]]).contains(0)
                   for i in range(19) for j in range(19))
    # 若只换软三元组但维持空间中的正向，反射等于对“钟也反向”的矩阵取共轭。
    other_edges = {ij:(pos[-1][ij][0], k, s, soft) for ij, (_,k,s,soft) in neg[-1].items()}
    fixed_orientation = (*neg[:-1], other_edges)
    reversed_clocks = [v.conjugate() for v in phi]
    other = matrix(fixed_orientation,1,reversed_clocks,delta)
    L = matrix(pos,1,phi,delta)
    assert all((L[i,j]-other[reflection[i],reflection[j]].conjugate()).contains(0)
               for i in range(19) for j in range(19))
    unchecked_clock_flip = matrix(fixed_orientation,1,phi,delta)
    assert any(not (L[i,j]-unchecked_clock_flip[reflection[i],reflection[j]].conjugate()).contains(0)
               for i,j in native)
    print("联合零点上逐项复核的正负荷与中性荷",charges)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
