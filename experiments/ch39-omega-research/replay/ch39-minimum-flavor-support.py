"""在第 11431 枚零点，逐个分类能达到参考顶底比的最小局部剖面支撑。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


WITNESSES = {
    (0,3,9,10,11,12,13,15): [(602281841,0),(-599199628,-50900631),(3599241531,-615354617),
        (-3406892841,-1038158703),(-4496962309,1089188271),(4172363167,1969019067),
        (3864773914,-579029134),(-3657393041,-1383764691)],
    (1,2,7,8,9,13,14,15): [(86357338,0),(135990876,79489511),(741878194,770305657),
        (1760014544,860530205),(-2718473167,-1517649071),(-3037963270,-1347611343),
        (7001180964,3418262345),(-3229538628,-1765278591)]
}


def adjoint(matrix):
    return matrix.conjugate().transpose()


def flavor_maps(shell, index):
    height = acb.zeta_zero(index).imag
    bases = {}
    for q in (1,4,-2):
        values, modes = shell["laplacian"](q,height).eig(right=True)
        order = sorted(range(19),key=lambda j:float(values[j].real))
        assert all(values[a].real < values[b].real for a,b in zip(order,order[1:]))
        columns = []
        for j in order[:3]:
            v = [modes[k,j] for k in range(19)]
            norm = sum(abs(z)**2 for z in v).sqrt()
            columns.append([z/norm for z in v])
        bases[q] = columns
    maps = []
    for q in (4,-2):
        matrix = acb_mat([[bases[1][i][k].conjugate()*bases[q][j][k] for k in range(19)]
                          for i in range(3) for j in range(3)])
        maps.append(matrix.conjugate() if q == 4 else matrix)
    return maps


def orbit(sites, rotation):
    result = []
    for _ in range(3):
        result.append(tuple(sites))
        sites = sorted(rotation[k] for k in sites)
    return set(result)


def submatrix(matrix, group):
    return acb_mat([[matrix[i,j] for j in group] for i in group])


def largest(vector):
    matrix = acb_mat([[vector[3*i+j,0] for j in range(3)] for i in range(3)])
    return max(v.real for v in (adjoint(matrix)*matrix).eig(multiple=True)).sqrt()


def charm_bound(up, down, group, threshold):
    u = acb_mat([[up[i,j] for j in group] for i in range(9)])
    d = acb_mat([[down[i,j] for j in group] for i in range(9)])
    k = adjoint(d)*d - 3*adjoint(u)*u/threshold**2
    values, vectors = k.eig(right=True)
    order = sorted(range(len(group)),key=lambda j:float(values[j].real))
    assert values[order[0]].real < 0 and all(values[j].real > 0 for j in order[1:])
    assert all(values[a].real < values[b].real for a,b in zip(order,order[1:]))
    columns = []
    for j in order:
        vec = [vectors[i,j] for i in range(len(group))]
        norm = sum(abs(x)**2 for x in vec).sqrt()
        columns.append([x/norm for x in vec])
    ref = u*acb_mat([[x] for x in columns[0]])
    center = acb_mat([[ref[3*i+j,0] for j in range(3)] for i in range(3)])
    singular = sorted([v.real.sqrt() for v in (adjoint(center)*center).eig()],key=float)
    scaled = acb_mat([[columns[a][i]/values[order[a]].real.sqrt() for a in range(1,len(group))]
                      for i in range(len(group))])
    error_map = u*scaled
    max_sq = max(v.real for v in (adjoint(error_map)*error_map).eig(multiple=True))
    radius = (-values[order[0]].real*max_sq).sqrt()
    bound = (singular[1]-radius)/(singular[2]+radius)
    return bound


def verify():
    ctx.prec = 256
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    up,down = flavor_maps(shell,11431)
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    gu,gd = adjoint(up)*up,adjoint(down)*down
    # 所有七阶及更小的主子式为正给正定性；循环等价组只算一次。
    form7 = (arb("18.83")**2/3)*gd-gu
    form8 = (arb(60)**2/3)*gd-gu
    unknown = set()
    for size in range(1,9):
        form = form7 if size <= 7 else form8
        for group in combinations(range(19),size):
            cyclic = orbit(group,rotation)
            if group != min(cyclic):
                continue
            minor = submatrix(form,group).det()
            assert minor.imag.contains(0)
            if size == 8 and minor.real < 0:
                unknown.update(cyclic)
            else:
                assert minor.real > 0
    assert len(unknown) == 6 and unknown == set.union(*(orbit(g,rotation) for g in WITNESSES))
    for group, coefficients in WITNESSES.items():
        assert submatrix(gd,group).det().real > 0
        bound = charm_bound(up,down,group,arb(".8573")/arb(".01410"))
        assert bound > (arb(".64") if group[0] == 0 else arb(".559"))
        print("达到顶底层级时的粲顶比下界",bound.str(12))
        profile = acb_mat(19,1)
        for k,(re,im) in zip(group,coefficients):
            profile[k,0] = acb(re,im)/(10**10)
        norm = sum(abs(row[0])**2 for row in profile.tolist()).sqrt()
        ratio = largest(up*profile)/largest(down*profile)
        assert (arb("88.22") < ratio < arb("88.23") if group[0] == 0
                else arb("63.11") < ratio < arb("63.12"))
        print("支撑代表、例证顶底比、归一剖面给出的上型最大值",
              group,ratio.str(14),(largest(up*profile)/norm).str(12))
    print("至多七站上界、其余八站上界、达到参考层级的最小支撑数与配置数",
          18.83,60,8,len(unknown))


if __name__ == "__main__":
    verify()
