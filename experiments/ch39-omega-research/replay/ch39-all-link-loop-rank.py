"""全连接相位的一阶力，按极点起算的度量回路分类。"""
from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_poly, arb_mat, fmpz_mat

phase = run_path(str(Path(__file__).with_name("ch39-matter-number-three-spokes.py")))
soft = phase["soft"]
pole, orbit = soft["pole"], soft["star"]["orbit"]
rotation = soft["star"]["charged"]["rotation"]


def loop_table():
    ds = {(i,j):d for i,j,d,_ in soft["geo"]["edges"]}
    vertices = set(range(19))-{pole}
    reps = {}
    for i,j in sorted(ds):
        if pole in (i,j) or (i,j) in reps:
            continue
        u,v = i,j
        for _ in range(3):
            reps[u,v] = (i,j),1
            reps[v,u] = (i,j),-1
            u,v = rotation[u],rotation[v]
    rows = sorted({r for r,sign in reps.values()})
    tables = defaultdict(Counter)
    for i in orbit:
        for j in vertices-set(orbit):
            for k in range(4):
                for interior in permutations(sorted(vertices-{i,j}),k):
                    path = (i,*interior,j)
                    cycle = (pole,*path,pole)
                    shape = tuple(ds.get(tuple(sorted((u,v)))) for u,v in zip(cycle,cycle[1:]))
                    if None in shape:
                        continue
                    shape = min(shape,shape[::-1])
                    for u,v in zip(path,path[1:]):
                        rep,sign = reps[u,v]
                        tables[shape][rep] += sign
    assert all(x%3 == 0 for col in tables.values() for x in col.values())
    cols = sorted(tables, key=lambda shape:(len(shape),shape))
    return rows,cols,fmpz_mat([[tables[c][r]//3 for c in cols] for r in rows])


def free_forces(matrices, weight, c, rows):
    sine = (1-c*c).sqrt()
    def determinant(m):
        a,b,d = phase["spoke_polynomials"](m,pole,orbit,weight,acb(1))
        return a+b*c+d*sine
    denominators = {q:determinant(m) for q,m in matrices.items()}
    values = []
    for i,j in rows:
        total = arb(0)
        for q,m in matrices.items():
            moved = []
            for sgn in (1,-1):
                twist = acb_mat(m)
                twist[i,j] *= acb(0,sgn)
                twist[j,i] = twist[i,j].conjugate()
                moved.append(determinant(twist))
            # Fourier degree one in a single oriented edge.
            numerator = (moved[1]-moved[0])*arb(".5")
            assert all(numerator[k].contains(0) for k in (17,18,19))
            numerator = arb_poly(numerator.coeffs()[:17])
            total += {0:33,5:18,10:3}[abs(q)]*phase["integrate_force"](
                numerator,denominators[q],"2e-8")
        values.append(total/(8*arb.pi()**2*sine))
    return values


def verify():
    matrices, weight, c, radial_weights, _ = phase["verify"]()
    rows, cols, mat = loop_table()
    n4 = sum(len(shape)<=4 for shape in cols)
    a4 = fmpz_mat([[mat[i,j] for j in range(n4)] for i in range(len(rows))])
    assert len(rows)==39 and n4==92 and a4.rank()==34
    n5 = sum(len(shape)<=5 for shape in cols)
    a5 = fmpz_mat([[mat[i,j] for j in range(n5)] for i in range(len(rows))])
    assert n5 == 871 and a5.rank() == 34 and len(cols) == 8039 and mat.rank() == 39
    print("Ranks through four, five and six links",34,34,39)
    values = free_forces(matrices,weight,c,rows)
    for i in ((6,7),(9,15),(10,13),(11,12)):
        p = rows.index(i)
        assert all(a5[p,j] == 0 for j in range(n5))
        print("Blind directed edge, free derivative / sine",i,values[p].str(9))
    assert values[rows.index((9,15))] < -arb(".0001")
    chosen = [1020,1448,1584,1585,1616,1847,1932,1956,1978,2044,2964,3008,3379,
              4099,4104,4927,5133,5146,5523,5524,5550,5659,5667,5697,5782,5860,
              5867,5894,6012,6532,7155,7205,7239,7391,7445,7580,7686,7790,7872]
    basis = arb_mat([[int(mat[i,j]) for j in chosen] for i in range(len(rows))])
    weights = basis.solve(arb_mat([[-v] for v in values]))
    floor = min(v.lower() for v in weights.entries())
    print("Least certified six-link basis weight",floor.str(10))
    assert floor > arb("6e-6") and all(len(cols[j]) == 6 for j in chosen)
    radial = run_path(str(Path(__file__).with_name("ch39-radial-phase-stability.py")))
    radial["classify_faces"]()
    sites = soft["geo"]["sites"]
    mirror = radial["mirror"]
    reflection = [sites.index((x,z,y)) for x,y,z in sites]
    assert {reflection[i] for i in orbit} == mirror and reflection[pole] == pole
    ds = {(i,j):d for i,j,d,_ in soft["geo"]["edges"]}
    assert all(ds[tuple(sorted((reflection[i],reflection[j])))] == d for (i,j),d in ds.items())
    chords = [edge for edge in ds if pole not in edge]
    cycles = []
    for i,j,k in combinations(sorted(set(range(len(sites)))-{pole}),3):
        if all(e in ds for e in ((i,j),(j,k),(i,k))):
            row = [0]*len(chords)
            for e,v in (((i,j),1),((j,k),1),((i,k),-1)):
                row[chords.index(e)] = v
            cycles.append(row)
    assert len(chords)==117 and len(cycles)==412 and fmpz_mat(cycles).rank()==100

    # Stationarity fixes K_AB via the mirror group force and K_AB+K_AC via
    # the active group force. Reflection of the metric Wilson action gives
    # K_BC=K_AC. Evaluate the two remaining physical curvatures directly.
    active = radial["active"]
    kaa = 3*sum(radial_weights)
    kab = 3*radial_weights[0]
    kbulk = kaa-kab
    ua,ub,uab = (radial["free_curvature"](matrices,weight,c,v)
                 for v in (active,mirror,active|mirror))
    h1, h2 = ua+c*kaa, ub+c*kab+kbulk
    cross = (uab-ua-ub)/2-c*kab
    det = h1*h2-cross**2
    assert h1 > arb(".017") and det > arb(".0024")
    print("Flat nonpole triangle rank, residual physical directions and reduced determinant",
          100,2,det.str(8))
    return rows,cols,mat,values


if __name__ == "__main__":
    verify()
