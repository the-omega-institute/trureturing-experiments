"""两枚实际零点下不可见的四站点与最少十六个可辨识重层约束。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def sine_between(a,b):
    norm2 = lambda v: sum((z*z.conjugate()).real for z in v)
    minors = [a[i]*b[j]-a[j]*b[i] for i,j in combinations(range(4),2)]
    return (norm2(minors)/(norm2(a)*norm2(b))).sqrt()


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    lift = run_path(str(directory/"ch39-four-site-heavy-lift.py"))
    ctx.prec = 320
    heights = [acb.zeta_zero(n).imag for n in (11431,11435)]
    pole, tuned = shell["soft"]["pole"],shell["soft"]["star"]["orbit"]
    edge_data = shell["geo"]["edges"]
    edges = {(i,j) for i,j,_,_ in edge_data}
    clocks = {(i,j) for i,j,d,_ in edge_data if d in (724,804) and pole not in (i,j)}
    variable = clocks | {tuple(sorted((pole,i))) for i in tuned}
    cells = [s for s in combinations(range(19),4) if all(e in edges for e in combinations(s,2))]
    blind = [s for s in cells if not any(i in s and j in s for i,j in variable)]
    blind_stars = [s for s in cells if not any(i in s or j in s for i,j in variable)]
    assert (len(blind),len(blind_stars)) == (925,126)
    chain = [(0,2,6,7),(0,2,6,8),(0,2,7,9),(0,2,7,10),(0,2,7,11),(0,2,7,13),
             (0,2,7,14),(0,2,7,16),(1,6,7,8),(1,3,6,7),(1,3,6,12),(1,3,6,15),
             (1,3,6,17),(4,6,7,8),(4,5,6,8),(4,5,6,18)]
    seen = set(chain[0])
    for s in chain[1:]:
        assert len(seen.intersection(s)) == 3
        seen.update(s)
    assert len(seen) == 19 and all(s in cells for s in chain)
    triangle = list(combinations(range(19),3))
    for q, threshold in ((2,".00137"),(4,".00202"),(6,".00556")):
        matrices = [shell["laplacian"](q,h) for h in heights]
        response = {}
        for t in triangle:
            if not all(e in edges for e in combinations(t,2)):
                continue
            i,j,k = t
            curves = [-m[i,j]*m[j,k]*m[k,i] for m in matrices]
            response[t] = abs(curves[0]/abs(curves[0])-curves[1]/abs(curves[1]))
        assert all(max(response[t].lower() for t in combinations(s,3))>1 for s in chain)
        rows = []
        for matrix in matrices:
            values, vectors = matrix.eig(right=True)
            order = sorted(range(19),key=lambda j:float(values[j].real))
            assert values[order[2]].real < 9 < values[order[3]].real
            basis = acb_mat([[vectors[i,j] for j in order[:3]] for i in range(19)])
            coeff = lift["light_minors"](basis,triangle)
            assert all(abs(z)>0 for z in coeff.values())
            rows.append({s:lift["cofactor_row"](s,coeff) for s in blind})
        smallest = min(sine_between(rows[0][s],rows[1][s]).lower() for s in blind)
        assert smallest > arb(threshold)
        print("荷、925 个局部相同四点组的湮灭行夹角正弦下界",q,smallest.str(11))


if __name__ == "__main__":
    verify()
