"""五荷总振幅的十五站局部反演及六个全域反例。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


CHARGES = (1, -4, 2, -3, 6)
# Ordered by omitted orbit; angles are physical radians before multiplication by charge.
CENTERS = (
    ("3.7213361701883096", "2.1550811181642398"),
    ("2.5618491369912753", "4.128104189015347"),
    ("4.0781832277397205", "3.7525702131162912"),
    ("2.476468441912506", "1.1734257799043757"),
    ("3.806716865267083", "5.109759527275216"),
    ("5.191272309661825", "1.4316311828028492"),
)


def adj(matrix):
    return matrix.conjugate().transpose()


def frobenius(matrix):
    return sum(z.abs_upper()**2 for z in matrix.entries()).sqrt()


def left_inverse_bound(matrix):
    mid = matrix.mid()
    left = (adj(mid)*mid).inv()*adj(mid)
    size = matrix.ncols()
    eye = acb_mat([[int(i == j) for j in range(size)] for i in range(size)])
    return (1-frobenius(eye-left*matrix)).lower()/frobenius(left).upper()


def eigenline_derivatives(matrix, slopes, first_line, pivot=5):
    entries, _ = first_line(matrix)
    n = len(entries)
    assert abs(entries[pivot]) > 0
    column = acb_mat([[z/entries[pivot]] for z in entries])
    column[pivot,0] = 1
    energy = (matrix*column)[pivot,0].real
    bordered = acb_mat(n+1,n+1)
    for i in range(n):
        for j in range(n):
            bordered[i,j] = matrix[i,j] - int(i == j)*energy
        bordered[i,n] = -column[i,0]
    bordered[n,pivot] = 1
    rhs = acb_mat(n+1,len(slopes))
    for j, slope in enumerate(slopes):
        product = slope*column
        for i in range(n):
            rhs[i,j] = -product[i,0]
    derivative = bordered.solve(rhs)
    return [column.entries(), *([derivative[i,j] for i in range(n)] for j in range(len(slopes)))]


def determinant_system(angles, rows, operators, first_line):
    columns = []
    for charge, (base, first, second) in zip(CHARGES, operators, strict=True):
        factor = acb(0,charge)
        a = first*(factor*angles[0]).exp()
        b = second*(factor*angles[1]).exp()
        slopes = [(item-adj(item))*factor for item in (a,b)]
        columns.append(eigenline_derivatives(base+a+adj(a)+b+adj(b), slopes, first_line))

    value, slopes = determinant_derivatives(columns, rows)
    return arb_mat([[value.real],[value.imag]]), arb_mat([[z.real for z in slopes], [z.imag for z in slopes]])


def determinant_derivatives(columns, rows):
    def minor(direction=0, replace=None):
        return acb_mat([[col[direction if j == replace else 0][i]
                         for j,col in enumerate(columns)] for i in rows]).det()

    slopes = [sum(minor(direction,j) for j in range(len(columns)))
              for direction in range(1,len(columns[0]))]
    return minor(), slopes


def verify(root):
    ctx.prec = 200
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    sites = shell["geo"]["sites"]
    rotation = [sites.index((b,c,a)) for a,b,c in sites]
    orbits = sorted({tuple(sorted({j,rotation[j],rotation[rotation[j]]})) for j in range(19)})
    assert [orbit for orbit in orbits if len(orbit) == 1] == [(14,)]
    triples = [orbit for orbit in orbits if len(orbit) == 3]
    assert len(triples) == 6
    height = acb.zeta_zero(2_306_888_100_085).imag
    laplacians = [shell["laplacian"](q,height) for q in CHARGES]
    frames = list(pump["rotation_bases"](shell))
    row_lists = []
    bounds = []
    for frame in frames:
        rows = [next(j for j in range(frame.ncols()) if abs(frame[orbit[0],j]) > 0) for orbit in triples]
        row_lists.append(rows)
        vectors = []
        for laplacian in laplacians:
            vector, gap = pump["first_line"](adj(frame)*laplacian*frame)
            assert gap > 1
            vectors.append(vector)
        synthesis = acb_mat(vectors).transpose()
        for selected in combinations(rows,5):
            part = acb_mat([[synthesis[i,j] for j in range(5)] for i in selected])
            bounds.append(left_inverse_bound(part))
    central = min(v.lower() for v in bounds)
    assert central > arb(".0002032")
    change = 4*arb(2).sqrt()*6*arb("1e-7")+6*arb("1e-8")
    local = central - arb(10).sqrt()*change/(1-change)
    assert local > arb(".00019")
    print("中央背景以及角半宽 1e-7、全软化区间的统一奇异值下界", central.str(12), local.str(12))

    frame = frames[0]
    operators = []
    for charge in CHARGES:
        base, links = pump["clock_matrices"](shell,charge)
        operators.append(tuple(adj(frame)*m*frame for m in (base,*links)))
    radius = arb("1e-11")
    eye = arb_mat([[1,0],[0,1]])
    for omitted, point in enumerate(CENTERS):
        rows = [j for k,j in enumerate(row_lists[0]) if k != omitted]
        center = list(map(arb,point))
        residual, jacobian = determinant_system(center, rows, operators, pump["first_line"])
        pre = jacobian.mid().inv()
        assert not pre.det().contains(0)
        box = [x+arb(0,radius.upper()) for x in center]
        _, whole = determinant_system(box, rows, operators, pump["first_line"])
        offset = pre*residual
        defect = eye-pre*whole
        ratio = max(offset[i,0].abs_upper()/radius + sum(z.abs_upper() for z in row)
                    for i,row in enumerate(defect.tolist()))
        assert ratio < arb(".0005")
        print("省去的轨道与失秩角盒自映比", triples[omitted], ratio.str(10))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())
