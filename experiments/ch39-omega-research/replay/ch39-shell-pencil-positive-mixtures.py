"""检验壳铅笔的实、复单本征线见证和重根反例的有限维子式。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat


def has_minor(rows, size):
    for chosen_rows in combinations(rows,size):
        for cols in combinations(range(len(rows[0])),size):
            minor = arb_mat([[row[j] for j in cols] for row in chosen_rows])
            if not minor.det().contains(0):
                return True
    return False


def eigen_constraints(generators, vector, pivot, complex_line):
    column = acb_mat([[v] for v in vector])
    rows = []
    for k in range(3):
        if k == pivot:
            continue
        # 下述见证的 pivot 分量均为 1。
        differences = []
        for matrix in generators:
            image = matrix*column
            differences.append(image[k,0]-vector[k]*image[pivot,0])
        rows.append([z.real for z in differences])
        if complex_line:
            rows.append([z.imag for z in differences])
    return rows


def verify(replay):
    shell = run_path(str(Path(__file__).with_name("ch39-shell-pencil-six-mass-window.py")))
    orbit = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))
    model, _, _ = orbit["calibration"]()
    assert all(z.imag==0 for z in model["sources"].entries())
    _, h = shell["shell_frame"](model)
    a,b = h[0,1].real, h[1,2].real
    assert 0 < a < arb(".001") < arb(".9") < b
    assert all(h[i,i].contains(0) for i in range(3)) and h[0,2].contains(0)
    zero = acb_mat([[1,0,0],[0,0,0],[0,0,0]])
    cycle = acb_mat([[0,a,0],[a,0,b],[0,b,0]])
    square = acb_mat([[0,0,a*b],[0,a*a+b*b,0],[a*b,0,b*b]])
    real01 = acb_mat([[0,1,0],[1,0,0],[0,0,0]])
    j = acb(0,1)
    imag01 = acb_mat([[0,j,0],[-j,0,0],[0,0,0]])
    three = [zero,cycle,square]
    four = three+[real01]
    five = four+[imag01]
    for gens, line, pivot, complex_line in (
        (three,(1,0,1),0,False),(three,(0,1,1),1,False),
        (five,(1,j,1),0,True),(five,(0,1,j),1,True)):
        rows = eigen_constraints(gens,line,pivot,complex_line)
        assert has_minor(rows,len(gens)-1)
    # 两个稠密性见证均有 t=ab+b²，末二阶子式相同且非零。
    assert not ((a*a-a*b)*(-a*b)).contains(0)
    diagonal = acb_mat([[1,0,0],[0,2,0],[0,0,3]])
    generic = acb_mat([[0,2*a,0],[2*a,0,b],[0,b,0]])
    for target in (diagonal, generic):
        brackets = [v*target-target*v for v in four]
        rows = [[m[i,j].real for m in brackets] for i,j in ((0,1),(0,2),(1,2))]
        assert has_minor(rows,3)
    assert 4*a*a+b*b > 0
    print("三项/五项两组单谱见证与边界反例，四项两组交换子满秩见证均通过。")


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
