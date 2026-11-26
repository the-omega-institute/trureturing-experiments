"""用固定有限高度的一百个作用差唯一辨识共用同一真空的105个绕数权。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    lib = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    fixed = run_path(str(root / "ch39-height-independent-flavor-vacuum.py"))
    angles = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 320
    tri = lib["terms"](shell)
    x = arb("-.01046")
    y = x-angles["reference_angles"](flavor)[0]
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    orders = range(1,106)
    columns = [[arb(1), *fixed["gradient_columns"](lib, tri, q, x, y)] for q in orders]
    intervals = [lib["coeff"](tri, q, delta) for q in orders]
    included = [i for i in range(1,106) if i not in (1,3,5,6,53)]
    for i in included:
        for q, coeff, col in zip(orders,intervals,columns,strict=True):
            c = coeff[0]+coeff[2]
            phase = acb(0,arb.pi()*((q*i)%210)/105).exp()
            col.append((c*(1-phase)).real)
    controlled_bound = amplification(columns)
    print("100个作用差对各绕数权的最大绝对误差放大上界", controlled_bound.str(13))
    assert controlled_bound < arb("11.584")

    # 原作用的零相位常数也随 q 和软边变化，按每个高度重新计算全部529张三角面。
    columns = [column[:5] for column in columns]
    consecutive = [column[:] for column in columns]
    numbered_reads(lib, shell, tri, columns)
    numbered_bound = amplification(columns)
    print("100个实际零点处完整作用读数的最大绝对误差放大上界", numbered_bound.str(13))
    assert numbered_bound < arb("2.49")
    numbered_reads(lib, shell, tri, consecutive, consecutive_prefix=True)
    prefix_bound = amplification(consecutive)
    assert prefix_bound < arb("9.849")
    print("前100个连续编号的误差放大因子", prefix_bound.str(15))


def amplification(columns):
    matrix = arb_mat(columns).transpose()
    assert matrix.ncols()==matrix.nrows()==105
    inverse = matrix.inv()
    return max(sum(abs(inverse[i,j]) for j in range(5,105)).upper() for i in range(105))


def numbered_reads(lib, shell, tri, columns, *, consecutive_prefix=False):
    numbers = (
        2, 3, 11, 12, 15, 17, 18, 20, 25, 27, 28, 30, 34, 37, 39,
        41, 43, 44, 45, 46, 48, 49, 50, 51, 53, 54, 59, 62, 72, 73,
        75, 76, 93, 94, 95, 96, 97, 99, 104, 108, 114, 115, 116, 117, 119,
        120, 121, 122, 125, 127, 130, 131, 133, 135, 138, 140, 143, 145, 147, 150,
        152, 153, 155, 159, 162, 163, 169, 170, 171, 172, 174, 175, 180, 181, 182,
        184, 188, 189, 194, 196, 197, 200, 202, 205, 207, 212, 214, 217, 220, 221,
        223, 230, 234, 237, 240, 243, 245, 246, 248, 250,
    )
    if consecutive_prefix:
        numbers = range(1, 101)
    assert len(numbers)==len(set(numbers))==100 and min(numbers)>0
    base = []
    for q in range(1,106):
        polynomial = [arb(0), arb(0), arb(0)]
        for faces in tri.values():
            for z, softened in faces:
                assert softened in (0,1,2)
                value = 1-(z**q).real
                polynomial[0] += value
                polynomial[1] -= softened*value
                if softened==2:
                    polynomial[2] += value
        base.append(polynomial)
    unit = [lib["coeff"](tri,q,arb(0)) for q in range(1,106)]
    linear = [[sum((z**q for z,n in tri[k] if n==1),acb(0))
               for k in ((1,0),(0,1),(1,-1))] for q in range(1,106)]
    for number in numbers:
        h = acb.zeta_zero(number).imag
        phases = shell["geo"]["phases"](h)
        x, y = (phases[d].arg() for d in (724,804))
        delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
        for q,column in enumerate(columns,1):
            c = [u-delta*v for u,v in zip(unit[q-1],linear[q-1],strict=True)]
            a, b, d = base[q-1]
            column.append(a+delta*(b+delta*d)+lib["energy"](q,c,x,y))


if __name__=="__main__":
    verify()
