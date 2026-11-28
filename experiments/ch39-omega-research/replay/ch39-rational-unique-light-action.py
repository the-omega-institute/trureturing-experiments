"""用有限的七谐波消元与零点计数核验十二次有理多项式的唯一全局极小。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_poly, ctx


ROOTS = (
    "0.5090805312398923001015474683573496538267441645",
    "8.5982764736122745586265119637861770414543056840",
    "7.6724831982536457143459630603755122322333621347",
    "0.2552920596195782731987106532105009676172226764",
    "7.9838911542018547498230084673217413920774922646",
    "8.4486910135327225890032303515917146925000896755",
)
CENTERS = ((-4, 2345666), (2, 964890))


def det3(rows):
    a, b, c = rows
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def eliminate(rows, sines):
    # 三个余弦的 Cramer 分子与分母，满足环面 Gram 行列式恒等式。
    matrix = [row[1:4] for row in rows]
    denominator = det3(matrix)
    rhs = [-row[0]-sum(a*b for a,b in zip(row[4:], sines, strict=True)) for row in rows]
    numerators = [
        det3([[rhs[i] if j == axis else matrix[i][j] for j in range(3)] for i in range(3)])
        for axis in range(3)
    ]
    x, y, z = numerators
    polynomial = denominator*(x*x+y*y+z*z) - 2*x*y*z - denominator**3
    return polynomial, denominator, numerators


def fourier_row(matrix, edges, adj):
    def sample(a, b):
        link = a*edges[0] + b*edges[1]
        value = (matrix + link + adj(link)).det()
        assert value.is_finite() and value.imag.contains(0)
        return value.real

    corners = list(product((1, -1), repeat=2))
    values = [sample(x, y) for x, y in corners]
    even = [sum(x**a*y**b*v for (x,y),v in zip(corners,values,strict=True))/4
            for a,b in ((0,0), (1,0), (0,1), (1,1))]
    constant, cos_x, cos_y, cos_xy = even
    i = acb(0, 1)
    first = sample(i, 1) - constant - cos_y
    second = sample(1, i) - constant - cos_x
    third = (first-second-(sample(i, -i)-constant+cos_xy))/2
    return [*even, first-third, second+third, third]


def determinant_rows(pump, shell, radial, frames, charge, targets):
    base, links = pump["clock_matrices"](shell, charge)
    adj = pump["adjoint"]
    result = []
    for frame in frames:
        fixed = adj(frame)*base*frame
        compressed_radial = adj(frame)*radial*frame
        edges = [adj(frame)*link*frame for link in links]
        identity = acb_mat(fixed.nrows(), fixed.nrows())
        for i in range(identity.nrows()):
            identity[i, i] = 1
        rows = []
        for target in targets:
            matrix = fixed-target*identity
            at_zero = fourier_row(matrix, edges, adj)
            at_one = fourier_row(matrix-compressed_radial, edges, adj)
            # 每字符的径向项秩一。按行乘同一个正数只作数值调理。
            scale = max(abs(v).upper() for v in at_zero)
            rows.append([arb_poly([x/scale, (y-x)/scale])
                         for x,y in zip(at_zero, at_one, strict=True)])
        result.append(rows)
    return result


def contract_height(rows, sign, full, kappa, offset):
    interval = full
    sines = [arb(-1).union(arb(1))]*3
    for _ in range(5):
        _, denominator, numerators = eliminate(rows, sines)
        assert abs(denominator(interval)) > 0
        x,y,z = [p(interval)/denominator(interval) for p in numerators]
        assert -1 < x < 1
        sine_x = sign*(1-x*x).sqrt()
        sine_y = (z-x*y)/sine_x
        sines = [sine_x, sine_y, sine_x*y-sine_y*x]
        polynomial, _, _ = eliminate(rows, sines)
        mid = interval.mid()
        derivative = polynomial.derivative()(interval)
        assert abs(derivative) > 0
        interval = interval.intersection(mid-polynomial(mid)/derivative)
    assert interval > 0
    return (kappa/interval-offset).sqrt()


def verify():
    replay = Path(__file__).resolve().parent
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    selector = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = selector["light_spectra"](replay)
    levels = [e for q,n in CENTERS for e in light(q,n)]
    ctx.prec = 256
    roots = [arb(r) for r in ROOTS]
    eps = arb("1e-32")
    assert all(abs(x-y) < arb("1e-45") for x,y in zip(levels, roots, strict=True))
    assert all(abs(x-y) > arb(".14958") for i,x in enumerate(roots) for y in roots[i+1:])
    assert all(0 < x < 9 for x in [*levels, *roots])
    targets = [r+arb(0, eps.upper()) for r in roots]
    frames = list(pump["rotation_bases"](shell))
    radial = acb_mat(19, 19)
    pole = shell["soft"]["pole"]
    for site in shell["soft"]["star"]["orbit"]:
        radial[pole, pole] += 1
        radial[site, site] += 1
        radial[pole, site] -= 1
        radial[site, pole] -= 1
    sites = shell["geo"]["sites"]
    rotate = lambda i: sites.index((sites[i][1], sites[i][2], sites[i][0]))
    clock_targets = {j for i,j,d,_ in shell["geo"]["edges"] if d in (724,804)}
    assert clock_targets == {6,7,8} == {6,rotate(6),rotate(rotate(6))}
    assert all(i not in clock_targets for i,j,d,_ in shell["geo"]["edges"] if d in (724,804))
    tuned = set(shell["soft"]["star"]["orbit"])
    anchor = min(tuned)
    assert tuned == {anchor, rotate(anchor), rotate(rotate(anchor))} and len(tuned) == 3

    kappa = arb(23)/25*arb(2)**-20
    offset = shell["soft"]["point"]["study"]["offset"]
    assert kappa/(acb.zeta_zero(1).imag**2+offset) < arb("4e-9")
    full = arb(0).union(arb("4e-9"))
    wide = arb(-1).union(arb(1))
    for charge, index in CENTERS:
        atlas = determinant_rows(pump, shell, radial, frames, charge, targets)
        surviving = []
        for indices in product(range(6), repeat=3):
            rows = [block[i] for block,i in zip(atlas,indices,strict=True)]
            polynomial, _, _ = eliminate(rows, [wide]*3)
            if polynomial(full).contains(0):
                surviving.append(indices)
        expected = (0,1,2) if charge == -4 else (3,4,5)
        assert surviving == [expected]
        print("荷与 216 组中唯一保留的有序三元组", charge, expected)
        rows = [block[i] for block,i in zip(atlas, expected, strict=True)]
        for sign in (-1,1):
            height = contract_height(rows, sign, full, kappa, offset)
            lo = height.lower().zeta_nzeros().unique_fmpz()
            hi = height.upper().zeta_nzeros().unique_fmpz()
            assert lo is not None and hi is not None
            assert acb.zeta_zero(int(lo)).imag < height
            if (charge,sign) in ((-4,-1),(2,1)):
                assert (int(lo),int(hi)) == (index-1,index)
                assert height.contains(acb.zeta_zero(index).imag)
            else:
                assert lo == hi
            print("荷、正弦分支、高度与端点零点数", charge, sign, height.str(17), lo, hi)

    central_bound = 3 * arb("1e-45")**2 * 9**10
    competitor_bound = eps**2 * (arb(".14958")/2)**10
    assert central_bound < arb("1.047e-80") and competitor_bound-central_bound > arb("5.47e-76")
    print("中央作用上界与所有竞争作用下界", central_bound.str(9), competitor_bound.str(9))


if __name__ == "__main__":
    verify()
