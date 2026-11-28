"""以三个边界站点的五步子空间校正四步闭路迹，对近迹重谱作前向检验。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

ROOT = Path(__file__).resolve().parent
NEW = tuple(range(2_345_868, 2_345_963, 2))


def ritz_light_trace(matrix, pivot, sites, depth):
    # 每一步从前一层三个正交向量生成新列，逐列验证其范数正。
    layer = [acb_mat([[int(i==j)] for i in range(19)]) for j in sites]
    basis = []
    for step in range(depth+1):
        current = []
        for vec in layer:
            for previous in basis:
                vec -= previous*(previous.conjugate().transpose()*vec)[0,0]
            norm = sum(abs(x)**2 for x in vec.entries()).sqrt()
            assert norm > 0
            vec = vec*(1/norm)
            basis.append(vec)
            current.append(vec)
        if step < depth:
            layer = [matrix*v for v in current]
    assert len(basis) == 3*(depth+1)
    frame = acb_mat([[v[i,0] for v in basis] for i in range(19)])
    image = (matrix*matrix-matrix*pivot)*frame
    roots = (image.conjugate().transpose()*image).eig(multiple=True)
    assert len(roots)==len(basis) and all(z.imag.contains(0) for z in roots)
    return sum(sorted(z.real.lower() for z in roots)[:3]).union(
        sum(sorted(z.real.upper() for z in roots)[:3]))


def corrections(geometry, number, depth):
    h = acb.zeta_zero(number).imag
    offset = geometry["soft"]["point"]["study"]["offset"]
    pivot = 15-arb(23)/25*arb(2)**-20/(3*(h*h+offset))
    unit = arb("2.435e18")/(724**2*arb("1.22089e19"))
    sites = geometry["boundary_sites"]()
    total = arb(0)
    for q,g in ((1,-24),(4,-12),(2,-12),(3,-4),(6,-4)):
        moment = ritz_light_trace(geometry["laplacian"](q,h),pivot,sites,depth)
        total += g*((32 if q%2==0 else 1)*unit)**2*moment
    return total


def verify_window(load, numbers, depths):
    near = load("fifth-shell-near-traces")
    prediction = near["nearest"](load,numbers)
    involved = sorted({n for pair in prediction for n in pair})
    geometry = load("three-ray-chiral-boundary")
    ctx.prec = 320
    reads = [{n:corrections(geometry,n,d) for n in involved} for d in depths]
    matches = [0]*(len(depths)+1)
    intervals = load("fifth-shell-vacuum-order")["vacuum_intervals"]
    for a,b,_,change,_ in intervals(load, involved, selected_pairs=prediction):
        uncorrected = prediction[a,b]
        estimates = [uncorrected, *[uncorrected-(c[b]-c[a])/(-arb(103)/3) for c in reads]]
        assert not change.contains(0) and all(not v.contains(0) for v in estimates)
        for j,estimate in enumerate(estimates):
            matches[j] += int(estimate*change > 0)
        print(a,b,"全谱与修正预报、真实漂移区间",
              *[v.str(8) for v in estimates],change.lower().str(8),change.upper().str(8),flush=True)
    print("各方法命中数",matches,flush=True)
    return matches


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    near = load("fifth-shell-near-traces")
    assert verify_window(load, near["NUMBERS"], range(6)) == [3,4,1,4,4,3,5]
    assert verify_window(load, NEW, [5]) == [4,5]


if __name__=="__main__":
    verify()
