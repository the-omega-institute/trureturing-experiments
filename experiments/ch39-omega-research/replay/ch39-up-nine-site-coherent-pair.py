"""九站特殊支撑上，以第 16、18 目的层的相干差约束轻上夸克。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


SUPPORT = (0, 1, 2, 3, 4, 5, 11, 12, 14)


def adj(matrix):
    return matrix.conjugate().transpose()


def column(matrix, j):
    return acb_mat([[matrix[i,j]] for i in range(matrix.nrows())])


def verify(root):
    local = run_path(str(root / "ch39-up-eight-site-cover.py"))
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    positive = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["positive"]
    weighted = run_path(str(root / "ch39-common-unitary-quark-right-trace.py"))["weighted"]
    n = 2_306_888_100_086
    _, labels, _ = training(root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    source = spectrum(shell,-4,n)[1]
    target = spectrum(shell,-4,n+1)[1]
    assert min(labels[n,-4]) > 0
    h = weighted(source, labels[n,-4])
    block = acb_mat([[h[i,j] for j in SUPPORT] for i in SUPPORT])
    positive(block)
    rows = local["rows"]
    inside = rows(target, SUPPORT)
    outside = acb_mat(target)
    for i in SUPPORT:
        for j in range(19):
            outside[i,j] = 0
    mixed = h*outside
    rhs = rows(mixed, SUPPORT)

    def bound(z, out, hout, b, gamma):
        shifted = acb_mat(block)
        for i in range(len(SUPPORT)):
            shifted[i,i] += gamma
        q = (adj(z)*z)[0,0].real
        return (adj(out)*hout-adj(b)*shifted.solve(b))[0,0].real - gamma*q

    # 取支撑内两个目的矢的相干差，使其在九站内的范数最小。
    inner = (adj(column(inside,15))*column(inside,17))[0,0]
    assert abs(inner) > 0
    phase = inner.conjugate()/abs(inner)

    def difference(matrix):
        return (column(matrix,15)-phase*column(matrix,17))*(1/arb(2).sqrt())

    delta = bound(difference(inside),difference(outside),difference(mixed),
                  difference(rhs),arb(".015"))
    assert delta > arb(".00022")
    print("第16、18层相干差在九站内的概率及振幅和下界",
          (adj(difference(inside))*difference(inside))[0,0].real.str(15), delta.str(15),flush=True)

    single = []
    for j in range(19):
        gamma = arb(".1") if j == 9 else arb(0)
        value = bound(column(inside,j),column(outside,j),column(mixed,j),column(rhs,j),gamma)
        single.append(max(arb(0),value.lower()))
    four = rows(target,(14,8,7,6))
    for j, floor in ((15,arb(".000135")), (17,arb(".000085"))):
        weights = single.copy()
        weights[j] = max(floor,weights[j])
        lower = weighted(four,weights)
        for i in range(4):
            lower[i,i] -= arb(".000048")
        positive(lower)
        print("分支中的层号、该层振幅下界、四站最小本征值严格下界",
              j+1, floor.str(9), "0.000048",flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
