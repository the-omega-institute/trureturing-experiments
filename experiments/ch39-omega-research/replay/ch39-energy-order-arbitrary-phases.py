"""两壳荷六谱线按全谱能序一一配对，无论如何定相时的夸克振幅下界。"""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def polygon_floor(coefficients):
    lengths = [abs(z) for z in coefficients]
    return max(arb(0), (2*max(z.lower() for z in lengths)-sum(lengths)).lower())


def layer_floors(left, middle, right, source):
    low = min(w.lower() for w in source)
    floors = []
    pairs = list(combinations_with_replacement(range(19), 2))
    for j in range(19):
        floor = low
        for i, weight in enumerate(source):
            if middle is None:
                coefficients = [left[i,k]*right[k,j] for k in range(19)]
            else:
                coefficients = []
                for k,l in pairs:
                    first = left[i,k]*middle[k,l]*right[l,j]
                    if l != k:
                        first += left[i,l]*middle[l,k]*right[k,j]
                    coefficients.append(first)
            floor += (weight-low).lower()*polygon_floor(coefficients)**2
        floors.append(floor.lower())
    return floors


def verify(root):
    train = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    weighted = run_path(str(root / "ch39-common-unitary-quark-right-trace.py"))["weighted"]
    positive = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["positive"]
    n = 2_306_888_100_086
    _, labels, _ = train["verify_nodes"](root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 320
    s,t = [spectrum(shell,6,k)[1] for k in (n,n+1)]
    for charge in (-4,2):
        a,b = [spectrum(shell,charge,k)[1] for k in (n,n+1)]
        if charge == -4:
            # A* X² B = A* T D (S* T) D S* B.
            left, mid, right = adj(a)*t, adj(s)*t, adj(s)*b
        else:
            # A* X* B = A* S D* T* B.
            left, mid, right = adj(a)*s, None, adj(t)*b
        floors = layer_floors(left, mid, right, labels[n,charge])
        four = acb_mat([[b[i,j] for j in range(19)] for i in (14,8,7,6)])
        lower = weighted(four, floors)
        values = sorted(v.real.lower() for v in lower.eig(algorithm="rump"))
        print("荷、逐层界给出的四站升序谱下端", charge, [x.str(14) for x in values], flush=True)
        if charge == 2:
            ident = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
            positive(lower-ident*arb(".002627"))
        if charge == -4:
            limit = arb(".004976")
            ident = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
            positive(lower - ident*limit)
            distance = (limit.sqrt()-arb(".00000629").sqrt())/(2*max(w.upper() for w in labels[n,charge]).sqrt())
            assert distance>arb(".0349")
            print("达到轻上夸克窗所需的生成元距离下界",distance.lower().str(12),flush=True)


    # 一对调换时先检验较便宜的下型界，仅对余项计算上型的二次式。
    ud, vd = [spectrum(shell, 2, k)[1] for k in (n, n+1)]
    uu, vu = [spectrum(shell, -4, k)[1] for k in (n, n+1)]
    identity = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    pick = (14,8,7,6)
    fd = acb_mat([[vd[i,j] for j in range(19)] for i in pick])
    fu = acb_mat([[vu[i,j] for j in range(19)] for i in pick])
    down_count = 0
    up_count = 0
    for a,b in combinations(range(19),2):
        perm = list(range(19))
        perm[a],perm[b] = b,a
        swapped = acb_mat([[t[i,j] for j in perm] for i in range(19)])
        bounds = layer_floors(adj(ud)*s, None, adj(swapped)*vd, labels[n,2])
        constrained = weighted(fd, bounds)-identity*arb(".0000137")
        leading = [acb_mat([[constrained[i,j] for j in range(k)] for i in range(k)]).det().real for k in range(1,5)]
        if all(d>0 for d in leading):
            down_count += 1
            continue
        bounds = layer_floors(adj(uu)*swapped, adj(s)*swapped, adj(s)*vu, labels[n,-4])
        positive(weighted(fu,bounds)-identity*arb(".00000629"))
        up_count += 1
    assert down_count+up_count==171
    print("单次调换后由最轻下型、最轻上型排除的配对数",down_count,up_count,flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
