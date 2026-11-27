"""在十站核心外的三个循环轨道中各取一站，排除全部二十七个十三站支撑及其子集。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


CORE = (0,1,2,3,4,5,11,12,14,18)
REST = ((6,7,8),(9,15,16),(10,13,17))


def verify(root):
    api = run_path(str(root/"ch39-up-eight-site-cover.py"))
    pair = run_path(str(root/"ch39-up-ten-site-coherent-exceptions.py"))
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    spectra = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    weighted = run_path(str(root/"ch39-common-unitary-quark-right-trace.py"))["weighted"]
    positive = run_path(str(root/"ch39-up-quark-arbitrary-site-phase-floor.py"))["positive"]
    rows, cols, adj = api["rows"], pair["columns"], pair["adj"]
    n = 2_306_888_100_086
    _,labels,_ = train["verify_nodes"](root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    source, target = [spectra(shell,-4,k)[1] for k in (n,n+1)]
    assert min(labels[n,-4]) > 0
    h = weighted(source,labels[n,-4])
    four = rows(target,(14,8,7,6))
    sets = [set() for _ in (11,12,13)]
    least_level, least_gain = arb(1), arb(1)
    for extra in product(*REST):
        support = tuple(sorted(CORE+extra))
        inside = rows(target,support)
        outside = acb_mat(target)
        for site in support:
            for j in range(19):
                outside[site,j] = 0
        mixed = h*outside
        block = acb_mat([[h[i,j] for j in support] for i in support])
        rhs = rows(mixed,support)
        gram, constant = adj(inside)*inside, adj(outside)*mixed

        def schur(ids, gamma):
            shifted = acb_mat(block)
            for i in range(len(support)):
                shifted[i,i] += gamma
            positive(shifted)
            coupling = cols(rhs,ids)
            direct = acb_mat([[constant[i,j]-gamma*gram[i,j] for j in ids] for i in ids])
            return direct - adj(coupling)*shifted.solve(coupling)

        gammas = api["multipliers"](block,rhs,[gram[j,j].real for j in range(19)])
        floors = [max(arb(0),schur([j],gamma)[0,0].real.lower()) for j,gamma in enumerate(gammas)]

        def alternatives(parent, j, k, gamma, high_j, high_k):
            nonlocal least_gain
            form = schur((j-1,k-1),arb(gamma))
            r, s = arb(high_j),arb(high_k)
            a, b = r-form[0,0].real,s-form[1,1].real
            assert a > 0 and b > 0
            gain = abs(form[0,1])-(r*s).sqrt()-(a*b).sqrt()
            assert gain > 0
            least_gain = min(least_gain,gain.lower())
            first,second = parent.copy(),parent.copy()
            first[j-1] = max(first[j-1],r)
            second[k-1] = max(second[k-1],s)
            return first,second

        a,b = alternatives(floors,1,10,".09",".001",".003")
        leaves = (alternatives(a,10,16,".016",".00075",".000025")
                  + alternatives(b,16,18,".0039",".00002",".000013"))
        for guaranteed in leaves:
            matrix = weighted(four,guaranteed)
            values = matrix.eig(algorithm="rump")
            assert all(v.imag.contains(0) for v in values)
            level = min(v.real.lower() for v in values)
            assert level > arb(".00000694")
            least_level = min(least_level,level)
        for size, seen in zip((11,12,13),sets,strict=True):
            seen.update(combinations(support,size))
    counts = [len(seen) for seen in sets]
    assert counts == [1494,297,27]
    print("不同十一、十二、十三站支撑数",counts,flush=True)
    print("三个两列判据的最小余量、四站谱严格下界",least_gain.str(14),least_level.str(14),flush=True)
    return sets


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
