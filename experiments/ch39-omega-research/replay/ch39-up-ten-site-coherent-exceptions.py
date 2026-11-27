"""用共同酉下的两列相干约束排除十站逐层覆盖所余的 32 个支撑。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


THRESHOLD = ".00000629"
BASE = (0,1,2,3,4,5,14)


def adj(a):
    return a.conjugate().transpose()


def columns(matrix, indices):
    return acb_mat([[matrix[i,j] for j in indices] for i in range(matrix.nrows())])


def test_sets():
    pairs = []
    for new, params in ((6,(1,16,".004",".000011",".000013")),
                        (9,(16,18,".012",".000013",".0000043")),
                        (10,(16,18,".012",".000013",".0000043"))):
        for duo in combinations((11,12,18),2):
            pairs.append((tuple(sorted(BASE+(new,)+duo)),params))
    pairs.extend(((tuple(sorted(BASE+(11,12,18))),(1,10,".3",".001",".003")),
                  ((0,3,4,9,11,12,14,15,16,18),(16,18,".000014",".0000076",".0000065")),
                  ((9,10,11,12,13,14,15,16,17,18),(16,18,".000017",".0000024",".0000039"))))
    restricted = [tuple(sorted(BASE+orbit)) for orbit in ((9,15,16),(10,13,17))]
    return pairs, restricted


def verify(root):
    local = run_path(str(root/"ch39-up-eight-site-cover.py"))
    model = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    positive = run_path(str(root/"ch39-up-quark-arbitrary-site-phase-floor.py"))["positive"]
    weighted = run_path(str(root/"ch39-common-unitary-quark-right-trace.py"))["weighted"]
    read_rows = local["rows"]
    n = 2_306_888_100_086
    _, labels, _ = model["verify_nodes"](root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    source, target = [spec(shell,-4,k)[1] for k in (n,n+1)]
    assert min(labels[n,-4]) > 0
    h = weighted(source,labels[n,-4])
    four = read_rows(target,(14,8,7,6))
    tau = arb(THRESHOLD)
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    assert {rotation[i] for i in (14,8,7,6)} == {14,8,7,6}
    easy, hard = test_sets()
    certificates = easy + [(s,None) for s in hard]
    covered = set()
    for support, params in certificates:
        inside = read_rows(target,support)
        outside = acb_mat(target)
        for site in support:
            for j in range(19):
                outside[site,j] = 0
        mixed = h*outside
        rhs = read_rows(mixed,support)
        block = acb_mat([[h[i,j] for j in support] for i in support])
        gram, constant = adj(inside)*inside, adj(outside)*mixed
        norm = [gram[j,j].real for j in range(19)]
        gammas = local["multipliers"](block,rhs,norm)

        def schur(ids, gamma, coordinates):
            matrix = acb_mat([[block[i,j] for j in coordinates] for i in coordinates])
            for i in range(len(coordinates)):
                matrix[i,i] += gamma
            positive(matrix)
            data = columns(read_rows(rhs,coordinates),ids)
            main = acb_mat([[constant[i,j]-gamma*gram[i,j] for j in ids] for i in ids])
            return main-adj(data)*matrix.solve(data)

        all_sites = list(range(len(support)))
        floors = [max(arb(0),schur([j],g,all_sites)[0,0].real.lower()) for j,g in enumerate(gammas)]
        def alternatives(prior, row):
            j,k = row[0]-1,row[1]-1
            gamma, r, s = map(arb,row[2:])
            form = schur((j,k),gamma,all_sites)
            a, b = r-form[0,0].real, s-form[1,1].real
            assert a > 0 and b > 0
            slack = abs(form[0,1])-(r*s).sqrt()-(a*b).sqrt()
            assert slack > 0
            print(support,"两列阈值矛盾的余量",slack.lower().str(10),flush=True)
            result = []
            for index, cap in ((j,r),(k,s)):
                raised = list(prior)
                raised[index] = max(raised[index],cap)
                result.append(raised)
            return result

        if params is not None:
            branches = alternatives(floors,params)
            if params[:2] == (1,10):
                left, right = branches
                branches = (alternatives(left,(10,16,".03",".00057",".000018"))
                            + alternatives(right,(16,18,".0095",".00002",".000013")))
            for raised in branches:
                matrix = weighted(four,raised)
                for i in range(4):
                    matrix[i,i] -= tau
                positive(matrix)
        else:
            # 用四站最小本征线限制候选输出态在第 16、18 层上的读出。
            matrix = weighted(four,floors)
            values,vectors = matrix.eig(right=True,algorithm="rump")
            order = sorted(range(4),key=lambda j: values[j].real.mid())
            levels = [values[j].real for j in order]
            assert all(a < b for a,b in zip(levels,levels[1:]))
            ratio = (tau-levels[0])/(levels[1]-levels[0])
            assert 0 < ratio < 1
            y = columns(vectors,[order[0]])
            y *= 1/(adj(y)*y)[0,0].real.sqrt()
            ids = (15,17)
            readout = [arb(".32"),arb(".51")]
            for j,u in zip(ids,readout,strict=True):
                v = columns(four,[j])
                size = abs((adj(y)*v)[0,0])
                perp = ((adj(v)*v)[0,0].real-size**2).sqrt()
                projection = (1-ratio).sqrt()*size-ratio.sqrt()*perp
                assert projection > u.sqrt()

            # 存在单位二列组合，其活动部分经共同酉映到站 14 的正交补。
            keep = [i for i,site in enumerate(support) if site != 14]
            form = schur(ids,arb("-.0003"),keep)
            bound = tau-levels[0]+sum(u*floors[j] for j,u in zip(ids,readout,strict=True))
            scaled = acb_mat([[form[i,j]*(readout[i]*readout[j]).sqrt()
                               for j in range(2)] for i in range(2)])
            for i in range(2):
                scaled[i,i] -= bound
            positive(scaled)
            gaps = [z.real for z in scaled.eig(algorithm="rump")]
            print(support,"消去极点后的加权振幅余量",min(gaps).lower().str(10),flush=True)
        current = support
        for _ in range(3):
            covered.add(current)
            current = tuple(sorted(rotation[i] for i in current))
        assert current == support
    assert len(certificates) == 14 and len(covered) == 32
    print("共同酉证书覆盖的十站支撑",len(covered),flush=True)
    return covered


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
