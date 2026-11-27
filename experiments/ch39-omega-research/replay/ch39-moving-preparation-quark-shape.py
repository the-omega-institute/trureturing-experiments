"""在十个相位分区验证共用制备背景的两部门相对奇异值比。"""
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).resolve().parent


def lib(name):
    return run_path(str(ROOT / ("ch39-"+name+".py")))


def orthonormal_sources(graph, adjoint):
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](ROOT)
    ports = [graph["soft"]["pole"], *graph["boundary_sites"]()]
    columns = []
    for row in gram[1:]:
        vector = acb_mat([[row[ports.index(i)] if i in ports else 0] for i in range(19)])
        for column in columns:
            vector -= column*(adjoint(column)*vector)[0,0]
        columns.append(vector/(adjoint(vector)*vector)[0,0].real.sqrt())
    sources = acb_mat([v.entries() for v in columns]).transpose()
    overlap = adjoint(sources)*sources
    assert all(overlap[i,j].contains(int(i==j)) for i in range(3) for j in range(3))
    return sources


def certify(worker):
    graph = lib("three-ray-chiral-boundary")
    pump = lib("native-light-projector-pump")
    cayley = lib("first-ray-rational-windows")["cayley"]
    values = lib("first-zero-relative-quark-spectrum")["singular_values"]
    ctx.prec = 160
    adjoint = pump["adjoint"]
    sources = orthonormal_sources(graph,adjoint)
    models = []
    for charge in (-4,2):
        base,links = pump["clock_matrices"](graph,charge)
        sectors = []
        for frame in pump["rotation_bases"](graph):
            a,b = [adjoint(frame)*m*frame for m in links]
            squares = [adjoint(m)*m+m*adjoint(m) for m in (a,b)]
            sectors.append((adjoint(frame)*base*frame,a,b,*squares,adjoint(frame)*sources))
        models.append(sectors)

    @lru_cache(maxsize=8192)
    def node(n,i):
        quarter,k = divmod(i,n)
        center = cayley(arb(2*k+1)/(2*n))
        radius = max(abs(center-cayley(arb(a)/n)).upper() for a in (k,k+1))
        z = acb(0,1)**quarter*center
        # 以原钟的平方为坐标。上型取负二次幂，下型取一次幂。
        return (z**-2,2*radius),(z,radius)

    def enclose(n,i,j):
        matrices,errors = [],[]
        for (u,ru),(v,rv),sectors in zip(node(n,i),node(n,j),models,strict=True):
            rows = []
            mismatch = arb(0)
            change = arb(2).sqrt()*(ru+rv)+6*arb("1e-8")
            for base,a,b,sa,sb,inputs in sectors:
                center = base+a*u+adjoint(a*u)+b*v+adjoint(b*v)
                vector,gap = pump["first_line"](center)
                if not gap > 2*change:
                    return None
                column = acb_mat([[z] for z in vector])
                row = adjoint(column)
                residual = 6*arb("1e-8")
                for radius,square in ((ru,sa),(rv,sb)):
                    residual += radius*((row*square*column)[0,0].real.max(arb(0))).sqrt()
                # 三字符正交。逐线配相后的标架差取三个范数的最大者。
                mismatch = max(mismatch,(arb(2).sqrt()*residual/(gap-change)).upper())
                rows.append((row*inputs).entries())
            matrices.append(acb_mat(rows))
            errors.append(mismatch+arb("2e-8"))
        return matrices,errors

    def passes(data):
        if data is None:
            return False
        matrices,errors = data
        for a,b in ((0,1),(1,0)):
            denominator = values(matrices[b])[0].lower()
            if not denominator > errors[b]:
                continue
            ratios = values(matrices[a]*matrices[b].inv())
            assert ratios[0] < ratios[1] < ratios[2]
            low = (ratios[1]-errors[a]/denominator)/(1+errors[b]/denominator)
            high = (ratios[2]+errors[a]/denominator)/(1-errors[b]/denominator)
            if low > 0 and high < 32*low:
                return True
        return False

    cells = [(4,i,j) for k,(i,j) in enumerate(product(range(16),repeat=2)) if k%10==worker]
    accepted = finest = 0
    while cells:
        n,i,j = cells.pop()
        if passes(enclose(n,i,j)):
            accepted += 1
            finest = max(finest,n)
        else:
            assert n < 16384,(n,i,j)
            cells.extend((2*n,2*i+a,2*j+b) for a,b in product(range(2),repeat=2))
    return worker,accepted,finest


def target_gaps():
    ctx.prec = 192
    def box(m,r):
        return arb(m)+arb(0,arb(r).upper())
    up = [box(*x) for x in ((".00000615",".00000014"),(".00311",".00005"),(".8616",".0043"))]
    down = [box(*x) for x in ((".0000135",".0000002"),(".000268",".000003"),(".01401",".00009"))]
    angles = [box(*x) for x in ((".22704",".00082"),(".04275",".00042"),(".00377",".00008"),("1.139",".023"))]
    ckm = lib("shell-clock-ckm-blocks")["standard_matrix"]
    v = acb_mat(ckm(*(z.sin() for z in angles[:3]),angles[3]))
    target = acb_mat([[up[i]*v[i,j]/down[j] for j in range(3)] for i in range(3)])
    inverse = acb_mat([[down[i]*v[j,i].conjugate()/up[j] for j in range(3)] for i in range(3)])
    det = up[0]*up[1]*up[2]/(down[0]*down[1]*down[2])
    def frob(matrix):
        return sum(abs(x)**2 for x in matrix.entries()).sqrt().upper()
    highest_squared = sum(abs(target[2,i])**2 for i in range(3)).lower()
    inverse_squared = sum(abs(inverse[i,0])**2 for i in range(3)).lower()
    high_gap = highest_squared/(frob(inverse)*det.upper())
    low_gap = inverse_squared*det.lower()/frob(target)
    assert high_gap > arb("61.32") and low_gap > arb("77.12")
    assert (arb(9)/7)*arb("1.197e-8") < arb("2e-8")
    assert 32*(arb(9)/7)**2 < arb("52.90") < high_gap < low_gap
    print("1 TeV 参考窗上、下相邻奇异值比的下界",high_gap.str(10),low_gap.str(10))


if __name__ == "__main__":
    target_gaps()
    with ProcessPoolExecutor(max_workers=10) as pool:
        for worker,count,depth in pool.map(certify,range(10)):
            print("分区、终格数、最细象限等分数",worker,count,depth,flush=True)
