"""八站不可见方向的三维试探空间平均能量超额，中央点与全相位证书。"""
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


CHARGES = (1, -4, 2, -3, 6)
CENTRAL_BOUNDS = {1: (".47613", ".50359"), -4: (".40561", ".51082"),
                  2: (".43913", ".49063"), -3: (".43536", ".52423"), 6: (".38024", ".52557")}


def lib(root, name):
    return run_path(str(root / ("ch39-"+name+".py")))


def adj(a):
    return a.conjugate().transpose()


def eye(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])


def eig(matrix):
    values, vectors = matrix.eig(right=True, algorithm="rump")
    assert all(v.is_finite() and v.imag.contains(0) for v in values)
    order = sorted(range(len(values)), key=lambda i:values[i].real.mid())
    values = [values[i].real for i in order]
    assert all(a < b for a,b in zip(values, values[1:]))
    return values, vectors, order


def column(matrix, j):
    v = acb_mat([[matrix[i,j]] for i in range(matrix.nrows())])
    return v / (adj(v)*v)[0,0].real.sqrt()


def principal(matrix, ids):
    return acb_mat([[matrix[i,j] for j in ids] for i in ids])


def restrict(v, ids):
    return acb_mat([[v[i,0]] for i in ids])


def central(graph, outside, h, q):
    H = graph["laplacian"](q, h)
    levels, vectors, order = eig(H)
    low = levels[2]
    pos = principal(H, outside)-low*eye(len(outside))
    for j in (0,1):
        x = restrict(column(vectors,order[j]),outside)
        pos += (low-levels[j])*x*adj(x)
    energies, inner, index = eig(pos)
    xlocal = column(inner,index[0])
    x = acb_mat(19,1)
    for i, full in enumerate(outside):
        x[full,0] = xlocal[i,0]
    projector = x*adj(x)
    orth = eye(19)-projector
    restricted = orth*H*orth + 40*projector
    trial,_,_ = eig(restricted)
    assert trial[1] < 40
    upper = ((adj(x)*H*x)[0,0].real+sum(trial[:2])-sum(levels[:3]))/3
    lower = energies[0]/3
    lo,hi = CENTRAL_BOUNDS[q]
    assert arb(lo) < lower < upper < arb(hi)
    print("中央荷、平均超额下界和见证上界",q,lower.str(12),upper.str(12),flush=True)


def positive(matrix):
    # Hermitian 矩阵不选主元作 Schur 消元，每个严格正枢轴证正定。
    size = matrix.nrows()
    work = matrix.tolist()
    for i in range(size):
        a = work[i][i].real
        assert a.is_finite() and a>0
        for j in range(i+1,size):
            factor = work[j][i]/a
            for k in range(i+1,j+1):
                work[j][k] -= factor*work[k][i].conjugate()


def torus(case):
    root, outside, q, half = case
    ctx.prec = 144
    graph = lib(root, "three-ray-chiral-boundary")
    pump = lib(root, "native-light-projector-pump")
    base, clocks = pump["clock_matrices"](graph, q)
    frames = list(pump["rotation_bases"](graph))
    cells = lib(root, "terminal-port-trial-energy")["circle_cells"]()
    phases = [z for z,r in cells]
    n = len(cells)
    blocks = []
    for frame in frames:
        A, B = [adj(frame)*link*frame for link in clocks]
        blocks.append((adj(frame)*base*frame, [A*z+adj(A*z) for z in phases],
                      [B*z+adj(B*z) for z in phases]))
    restrict_frames = [acb_mat([frame.tolist()[i] for i in outside]) for frame in frames]
    reduced = principal(base,outside)
    reduced_clocks = [principal(link,outside) for link in clocks]
    clock_rows = [[a*z+adj(a*z) for z in phases] for a in reduced_clocks]
    # 两钟分别占六条无向边。每格用有理 Cayley 端点的弦长半径。
    # 三条星边的图拉普拉斯扰动 Frobenius 范数用三条秩一边范数之和界定。
    theta = arb("8.685")
    target = arb(".24")
    identity = eye(len(outside))
    for i in range(half*n//2,(half+1)*n//2):
        for j in range(n):
            drift = (12*(cells[i][1]**2+cells[j][1]**2)).sqrt()+6*arb("1e-8")
            pos = reduced+clock_rows[0][i]+clock_rows[1][j]-(theta+3*target+drift)*identity
            for (H,A,B),frame in zip(blocks,restrict_frames,strict=True):
                block = H+A[i]+B[j]
                vector,_ = pump["first_line"](block)
                v = acb_mat([[z] for z in vector])
                val = (adj(v)*block*v)[0,0].real
                assert val < theta
                w = frame*v
                pos += (theta-val)*w*adj(w)
            positive(pos)
    return q,half,n*n//2


def verify(root):
    graph = lib(root, "three-ray-chiral-boundary")
    sites = lib(root, "terminal-port-support")["verify"](root)
    outside = [i for i in range(19) if i not in sites]
    ctx.prec = 144
    h = acb.zeta_zero(2_306_888_100_085).imag
    for q in CHARGES:
        central(graph,outside,h,q)
    cases = [(root,outside,q,half) for q,half in product(CHARGES,range(2))]
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        results = list(pool.map(torus,cases))
    assert {(q,h) for q,h,count in results} == set(product(CHARGES,range(2)))
    print("五荷十个半环面的方格总数、不可见所需的平均超额严格下界", sum(x[2] for x in results), "0.24")


if __name__ == "__main__":
    source = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(source)
