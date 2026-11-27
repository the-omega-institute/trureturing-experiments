"""右外壳补空间三次单位根下的轻上夸克见证与局部荷六反求。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


WITNESSES = [((14,17,8,16,4,13), tuple(j%3 for j in range(13)))]


def orthonormal(rows, dim):
    matrix = acb_mat(dim,dim)
    for i in range(dim):
        ints = next(rows)
        assert len(ints) == 2*dim
        for j in range(dim):
            matrix[i,j] = acb(ints[2*j],ints[2*j+1])
    vectors = []
    for j in range(dim):
        v = acb_mat([[matrix[i,j]] for i in range(dim)])
        for basis in vectors:
            v -= basis*(adj(basis)*v)[0,0]
        norm = (adj(v)*v)[0,0].real.sqrt()
        assert norm > 0
        vectors.append(v*(1/norm))
    return acb_mat([[vectors[j][i,0] for j in range(dim)] for i in range(dim)])


def adj(a):
    return a.conjugate().transpose()


def previous_three_columns(root, target):
    original = run_path(str(root/"ch39-up-gravity-block-stabilizers.py"))
    direct = original["build_inputs"](root)[0]
    rest = [j for j in range(19) if j not in (14,6,7)]
    v = acb_mat([[target[i,j] for j in (15,17)] for i in rest])
    unitary = acb_mat([[direct[i,j] for j in rest] for i in rest])
    a, b = adj(v)*v, adj(v)*unitary*v
    roots = [acb(1),acb(-arb(1)/2,arb(3).sqrt()/2),acb(-arb(1)/2,-arb(3).sqrt()/2)]
    largest = arb(-1)
    # 相差 2π/3 的整体相位只会重排三个 Gram 矩阵。
    for i in range(1024):
        angle = (arb(i)+arb(1)/2+arb(0,arb(".5")))*2*arb.pi()/(3*1024)
        sin,cos = angle.sin_cos()
        rotated = b*acb(cos,sin)
        smallest = []
        for phase in roots:
            part = rotated*phase
            matrix = (a+part+adj(part))*(arb(1)/3)
            spread = (((matrix[0,0]-matrix[1,1]).real/2)**2+abs(matrix[0,1])**2).sqrt()
            smallest.append(((matrix[0,0]+matrix[1,1]).real/2-spread).upper())
        largest = max(largest,min(smallest))
    assert largest < arb("-.029")
    print("三站旧见证两列在全部相位上的最小 Gram 本征值上界",largest.str(13),flush=True)


def gravity_choices(local_adj, sites, charge_source, charge_target):
    # 局部生成元满足 X²=U*。每条单纯本征线上有两个平方根，因而有两个 G=X³。
    values, vectors = local_adj.eig(right=True,algorithm="rump")
    assert all(not z.overlaps(w) for j,z in enumerate(values) for w in values[j+1:])
    assert all(not z.imag.contains(0) and abs(z).contains(1) for z in values)
    parts = []
    for j,z in enumerate(values):
        v = acb_mat([[vectors[i,j]] for i in range(len(sites))])
        projector = v*adj(v)*(1/(adj(v)*v)[0,0].real)
        parts.append(projector*z*z.sqrt())
    inside = acb_mat([[charge_source[i,k] for k in (17,18)] for i in sites])
    measure = acb_mat([[charge_target[i,15].conjugate() for i in sites]])
    full = acb_mat([[(adj(charge_target)*charge_source)[15,k] for k in (17,18)]])
    tail = full-measure*inside
    least = arb(10)
    for signs in product((-1,1),repeat=len(sites)):
        action = sum((s*p for s,p in zip(signs,parts,strict=True)),acb_mat(len(sites),len(sites)))
        local = measure*action*inside
        amplitude = (local*adj(local)+tail*adj(tail))[0,0].real-2*abs((local*adj(tail))[0,0])
        least = min(least,amplitude.lower())
    # 同时容许一个局部相对于块外的相位，因此该界包括整体重相的上型见证。
    assert least > arb(".000001")
    return least


def verify(root):
    up = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    common = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spectra = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    weighted = run_path(str(root/"ch39-common-unitary-quark-right-trace.py"))["weighted"]
    cutoffs = run_path(str(root/"ch39-adjacent-ordered-charge-six-cutoffs.py"))["verify"]
    n = 2_306_888_100_086
    _,weights,_ = common["verify_nodes"](root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    _,positions = cutoffs(root)[2]
    high = min(z.sqrt().lower() for z in positions[17:])
    assert high > arb(".1")
    ctx.prec = 640
    source,target = [spectra(shell,-4,j)[1] for j in (n,n+1)]
    previous_three_columns(root,target)
    gsource,gtarget = [spectra(shell,6,j)[1] for j in (n,n+1)]
    h = weighted(source,weights[n,-4]); eye=up["eye"]
    four = acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    data = (root/"ch39-up-order-three-complement-roots.tsv").read_text().splitlines()
    rows = iter(tuple(map(int,s.split())) for s in data if s and not s.startswith("#"))
    for group,pattern in WITNESSES:
        rest = [i for i in range(19) if i not in group]
        local = orthonormal(rows,len(group))
        basis = orthonormal(rows,len(rest))
        free = weighted(basis,[omega**k for k in pattern])
        assert (free*free*free).contains(eye(len(rest)))
        transport_adj = acb_mat(19,19)
        for idx,block in ((group,local),(rest,free)):
            for a,i in enumerate(idx):
                for b,j in enumerate(idx):
                    transport_adj[i,j]=block[a,b]
        assert (adj(transport_adj)*transport_adj).contains(eye(19))
        image = transport_adj*target
        energy = adj(image)*h*image
        raw = weighted(four,[energy[j,j].real for j in range(19)])
        eigen = sorted((v.real for v in raw.eig(algorithm="rump")),key=lambda x:x.mid())
        assert eigen[0] < arb(".00000601") < arb(".00000629") < eigen[1]
        assert all(a<b for a,b in zip(eigen,eigen[1:]))
        gap = gravity_choices(local,group,gsource,gtarget)
        ceiling = 16*arb("91.1876")/arb("1.22089e19")
        assert high*gap > ceiling
        print("站位、四站升序谱、荷六第16层振幅下界",group,[x.str(13) for x in eigen],
              (high*gap).lower().str(13),flush=True)
    assert next(rows,None) is None


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
