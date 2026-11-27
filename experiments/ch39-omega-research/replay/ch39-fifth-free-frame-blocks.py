"""从共同三阶压缩求两组四站系数，允许共同三代平面改变。"""
from flint import acb, acb_mat, arb, arb_mat


def adj(matrix):
    return matrix.conjugate().transpose()


def solve(root, lib, blocks, name, integers, directions):
    emb = lib("common-embedding-central-right-ckm")
    _, _, intrinsic, base = emb["common_embedding"](*lib("independent-unitary-quark-right-holdout")["GEOMETRY"])
    rotations = lib("fifth-two-charge-full-quark")
    frame0 = (rotations["rotation"](rotations["QUATERNION"])*rotations["rotation"](
        lib("fourth-full-ckm-gravity")["QUATERNION"])*base*intrinsic)
    targets = [adj(frame0)*b*frame0 for b in blocks]
    identity = lib("rank-two-generator-up-quark")["eye"](4)
    make = lib("next-coherent-full-ckm")["block_from_coordinates"]
    h = acb_mat(4,4)
    seen = set()
    for line in (root/(name+"-frame.tsv")).read_text().splitlines():
        i,j,a,b = map(int,line.split())
        assert 0<=i<=j<4 and (i,j) not in seen and (i!=j or b==0)
        seen.add((i,j))
        z = acb(arb(a)/10**13,arb(b)/10**13)
        h[i,j],h[j,i]=z,z.conjugate()
    assert len(seen)==10
    frame = (2*(identity-acb(0,1)*h).inv()-identity)*frame0
    basis = [make([arb(int(i==j)) for j in range(6)]) for i in range(6)]
    coeff = [arb(v)/10**14 for v in integers]
    generators = []
    for i,j,part in directions:
        e = acb_mat(4,4)
        z = acb(0,1) if part=="i" else acb(1)
        e[i,j],e[j,i]=z,-z.conjugate()
        generators.append(e)

    def coords(b):
        return [v for i in range(3) for j in range(i+1)
                for v in ([b[i,i].real] if i==j else [b[i,j].real,b[i,j].imag])]

    def system(point):
        ts=[c+d for c,d in zip(coeff,point[:12],strict=True)]
        left=identity
        for u,e in zip(point[12:],generators,strict=True):
            left=left-u*e
        left=left.inv()
        f=(2*left-identity)*frame
        dfs=[2*left*e*left*frame for e in generators]
        residual,jac=[],[]
        for k, goal in enumerate(targets):
            b=make(ts[k*6:(k+1)*6])
            residual.extend(coords(adj(f)*b*f-goal))
            columns=[coords(adj(f)*v*f) if i//6==k else [arb(0)]*9 for i,v in enumerate(basis+basis)]
            columns.extend(coords(adj(d)*b*f+adj(f)*b*d) for d in dfs)
            jac.extend(arb_mat(columns).transpose().tolist())
        return arb_mat([[v] for v in residual]),arb_mat(jac)

    center=[arb(0)]*18
    for _ in range(3):
        f,j=system(center)
        q=j.mid().solve(f.mid()).entries()
        center=[(v-d).mid() for v,d in zip(center,q,strict=True)]
    f,j=system(center)
    pre=j.mid().inv()
    assert not pre.det().contains(0)
    offset=max(abs(v).upper() for v in (pre*f).entries())
    radius=2*offset
    box=[v+arb(0,radius.upper()) for v in center]
    _,wide=system(box)
    unit=arb_mat([[int(i==j) for j in range(18)] for i in range(18)])
    ratio=lib("next-common-ckm")["outer_sum"](unit-pre*wide)+offset/radius
    print("三代标架与四站系数的十八维自映比",ratio.upper().str(10),flush=True)
    assert ratio < 1
    ts=[c+d for c,d in zip(coeff,box[:12],strict=True)]
    return [make(ts[:6]),make(ts[6:])]
