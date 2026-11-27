"""全相干时，允许任意重定相的对易第二步仍排除轻上夸克。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def positive(h):
    for k in range(1,h.nrows()+1):
        minor = acb_mat([[h[i,j] for j in range(k)] for i in range(k)])
        assert minor.det().real>0,(k,minor.det())


def verify(root):
    get = lambda name: run_path(str(root/name))
    full = get("ch39-one-generator-gravity-full-ckm.py")
    center,radius = full["verify"](root)
    rank = get("ch39-rank-two-generator-up-quark.py")
    n = rank["N"]
    emb = get("ch39-common-embedding-central-right-ckm.py")
    _,labels,_ = emb["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    w = labels[n+1,-4]
    max_weight = max(v.upper() for v in w)
    shell = get("ch39-three-ray-chiral-boundary.py")
    spec = get("ch39-adjacent-pvm-transport.py")["spectrum"]
    source,target = [spec(shell,-4,node)[1] for node in (n+1,n+3)]
    anchor = get("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](
        root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    eig,vecs = anchor.eig(right=True,algorithm="rump")
    order = sorted(range(19),key=lambda j:eig[j].arg().mid())
    assert all(eig[a].arg()<eig[b].arg() for a,b in zip(order,order[1:]))
    norm = [sum(abs(vecs[i,j])**2 for i in range(19)).sqrt() for j in order]
    q = acb_mat([[vecs[i,j]/h for j,h in zip(order,norm,strict=True)] for i in range(19)])
    gap = min(abs(a-b).lower() for i,a in enumerate(eig) for b in eig[i+1:])
    dx = 2*sum(abs(z+arb(0,radius.upper())).upper()**2*(1 if i==j else 2)
               for z,(i,j,_) in zip(center,full["DIRECTIONS"],strict=True)).sqrt()
    assert gap>2*dx
    op_error = 2*arb(38).sqrt()*dx/(gap-dx)
    l,m = adj(target)*q,adj(q)*source
    duals = [[arb(int(cell))/10**12 for cell in s.split()] for s in
             (root/"ch39-next-shell-commuting-phase-duals.tsv").read_text().splitlines()]
    assert len(duals)==19 and all(len(row)==19 for row in duals)
    floors = []
    for j,row in enumerate(duals):
        k = acb_mat([[l[j,h]*m[h,i] for i in range(19)] for h in range(19)])
        gram = acb_mat([[k[h,i]*w[i] for i in range(19)] for h in range(19)])*adj(k)
        for i,d in enumerate(row):
            gram[i,i] -= d
        positive(gram)
        floors.append(max(min(z.lower() for z in w),sum(row).lower()))
    f = acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
    four = acb_mat([[f[i,j]*floors[j] for j in range(19)] for i in range(4)])*adj(f)
    raw = arb(".0000250")
    positive(four-rank["eye"](4)*raw)
    low_root = raw.sqrt()-op_error*max_weight.sqrt()
    assert low_root>arb(".00000629").sqrt()
    low = low_root**2
    distance = (low_root-arb(".00000629").sqrt())/(2*max_weight.sqrt())
    commutator = distance*(gap-2*dx)/arb(2).sqrt()
    print("简单谱隙、全盒相位算子误差、轻上界的下端",gap.str(12),op_error.upper().str(12),low.lower().str(12),flush=True)
    print("重选第二步生成元到对易族的距离及交换子 Frobenius 范数下界",
          distance.lower().str(12),commutator.lower().str(12),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
