"""第四壳所有独立换基协变的第三步共同生成元的下型质量整盒下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def positive(matrix):
    # 中点 Cholesky 仅选取可逆合同变换，最后对整个区间矩阵验证严格正对角占优。
    n = matrix.nrows()
    mid = matrix.mid()
    lower = acb_mat(n,n)
    for i in range(n):
        for j in range(i+1):
            residual = mid[i,j]-sum(lower[i,k]*lower[j,k].conjugate() for k in range(j))
            if i==j:
                assert residual.real>0
                lower[i,j] = residual.real.sqrt().mid()
            else:
                lower[i,j] = (residual/lower[j,j]).mid()
    q = lower.inv()
    checked = q*matrix*adj(q)
    margin = min((checked[i,i].real-sum(abs(checked[i,j]) for j in range(n) if i!=j)).lower()
                 for i in range(n))
    assert margin>0,margin
    return margin


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    x,y,*_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec=384
    eye=rank["eye"](19)
    spectrum=lib("adjacent-pvm-transport")["spectrum"]
    shell=lib("three-ray-chiral-boundary")
    source,target=[spectrum(shell,2,node)[1] for node in (n+1,n+4)]
    eig,axes=(y*adj(x)).eig(right=True,algorithm="rump")
    angles=[z.log().imag for z in eig]
    order=sorted(range(19),key=lambda j:angles[j].mid())
    assert all(angles[i]<angles[j] for i,j in zip(order,order[1:]))
    gap=min(abs(a-b).lower() for i,a in enumerate(eig) for b in eig[i+1:])
    assert gap>arb(".04")
    sizes=[sum(abs(axes[i,j])**2 for i in range(19)).sqrt() for j in range(19)]
    beam=acb_mat([[axes[i,j]/sizes[j] for j in order] for i in range(19)])
    assert (adj(beam)*beam).contains(eye)
    w=labels[n+1,2]
    z=adj(beam)*y*y*x*source
    density=acb_mat([[z[i,j]*w[j] for j in range(19)] for i in range(19)])*adj(z)
    q=adj(beam)*target
    kernels=[acb_mat([[density[k,l]*q[k,j].conjugate()*q[l,j] for l in range(19)] for k in range(19)])
             for j in range(19)]
    r3=arb(3).sqrt()
    omega=acb(-1,r3)/2
    rays=[[target[14,j],sum(target[i,j] for i in (8,7,6))/r3] for j in range(19)]
    chars=[[abs(target[8,j]+phase*target[7,j]+phase**2*target[6,j])**2/3 for j in range(19)]
           for phase in (omega.conjugate(),omega)]
    def weighted(coefficients):
        return sum((K*p for K,p in zip(kernels,coefficients,strict=True)),acb_mat(19,19))
    pieces=[[weighted([r[a]*r[b].conjugate() for r in rays]) for b in range(2)] for a in range(2)]
    tensor=acb_mat([[pieces[i//19][j//19][i%19,j%19] for j in range(38)] for i in range(38)])
    forms={"tensor":tensor,"average":pieces[1][1],"char1":weighted(chars[0]),"char2":weighted(chars[1])}
    cert={}
    for row in (root/"ch39-fourth-equivariant-dual.tsv").read_text().splitlines():
        key,*entries=row.split()
        assert key in forms and key not in cert and len(entries)==19
        cert[key]=[arb(int(v))/10**12 for v in entries]
    assert set(cert)==set(forms)
    for key,form in forms.items():
        dual=cert[key]
        shifted=acb_mat(form)
        for i in range(form.nrows()):
            shifted[i,i]-=dual[i%19]
        margin=positive(shifted)
        print("证书、对偶下界、合同后对角占优余量",key,sum(dual).str(15),margin.str(10),flush=True)
    floors={key:sum(vals) for key,vals in cert.items()}
    light=min(floors[s].lower() for s in ("tensor","char1","char2"))
    strange=min(floors[s].lower() for s in ("average","char1","char2"))
    assert light>arb(".00004232") and strange>arb(".00141168")
    print("单位圆最小谱隙、任意三平面轻下与奇质量下界",gap.str(12),light.str(14),strange.str(14),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
