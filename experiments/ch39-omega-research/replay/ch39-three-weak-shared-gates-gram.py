"""同长度同门时，平方弦长32的六条不交边无法实现三壳弱差同时为零。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def constraints():
    eye = acb_mat([[int(i==j) for j in range(9)] for i in range(9)])
    rows = [acb_mat([[eye[i,j] for j in range(9)]]) for i in range(9)]
    constant = adj(rows[8])*rows[8]
    forms = []
    for offset in (0,4):
        for transposed in (False,True):
            def term(i,j):
                return rows[offset+(2*j+i if transposed else 2*i+j)]
            for i,j in ((0,0),(1,1),(0,1)):
                block = adj(term(i,0))*term(j,0)+adj(term(i,1))*term(j,1)
                if i==j:
                    forms.append(block-constant)
                else:
                    forms.extend(((block+adj(block))/2,(block-adj(block))*acb(0,-arb(1)/2)))
    return eye,constant,forms


def experiment(root,fifth,columns,levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    rank = lib("three-shell-rank-one-obstruction")
    _,cap = rank["experiment"](root,fifth,columns,levels)
    edges = lib("three-weak-shortest-edge-radius")["EDGES"]
    graph = lib("three-ray-chiral-boundary")
    bases = [adj(lib("adjacent-pvm-transport")["spectrum"](graph,6,n)[1])
             for n in range(2_306_888_100_092,2_306_888_100_095)]
    ctx.prec = 256
    c = fifth*columns
    fixed = sorted(set(range(19))-{k for edge in edges for k in edge})
    eye,const,forms = constraints()
    target = arb(".00015")
    jitter = arb(".000002")
    for raw,ix in zip((root/"ch39-three-weak-shared-gates-gram.tsv").read_text().splitlines(),
                      product(*rank["CANDIDATES"]),strict=True):
        nums = list(map(int,raw.split()))
        assert tuple(nums[:3]) == ix and len(nums)-3 == len(forms)
        linear = []
        for base,k in zip(bases,ix,strict=True):
            for j in range(2):
                row = [sum(base[k,e[a]]*c[e[b],j] for e in edges[g:g+3])
                       for g in (0,3) for a,b in product(range(2),repeat=2)]
                linear.append(row+[sum(base[k,i]*c[i,j] for i in fixed)])
        h = acb_mat(linear)
        gram = adj(h)*h+eye*jitter-const*target
        for m,q in zip(nums[3:],forms,strict=True):
            gram += q*(arb(m)/10**8)
        # 先扣除 Hermitian 误差的行和范数，避免消元重复使用原输入误差。
        mid = (gram.mid()+adj(gram.mid()))/2
        error = max(sum(abs(z).upper() for z in row) for row in (gram-mid).tolist())
        gram = mid-eye*error.upper()
        # 对每个严格正的对角枢轴作 Schur 补，最后一维也要正。
        for i in range(9):
            pivot = gram[i,i].real
            assert pivot > 0,(ix,i,pivot)
            for a in range(i+1,9):
                for b in range(a,9):
                    updated = gram[a,b]-gram[a,i]*gram[i,b]/pivot
                    gram[a,b]=updated
                    if a!=b:
                        gram[b,a]=updated.conjugate()
    # 两个二阶酉块各贡献 Frobenius 范数平方 2，齐次坐标为1。
    bound = (target-5*jitter).sqrt()
    print("十八种低线组合的共同重投影范数下界与三弱必要上界",bound.str(14),cap.str(14),flush=True)
    assert bound > arb(".01183") > cap


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*src["verify"](root,return_operator="spectral"))
