"""三个四站块共用同一酉矩阵时，三壳弱差归零只容许两组低线选择。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

COPIES = (((6,11,9,13),(7,12,15,10),(8,18,16,17)), ((0,2),(1,3),(4,5)))
SURVIVORS = {(9,15,15),(15,9,15)}


def adj(a):
    return a.conjugate().transpose()


def constraints():
    # 每个 Hermitian 型的变量是两个方阵依行展开的条目，再附常数 1。
    n = 21
    eye = acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    rows = [acb_mat([[eye[i,j] for j in range(n)]]) for i in range(n)]
    const = adj(rows[-1])*rows[-1]
    offset = 0
    forms = []
    for size in (4,2):
        for transposed in (False,True):
            def r(i,j):
                return rows[offset+(size*j+i if transposed else size*i+j)]
            for i in range(size):
                for j in range(i,size):
                    g = sum((adj(r(i,k))*r(j,k) for k in range(size)),acb_mat(n,n))
                    if i==j:
                        forms.append(g-const)
                    else:
                        forms.extend(((g+adj(g))/2,(g-adj(g))*acb(0,-arb(1)/2)))
        offset += size*size
    return eye,const,forms


def experiment(root,fifth,columns,levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    rank = lib("three-shell-rank-one-obstruction")
    _,cap = rank["experiment"](root,fifth,columns,levels)
    graph = lib("three-ray-chiral-boundary")
    layers = lib("two-shell-ordered-bond-gates")["LAYERS"] + lib("three-shell-gravity-radius56")["EXTRA_LAYERS"]
    assert sorted(i for block in COPIES for copy in block for i in copy) == [i for i in range(19) if i!=14]
    for _,pairs in layers:
        assert any(tuple((copy[a],copy[b]) for copy in block)==pairs
                   for block in COPIES for a,b in product(range(len(block[0])),repeat=2))
    bases = [adj(lib("adjacent-pvm-transport")["spectrum"](graph,6,n)[1])
             for n in range(2_306_888_100_092,2_306_888_100_095)]
    ctx.prec = 256
    source = fifth*columns
    eye,const,forms = constraints()
    filtered = [ix for ix in product(*rank["CANDIDATES"]) if ix not in SURVIVORS]
    assert len(filtered)==16
    lower = []
    for line,selection in zip((root/"ch39-two-block-low-line-gram.tsv").read_text().splitlines(),filtered,strict=True):
        nums = list(map(int,line.split()))
        assert tuple(nums[:3]) == selection and len(nums)-5 == len(forms)
        target,jitter = (arb(v)/10**12 for v in nums[3:5])
        rows = []
        for base,k in zip(bases,selection,strict=True):
            for j in range(2):
                coeffs = [sum(base[k,s[a]]*source[s[b],j] for s in block)
                          for block in COPIES for a,b in product(range(len(block[0])),repeat=2)]
                rows.append(coeffs+[base[k,14]*source[14,j]])
        h = acb_mat(rows)
        hm = h.mid()
        error = sum(abs(x).upper()**2 for x in (h-hm).entries()).sqrt()
        gram = adj(hm)*hm+jitter*eye-target*const
        for num,q in zip(nums[5:],forms,strict=True):
            gram += q*(arb(num)/10**12)
        for i in range(21):
            pivot = gram[i,i].real
            assert pivot > 0,(selection,i,pivot)
            for a in range(i+1,21):
                for b in range(a,21):
                    value = gram[a,b]-gram[a,i]*gram[i,b]/pivot
                    gram[a,b] = value
                    if a!=b:
                        gram[b,a] = value.conjugate()
        # 酉矩阵的 Frobenius 范数平方是阶数；带齐次坐标的列范数平方是 4+2+1。
        lower.append(((target-7*jitter).sqrt()-arb(7).sqrt()*error).lower())
    print("十六组的实际投影范数一致下界、三弱条件给出的上界",min(lower).str(15),cap.str(14),flush=True)
    assert min(lower) > arb(".00050") > cap


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*src["verify"](root,return_operator="spectral"))
