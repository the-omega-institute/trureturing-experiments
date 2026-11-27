"""按 12、32、52、52 四轮逐边施同长共门时，三个原生第三枢轴弱差不能同时归零。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def frobenius(matrix):
    return sum(abs(z).upper()**2 for z in matrix.entries()).sqrt().upper()


def coordinates(previous):
    eye, const, forms = previous["constraints"]()
    # 三个结构零的位置不入列，两个相等的对角元共用一个坐标。
    indices = [i for i in range(21) if i not in (6, 7, 11, 15)]
    j = acb_mat([[int((10 if i == 15 else i) == k) for k in indices] for i in range(21)])
    rows = [acb_mat([[j[i, k] for k in range(len(indices))]]) for i in range(21)]

    def pair(i, k):
        return adj(rows[i]) * rows[k]

    reduced = [adj(j)*q*j for q in forms]
    const = adj(j)*const*j
    cross = pair(3, 2) + pair(10, 14)
    reduced.extend((pair(2, 2)+pair(14, 14)-pair(3, 3),
                    pair(3, 3)+pair(10, 10)-const,
                    pair(8, 8)+pair(9, 9)-pair(3, 3),
                    (cross+adj(cross))/2, (cross-adj(cross))*acb(0, -arb(1)/2)))
    return j, adj(j)*eye*j, const, reduced


def positive(matrix):
    for i in range(matrix.nrows()):
        pivot = matrix[i, i].real
        assert pivot > 0, (i, pivot)
        for a in range(i+1, matrix.nrows()):
            for b in range(a, matrix.ncols()):
                x = matrix[a,b]-matrix[a,i]*matrix[i,b]/pivot
                matrix[a,b] = x
                if a != b:
                    matrix[b,a] = x.conjugate()


def check(root, source, bases):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    previous = lib("two-block-low-line-gram")
    blocks = previous["COPIES"]
    four, two = blocks
    coordinates_used = ((0,four,0,1),(1,two,0,1),(2,four,0,2),(2,four,0,3))
    layout = tuple((gate, tuple((s[a],s[b]) for s in group))
                   for gate,group,a,b in coordinates_used)
    assert layout == lib("two-shell-ordered-bond-gates")["LAYERS"]
    j, metric, const, forms = coordinates(previous)
    lines = (root/"ch39-ordered52-weak-gram.tsv").read_text().splitlines()
    assert len(lines) == len(previous["SURVIVORS"])
    results = []
    for line, selection in zip(lines, sorted(previous["SURVIVORS"]), strict=True):
        nums = list(map(int, line.split()))
        assert tuple(nums[:3]) == selection and len(nums) == 5+len(forms)+72
        target, jitter = (arb(v)/10**14 for v in nums[3:5])
        assert jitter >= 0 and target > 7*jitter
        ps = iter(nums[5+len(forms):])
        pre = acb_mat([[acb(arb(next(ps))/10**12,arb(next(ps))/10**12)
                        for _ in range(6)] for _ in range(6)])
        # 先按三个壳的精确本征基合并，再乘带根盒宽度的源列。
        col = []
        for block in blocks:
            for a,b in product(range(len(block[0])),repeat=2):
                col.append([sum(sum(pre[row,2*k+h]*bases[k][selection[k],s[a]]
                                    for k in range(3))*source[s[b],h]
                                for s in block for h in range(2)) for row in range(6)])
        col.append([sum(sum(pre[row,2*k+h]*bases[k][selection[k],14] for k in range(3))
                        *source[14,h] for h in range(2)) for row in range(6)])
        mat = acb_mat(col).transpose()*j
        mid = mat.mid()
        gram = adj(mid)*mid + jitter*metric - target*const
        for value,q in zip(nums[5:5+len(forms)], forms, strict=True):
            gram += (arb(value)/10**14)*q
        positive(gram)
        error = frobenius(mat-mid)
        lower = ((target-7*jitter).sqrt()-arb(7).sqrt()*error)/frobenius(pre)
        print("低线组合、加权输入误差、实际两重列投影下界",
              tuple(k+1 for k in selection),error.str(12),lower.str(15),flush=True)
        assert lower > arb("4.6e-7" if selection == (9,15,15) else "8.2e-7")
        results.append(lower.lower())
    return min(results)


def verify(root):
    lib = lambda name: run_path(str(root/("ch39-"+name+".py")))
    fifth, columns, levels = lib("fifth-full-positive-joint")["verify"](root,return_operator="spectral")
    # 复用先前的十六种排除及从三壳弱差为零推出的投影上界。
    lib("two-block-low-line-gram")["experiment"](root,fifth,columns,levels)
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    bases = [adj(spectrum(shell,6,n)[1]) for n in range(2_306_888_100_092,2_306_888_100_095)]
    ctx.prec = 256
    return check(root,fifth*columns,bases)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
