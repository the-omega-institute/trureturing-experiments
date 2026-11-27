"""同块酉输运若在三个壳都满足弱条件，第九壳的两个引力目标便不相容。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


# 第一条线取下界，其余取上界。振幅以普朗克能标为单位，谱线编号从零计。
LIMITS = ((0, 1, ".19"), (9, -1, ".025"), (12, -1, ".13"),
          (16, -1, ".05"), (17, -1, ".005"))


def adj(matrix):
    return matrix.conjugate().transpose()


def positive_definite(matrix):
    for k in range(matrix.nrows()):
        pivot = matrix[k, k].real
        assert pivot > 0, (k, pivot)
        for i in range(k + 1, matrix.nrows()):
            for j in range(i, matrix.nrows()):
                value = matrix[i, j] - matrix[i, k]*matrix[k, j]/pivot
                matrix[i, j] = value
                if i != j:
                    matrix[j, i] = value.conjugate()


def verify_amplitudes(root, previous, fifth, columns, levels, bases):
    copies = previous["COPIES"]
    eye, const, forms = previous["constraints"]()
    source = fifth*columns
    floor, heavy, top = levels
    delta = arb(".065")
    input_cap = arb("2.853e-8")
    jitter = arb("1e-5")
    coeffs = iter((root/"ch39-shared-block-ninth-gravity.tsv").read_text().splitlines())

    def rows(pairs):
        result = []
        for basis, line in pairs:
            for j in range(2):
                row = [sum(basis[line, s[a]]*source[s[b], j] for s in block)
                       for block in copies for a, b in product(range(len(block[0])), repeat=2)]
                result.append(row + [basis[line, 14]*source[14, j]])
        return acb_mat(result)

    for selection in sorted(previous["SURVIVORS"]):
        h = rows(zip(bases, selection, strict=True))
        hm = h.mid()
        inverse = (hm*adj(hm)).inv()
        error = sum(abs(x).upper()**2 for x in (h-hm).entries()).sqrt()
        inverse_norm = max(sum(abs(x).upper() for x in row) for row in inverse.tolist())
        # 对齐次列 z，||z||²=7，||Hz||<input_cap。下式给出精确中点行空间投影的平方界。
        squared_cap = inverse_norm*(input_cap + arb(7).sqrt()*error)**2
        assert squared_cap < delta**2, (selection, squared_cap)
        weak_form = adj(hm)*inverse*hm - delta**2*const

        for i, sign, bound in LIMITS:
            nums = [int(t) for t in next(coeffs).split()]
            assert nums[:5] == [*selection, i, sign]
            multipliers = [arb(v)/10**12 for v in nums[5:]]
            assert multipliers[-1] >= 0
            v = rows([(bases[2], i)])
            diag = acb_mat([[heavy-floor, 0], [0, top-floor]])
            amplitude = floor*const + adj(v)*diag*v
            matrix = sign*(amplitude-arb(bound)*const) + jitter*(eye-7*const)
            for weight, form in zip(multipliers, [*forms, weak_form], strict=True):
                matrix += weight*form
            # 在每个可行 z 上，酉约束和正则项为零，带非负乘子的弱约束非正。
            # 用区间矩阵与显式 Hermitian 中点的 Frobenius 距离扣除输入误差。
            mid = matrix.mid()
            for a in range(mid.nrows()):
                mid[a, a] = mid[a, a].real
                for b in range(a):
                    mid[a, b] = mid[b, a].conjugate()
            radius = sum(abs(z).upper()**2 for z in (matrix-mid).entries()).sqrt().upper()
            positive_definite(mid - radius*eye)
        print("三线编号、预条件投影范数上界、通过的振幅界数",
              tuple(k+1 for k in selection), squared_cap.sqrt().upper().str(12), len(LIMITS), flush=True)
    assert next(coeffs, None) is None


def gravity(weights, field, heat, top):
    positive = [i for i, w in enumerate(weights) if w > 0]
    negative = [i for i, w in enumerate(weights) if w < 0]
    assert positive == [0, 9, 12, 15] and len(negative) == 15
    total = sum(map(abs, weights))
    x = top.upper()**2
    pi = arb.pi()
    scale_k = -sum(a for _, a, _ in field)/(2*pi)
    scale_w = -sum(g for g, _, _ in field)/(8*pi)
    assert scale_k > 0 and scale_w > 0
    # 积分中 x^j(1-exp(-t m/x)) 随 x 递增（j>=1），可统一在最大平方截止取界。
    err_k = total*x*sum(abs(a)*(1-heat(2, m/x)) for _, a, m in field)/(2*pi)
    err_w = total*x**2*sum(abs(g)*(arb(1)/2-heat(3, m/x)) for g, _, m in field)/(4*pi)
    c = -sum(weights[i] for i in negative if i not in (16, 17))
    limits = {i: arb(b) for i, _, b in LIMITS}
    limits[15] = arb(2241)/arb("1.22089e19")
    # 假定 W>=0。扣去重的两条负权线后，对其余负权应用带权 Cauchy–Schwarz。
    d = sum(weights[i]*limits[i]**4 for i in (9, 12, 15)) + err_w/scale_w
    b = -sum(weights[i]*limits[i]**2 for i in (16, 17))
    w = weights[0]
    assert w > c > 0 and d > 0
    # 对 t=a_1²，sqrt(c(w t²+d))-w t 的导数小于 sqrt(cw)-w<0。
    t = limits[0]**2
    upper = scale_k*(b+(c*(w*t**2+d)).sqrt()-w*t)+err_k
    print("第九壳两热核余项及真空项非负时的牛顿上界",
          *[v.upper().str(14) for v in (err_k, err_w, upper)], flush=True)
    assert err_k < arb(".138") and err_w < arb(".022") and upper < arb("-2.5")


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    previous = lib("two-block-low-line-gram")
    previous["experiment"](root, fifth, columns, levels)
    ctx.prec = 320
    graph = lib("three-ray-chiral-boundary")
    numbers = range(2_306_888_100_092, 2_306_888_100_095)
    bases = [adj(lib("adjacent-pvm-transport")["spectrum"](graph, 6, n)[1]) for n in numbers]
    verify_amplitudes(root, previous, fifth, columns, levels, bases)
    ninth = numbers[-1]
    pivots = next(row[2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, ninth) if row[0] == 6)
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [ninth])
    gravity([v.real for v in pivots[2]], field, heat, levels[-1])


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
