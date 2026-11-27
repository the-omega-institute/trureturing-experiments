"""长度 32 使用原生带相位边的两个谱模；任意循环复制的短边四站酉块仍排除九式。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def identity(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    rank = lib("three-shell-rank-one-obstruction")
    rank["experiment"](root, fifth, columns, levels)
    shared = lib("shared-block-ninth-gravity")
    old = lib("two-block-low-line-gram")
    graph = lib("three-ray-chiral-boundary")
    sites = graph["geo"]["sites"]
    rho = [sites.index((b, c, a)) for a, b, c in sites]
    copies = (old["COPIES"][0], ((0, 2), (3, 1), (4, 5)))
    assert all(tuple(rho[i] for i in group[j]) == group[(j+1)%3]
               for group in copies for j in range(3))
    numbers = range(2_306_888_100_092, 2_306_888_100_095)
    bases = [adj(lib("adjacent-pvm-transport")["spectrum"](graph, 6, n)[1]) for n in numbers]
    ctx.prec = 256
    # 对这三条定向边，原生边相位由球面面积给出。模型的零点相位仅作用于 724、804。
    geo = graph["geo"]
    t = geo["unit"](0, 2)**6
    for a, b in copies[1]:
        i, j = sorted((a, b))
        d, u = next((d, u) for k, l, d, u in geo["edges"] if (k, l) == (i, j))
        assert d == 32
        assert (u**6 if a == i else u.conjugate()**6).overlaps(t)
    # 循环置换保持内积与行列式，故三条有向边的相位在精确算术中相同。
    swap = acb_mat([[0, t], [t.conjugate(), 0]])
    lift = acb_mat(21, 19)
    for i in range(16):
        lift[i, i] = 1
    for j, sig in enumerate((1, -1)):
        proj = (identity(2) + sig*swap)/2
        for i, val in enumerate(proj.entries()):
            lift[16+i, 16+j] = val
    lift[20, 18] = 1
    eye = identity(19)
    const = acb_mat(19, 19)
    const[18, 18] = 1
    forms = [adj(lift)*q*lift for q in old["constraints"]()[2][:32]]
    for idx in (16, 17):
        f = -const
        f[idx, idx] = 1
        forms.append(f)
    src = fifth*columns
    def rows(selections):
        matrix = []
        for basis, k in selections:
            for j in range(2):
                terms = [sum(basis[k, copy[a]]*src[copy[b], j] for copy in group)
                         for group in copies for a, b in product(range(len(group[0])), repeat=2)]
                matrix.append([*terms, basis[k, 14]*src[14, j]])
        return acb_mat(matrix)*lift

    def spread(m):
        return sum(abs(z).upper()**2 for z in m.entries()).sqrt().upper()

    def validate(form, factors, constraints):
        for num, q in zip(factors, constraints, strict=True):
            form += (arb(num)/10**12)*q
        mid = form.mid()
        for i in range(mid.nrows()):
            mid[i, i] = mid[i, i].real
            for j in range(i):
                mid[i, j] = mid[j, i].conjugate()
        shared["positive_definite"](mid - spread(form-mid)*eye)

    records = iter((root/"ch39-cyclic-short-magnetic-gravity.tsv").read_text().splitlines())
    candidates = list(product(*rank["CANDIDATES"]))
    rest = old["SURVIVORS"]
    cap = arb("2.853e-8")
    for chosen in (x for x in candidates if x not in rest):
        nums = list(map(int, next(records).split()))
        assert tuple(nums[:3]) == chosen
        h = rows(zip(bases, chosen, strict=True))
        hm = h.mid()
        jitter = arb("1e-8")
        threshold = arb(".001")
        form = adj(hm)*hm - threshold**2*const + jitter*(eye-7*const)
        validate(form, nums[3:], forms)
        assert threshold - arb(7).sqrt()*spread(h-hm) > arb(".00099") > cap
    regs = arb("1e-5")
    delta = arb(".028")
    floor, heavy, top = levels
    for chosen in sorted(rest):
        h = rows(zip(bases, chosen, strict=True))
        hm = h.mid()
        projected = adj(hm)*(hm*adj(hm)+regs**2*identity(6)).inv()*hm
        cap_here = (cap + arb(7).sqrt()*spread(h-hm))/regs
        assert cap_here < delta
        weak = projected-delta**2*const
        for idx, sign, lim in shared["LIMITS"]:
            nums = list(map(int, next(records).split()))
            assert nums[:5] == [*chosen, idx, sign] and nums[-1] >= 0
            v = rows([(bases[-1], idx)])
            amp = floor*const + adj(v)*acb_mat([[heavy-floor, 0], [0, top-floor]])*v
            form = sign*(amp-arb(lim)*const) + arb("1e-5")*(eye-7*const)
            validate(form, nums[5:], [*forms, weak])
        print("原生边相位、低线、正则化投影范数上界", t.str(10), chosen, cap_here.str(13), flush=True)
    assert next(records, None) is None
    pivots = next(row[2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, numbers[-1]) if row[0] == 6)
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [numbers[-1]])
    shared["gravity"]([v.real for v in pivots[2]], field, heat, top)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *src["verify"](root, return_operator="spectral"))
