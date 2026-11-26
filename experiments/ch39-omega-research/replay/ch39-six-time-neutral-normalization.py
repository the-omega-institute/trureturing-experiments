"""由六项算术时间夹定中性壳平均响应，并限制单枚零点的算术和味概率。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx

neutral = run_path(str(Path(__file__).with_name("ch39-neutral-clock-response.py")))


def compact_moments(times):
    y = arb_poly([0, 1])
    factors = [1+(k*k-1)*y for k in range(2, 7)]
    denominator = prod(factors)
    functions = [denominator, (1-y)*denominator]
    functions.extend((1-y)*prod(factors[:i]+factors[i+1:]) for i in range(len(factors)))
    matrix = arb_mat([[f[k] for k in range(7)] for f in functions])
    return denominator, matrix.solve(arb_mat([[t] for t in times])).entries()


def radau(endpoint, moments, denominator, polynomial):
    # 端点与三个单根给出六次精确求积。固定端点后的三次多项式由正交性确定。
    shifted = [moments[i+1]-endpoint*moments[i] for i in range(6)]
    matrix = arb_mat([[shifted[i+j] for j in range(3)] for i in range(3)])
    cubic = matrix.solve(arb_mat([[-v] for v in shifted[3:]])).entries()
    roots = arb_poly([*cubic, 1]).complex_roots(tol=1e-50)
    assert len(roots) == 3 and all(z.imag.contains(0) and 0 < z.real < 1 for z in roots)
    inside = sorted((z.real for z in roots), key=float)
    assert all(a < b for a, b in zip(inside, inside[1:]))
    nodes = [arb(endpoint), *inside]
    matrix = arb_mat([[v**k if k else 1 for v in nodes] for k in range(4)])
    weights = matrix.solve(arb_mat([[v] for v in moments[:4]])).entries()
    for k in range(4, 7):
        assert sum(w*v**k for v, w in zip(nodes, weights)).overlaps(moments[k])
    probabilities = [w*denominator(v) for v, w in zip(nodes, weights)]
    assert all(w > 0 for w in probabilities)
    response = [arb(19) if endpoint == 1 and i == 0 else neutral["efficiency"](v/(1-v), polynomial)
                for i, v in enumerate(nodes)]
    return sum(w*h for w, h in zip(probabilities, response)), probabilities


def evaluation_kernel(moments, order):
    inverse = arb_mat([[moments[i+j] for j in range(order)] for i in range(order)]).inv()
    y = arb_poly([0, 1])
    return sum(inverse[i, j]*y**(i+j) for i in range(order) for j in range(order))


def cap_atoms(polynomials, numerator, cap):
    candidates = [cap*p-numerator for p in polynomials]
    pending = [(arb(0), arb(1))]
    while pending:
        left, right = pending.pop()
        mid = (left+right)/2
        shift = arb_poly([mid, 1])
        if any(f(shift)(arb(0, right-mid)) > 0 for f in candidates):
            continue
        assert right-left > arb("1e-12")
        pending.extend(((left, mid), (mid, right)))


def verify():
    ctx.prec = 384
    times = list(neutral["model"]["shared"]["clock"]["clock"]["times"]()[:7])
    times[6] += arb("1.360900010e-8")
    denominator, moments = compact_moments(times)
    polynomial = neutral["neutral_polynomial"]()
    assert all(v > 0 for v in polynomial(arb_poly([-2, 1])))
    assert neutral["shell"]["soft"]["point"]["study"]["offset"] > 1
    lo, _ = radau(0, moments, denominator, polynomial)
    hi, relaxed = radau(1, moments, denominator, polynomial)
    assert max(w.lower() for w in relaxed) > arb(".39156768")
    # 首个三窗零点序号为 8359。用 4δ 界定三条软化边带来的非零谱位移。
    first = arb(8462)
    assert acb.zeta_zero(8359).imag > first
    delta = arb(23)/25*arb(2)**-20/first**2
    assert 4*delta < 1
    # x/(λ+x) 对 λ 的导数绝对值至多 1/(4λ)，实际与单位辐边均有 λ>1。
    upper = hi + 18*delta
    assert arb("1.0020816770805") < lo < upper < arb("1.0020816771057")

    y = arb_poly([0, 1])
    localized = [moments[i+1]-moments[i+2] for i in range(5)]
    kernels = [evaluation_kernel(moments, 4), y*(1-y)*evaluation_kernel(localized, 3)]
    cap_atoms(kernels, denominator, arb(".391568"))
    flavor_cap = arb(".390755")
    cap_atoms(kernels, denominator*(1+18*y), flavor_cap*lo)
    rest = 1-2*flavor_cap
    entropy = -(2*flavor_cap*flavor_cap.log()+rest*rest.log())/arb(2).log()
    assert entropy > arb("1.53891")
    print("平均中性响应的下上界", lo.str(18), upper.str(18))
    print("单枚的算术概率、味概率上界和味熵下界", ".391568", flavor_cap.str(7), entropy.str(14))


if __name__ == "__main__":
    verify()
