"""严格三循环五壳共同解的候选低线与二维零空间中的纯态方向。"""

from bisect import bisect_left, bisect_right
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

S = 1 << 28


def floor(a):
    return int((a * S).lower().floor().fmpz())


def ceil(a):
    return int((a * S).upper().ceil().fmpz())


def divup(a, b):
    return -((-a) // b)


def adj(m):
    return m.conjugate().transpose()


def sq(m):
    return sum(abs(z) ** 2 for z in m.entries())


def row(m, i):
    return acb_mat([[m[i, k] for k in range(m.ncols())]])


def kernel_frame(rows, floor_energy, cutoff_gev):
    I = acb_mat([[int(i == j) for j in range(7)] for i in range(7)])
    B = acb_mat([v.tolist()[0] for v in rows])
    J = adj(B) * (B * adj(B)).inv()
    F = I - J * B
    basis = []
    for _ in range(2):
        idx = max(range(7), key=lambda i: F[i, i].real.lower())
        n = acb_mat([[F[i, idx]] for i in range(7)])
        n = n / sq(n).sqrt()
        basis.append(n)
        F = F - n * adj(n)
    eps = (arb(cutoff_gev) / arb("1.22089e19") - floor_energy).sqrt()
    cap = eps * sum(sum(abs(J[k, i]) ** 2 for k in range(7)).sqrt() for i in range(5))
    return basis, J, eps, cap


def shellprep(ids, model, t0, lam, floor_energy, T, cutoff_gev=3000):
    basis, J, eps, cap = kernel_frame(
        [row(s["zero"], i) for s, i in zip(model, ids, strict=True)],
        floor_energy,
        cutoff_gev,
    )
    out = []
    for s in model:
        terms = []
        for i in s["triv"]:
            f = row(s["zero"], i)
            a, b = [(f * n)[0, 0] for n in basis]
            cross = a.conjugate() * b
            # 圆球上的任意方向给出正交补的主特征投影。同一方向用于所有壳。
            arr = [
                (abs(a) ** 2 + abs(b) ** 2) / 2,
                cross.real,
                cross.imag,
                (abs(a) ** 2 - abs(b) ** 2) / 2,
            ]
            coeff = [floor(t0 * x) for x in arr]
            quant = sum(max(0, ceil(t0 * x) - c) for x, c in zip(arr, coeff))
            err = eps * sum(abs(v) for v in (f * J).entries())
            # 当次特征值至多 lam 时，用最高源迹 t0 代替投影主特征值所引入的算子范数误差至多 cap²+lam。
            delta = (
                ceil(cap**2 + lam + 2 * t0.sqrt() * err + err**2)
                + quant
                + ceil(floor_energy)
            )
            w = s["w"][i]
            terms.append((coeff, delta, floor(w), ceil(w)))
        wbig = [-s["w"][i] for i in s["big"]]
        out.append(
            (
                terms,
                ceil(sum(wbig)),
                ceil(T**2 * max(x.upper() for x in wbig)),
                ceil(T**4 * max(x.upper() for x in wbig)),
            )
        )
    return out


# 在有理格点包围各场的质量修正，按符号分开累加。
def make_lookup(shell):
    mesh = shell["defects"]
    points = 1024
    low = []
    hi = []
    for i in range(points + 1):
        x = arb(i) / points
        ka = max(0, bisect_left(mesh.x, x) - 1)
        kb = min(len(mesh.x) - 1, bisect_right(mesh.x, x))
        assert mesh.x[ka] <= x <= mesh.x[kb]
        low.append([floor(v) for v in mesh.parts[ka]])
        hi.append([ceil(v) for v in mesh.parts[kb]])
    return low, hi


def valid(bounds, coeffs, lookup, scales):
    terms, sw, h2, h4 = coeffs
    ll = [0, 0]
    hh = [0, 0]
    rad = S // (2 ** bounds[3])
    x = bounds[:3]
    upl = divup
    for c, delta, wl, wh in terms:
        mid = S * c[0] + sum(d * t for d, t in zip(c[1:], x))
        wid = rad * sum(abs(d) for d in c[1:])
        lo = max(0, (mid - wid) // S - delta)
        hi = min(S, upl(mid + wid, S) + delta)
        low, mh = lookup
        L = low[(lo * 1024) // S]
        H = mh[upl(hi * 1024, S)]
        for h in (0, 1):
            power = 2 + 2 * h
            sig = -1 if h == 0 else 1
            ind = 2 * h
            ek = L[ind] - H[ind + 1]
            ex = H[ind] - L[ind + 1]
            fl, fh = scales[h]
            a, b = sorted((sig * fl * lo**power, sig * fh * hi**power))
            a = a // (S**power) + ek
            b = upl(b, S**power) + ex
            vals = (wl * a, wl * b, wh * a, wh * b)
            ll[h] += min(vals)
            hh[h] += max(vals)
    A = S**2 - hh[0] - E2
    B = S**2 - ll[0] + E2
    C = ll[1] - E4 - 1
    D = hh[1] + E4
    # A、B 界住主字符的二次矩乘牛顿系数，C、D 界住四次矩乘真空系数，均以 S² 缩放。
    if min(B, D) < 0:
        return False
    k, ww = scales
    if A > k[1] * h2 or C > ww[1] * h4:
        return False
    return max(0, A) ** 2 * ww[0] <= sw * D * k[1] ** 2


E2 = divup(12001 * S * S, 1000000)
E4 = divup(2001 * S * S, 1000000)


def cover(data, table, monom, depthmax=8):
    todo = [(0, 0, 0, 0)]
    visited = 0
    while todo:
        state = todo.pop()
        a, b, c, depth = state
        rad = S >> depth
        visited += 1
        if (
            sum(max(0, abs(v) - rad) ** 2 for v in (a, b, c)) > S * S
            or sum((abs(v) + rad) ** 2 for v in (a, b, c)) < S * S
        ):
            continue
        if any(
            not valid(state, s, l, m)
            for s, l, m in zip(data, table, monom, strict=True)
        ):
            continue
        if depth >= depthmax:
            return False, visited
        todo.extend(
            (a + h * (rad // 2), b + k * (rad // 2), c + l * (rad // 2), depth + 1)
            for h, k, l in product((1, -1), repeat=3)
        )
    return True, visited


def experiment(root, fifth, columns, levels):
    choices = run_path(str(root / "ch39-seven-cyclic-weak-horizon.py"))["CANDIDATES"][
        :5
    ]
    sieve = run_path(str(root / "ch39-six-cyclic-line-sieve.py"))
    blocks, (t0, T, model) = sieve["prepare_shells"](
        root, fifth, columns, levels, choices, arb(0)
    )
    m = blocks[0]
    det = sum(
        abs(m[i, 0] * m[j, 1] - m[i, 1] * m[j, 0]) ** 2
        for i in range(7)
        for j in range(i)
    )
    lam = 2 * det / t0
    assert 0 < lam < arb("3e-11")
    table = [make_lookup(s) for s in model]
    monom = [[(floor(s[k]), ceil(s[k])) for k in ("ck", "cw")] for s in model]
    out = []
    for ids in product(*choices):
        data = shellprep(ids, model, t0, lam, levels[0], T)
        excluded, _ = cover(data, table, monom)
        if not excluded:
            out.append(ids)
    for line, cutoff in ((0, 10000), (15, 100000000)):
        data = shellprep((line,) * 5, model, t0, lam, levels[0], T, cutoff)
        excluded, _ = cover(data, table, monom, depthmax=9)
        assert excluded
        print(
            "五个引力双矩全部匹配时，至少一壳此谱线的截止超过给定 GeV",
            line + 1,
            cutoff,
            flush=True,
        )
    assert len(out) == 20 and all(set(ids) == {0, 15} for ids in out)
    height = lambda n: acb.zeta_zero(n).imag
    span = height(2_306_888_100_096) - height(2_306_888_100_092)
    assert 0 < span < arb("1.037551")
    assert (arb(100000) / 3).log() > arb("10.0374") * span
    print(
        "五壳候选、排除、未排除及一基选线",
        108,
        108 - len(out),
        len(out),
        [tuple(i + 1 for i in ids) for ids in out],
        flush=True,
    )


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
