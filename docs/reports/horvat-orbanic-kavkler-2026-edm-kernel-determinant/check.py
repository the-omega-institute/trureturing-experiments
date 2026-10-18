"""Horvat–Orbanić–Kavkler arXiv:2607.09671v1, Section 4 open problem: closed form for det M,
where q(x - y)^n = phi(x)^T M phi(y) with the source's basis of V (Definition 1):
  A: x^alpha, |alpha| <= n-1;  B: x^alpha, |alpha| = n;  C: |x|^{2r} x^alpha, r = 1..n, |alpha| = n - r.
M is computed exactly (sympy rationals) from the monomial expansion, then det M is compared with the
proposed formula
  det M = (-1)^{C(n+d-1,d) + n C(n+d-1,d-1)} P_n prod_{k<n} P_k^2,
  P_k = D_k^{-1} prod_{l=0}^{floor(k/2)} [2^k n! (n-l+d/2)_l / (n-k+l)!]^{h_{d,k-2l}},
  r_k = C(k+d-1,d-1), h_{d,s} = r_s - r_{s-2}, D_k = prod_{|alpha|=k} alpha!.
Also checked: d = 1 reproduces prod_a C(2n,a) (source Remark 5); n = 1 gives (-1)^{d+1} 2^d.
"""
import itertools, sys
import sympy as sp


def multi(d, k):
    return [a for a in itertools.product(range(k + 1), repeat=d) if sum(a) == k]


def det_M(d, n):
    xs = sp.symbols(f"x0:{d}"); ys = sp.symbols(f"y0:{d}")
    q = lambda v: sum(t ** 2 for t in v)
    basis = []
    for k in range(n):
        for a in multi(d, k): basis.append((0, a))
    for a in multi(d, n): basis.append((0, a))
    for r in range(1, n + 1):
        for a in multi(d, n - r): basis.append((r, a))
    mons = [a for k in range(2 * n + 1) for a in multi(d, k)]
    midx = {a: i for i, a in enumerate(mons)}
    N, Nm = len(basis), len(mons)
    T = sp.zeros(N, Nm)
    for i, (r, a) in enumerate(basis):
        poly = sp.Poly(sp.expand(q(xs) ** r * sp.prod([xs[j] ** a[j] for j in range(d)])), *xs)
        for mon, c in poly.terms(): T[i, midx[tuple(mon)]] += c
    K = sp.zeros(Nm, Nm)
    poly = sp.Poly(sp.expand(q([xs[i] - ys[i] for i in range(d)]) ** n), *xs, *ys)
    for mon, c in poly.terms(): K[midx[tuple(mon[:d])], midx[tuple(mon[d:])]] += c
    G = (T * T.T).inv()
    M = G * T * K * T.T * G
    assert T.T * M * T == K, "kernel not in V (x) V"
    return M.det(), N


def formula(d, n):
    C = sp.binomial
    r = lambda k: C(k + d - 1, d - 1) if k >= 0 else 0
    h = lambda s: r(s) - r(s - 2)
    def P(k):
        Dk = sp.prod([sp.prod([sp.factorial(ai) for ai in a]) for a in multi(d, k)])
        val = sp.Integer(1)
        for l in range(k // 2 + 1):
            base = 2 ** k * sp.factorial(n) * sp.rf(sp.Rational(n - l) + sp.Rational(d, 2), l) / sp.factorial(n - k + l)
            val *= base ** h(k - 2 * l)
        return val / Dk
    sign = (-1) ** (C(n + d - 1, d) + n * C(n + d - 1, d - 1))
    return sign * P(n) * sp.prod([P(k) ** 2 for k in range(n)])


def main():
    fails = []
    cases = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (5, 2)]
    for d, n in cases:
        dm, N = det_M(d, n)
        f = formula(d, n)
        ok = sp.simplify(dm - f) == 0
        print(f"d={d} n={n} N={N}: det M = {dm} formula = {f} {'OK' if ok else 'MISMATCH'}")
        if not ok: fails.append((d, n))
        if d == 1 and dm != sp.prod([sp.binomial(2 * n, a) for a in range(2 * n + 1)]): fails.append(("d1", n))
        if n == 1 and dm != (-1) ** (d + 1) * 2 ** d: fails.append(("n1", d))
    # control: a wrong variant (drop the factorial normalization D_k) must fail somewhere
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
