"""Exact check of Hemmer-Straub-Westrem, arXiv:2509.02796v1 (EJC 33(1) P1.57), Conjecture 2.3.

For lambda = (a_1, ..., a_r) |- n, Ev(lambda) is the multiset of 2^r cycle types
obtained by replacing each a_i by (2a_i) or (a_i, a_i).  The conjecture:
  sum_{t in Ev} sum_{mu in R3(2n)} (-1)^{len t} chi^mu_t = sum_{t in Ev} chi^{(n,n)}_t,
R3(2n) = partitions of 2n with at most three parts, all even.

Characters are computed twice:
  (1) Murnaghan-Nakayama rule on beta-sets (rim-hook removal), any mu;
  (2) the three-variable Frobenius formula chi^mu_t = [x^{mu+delta}] (a_delta * p_t)
      for mu with at most three parts (the convention proposed for Lean).
Prediction registered before the run: (1) = (2) on every pair used, and the
identity holds for every lambda |- n, 1 <= n <= N (default N = 11).
Also reports the reduction of the proposed proof: both sides equal
2^r * CT((1 - t) prod_i (t^{a_i} + 1 + t^{-a_i})).

Usage: python3 check.py [N]   (standard library only)
"""
import sys
from functools import lru_cache
from collections import defaultdict


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


@lru_cache(maxsize=None)
def mn_char(mu, rho):
    """Murnaghan-Nakayama: chi^mu at cycle type rho (tuples of positive parts)."""
    if sum(mu) != sum(rho):
        raise ValueError
    if not rho:
        return 1
    k, rest = rho[0], rho[1:]
    L = len(mu)
    beta = [mu[i] + (L - 1 - i) for i in range(L)]
    bset = set(beta)
    total = 0
    for b in beta:
        nb = b - k
        if nb < 0 or nb in bset:
            continue
        sign = (-1) ** sum(1 for c in beta if nb < c < b)
        newbeta = sorted([c for c in beta if c != b] + [nb], reverse=True)
        newmu = tuple(x for x in (newbeta[i] - (L - 1 - i) for i in range(L)) if x > 0)
        total += sign * mn_char(newmu, rest)
    return total


def polymul(p, q):
    r = defaultdict(int)
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            r[tuple(a + b for a, b in zip(e1, e2))] += c1 * c2
    return {e: c for e, c in r.items() if c}


def power_sum(k):
    return {(k, 0, 0): 1, (0, k, 0): 1, (0, 0, k): 1}


VANDERMONDE = polymul(polymul({(1, 0, 0): 1, (0, 1, 0): -1}, {(1, 0, 0): 1, (0, 0, 1): -1}),
                      {(0, 1, 0): 1, (0, 0, 1): -1})


@lru_cache(maxsize=None)
def frob_poly(rho):
    p = VANDERMONDE
    for k in rho:
        p = polymul(p, power_sum(k))
    return p


def frob_char(mu3, rho):
    """[x^{mu+delta}] (a_delta * p_rho) in three variables, delta = (2, 1, 0)."""
    e = (mu3[0] + 2, mu3[1] + 1, mu3[2])
    return frob_poly(rho).get(e, 0)


def ev(lam):
    out = [()]
    for a in lam:
        out = [t + (2 * a,) for t in out] + [t + (a, a) for t in out]
    return [tuple(sorted(t, reverse=True)) for t in out]


def r3(m):
    return [mu for mu in partitions(m) if len(mu) <= 3 and all(x % 2 == 0 for x in mu)]


def pad3(mu):
    return tuple(mu) + (0,) * (3 - len(mu))


def ct_reduction(lam):
    q = {0: 1}
    for a in lam:
        f = {a: 1, 0: 1, -a: 1}
        nq = defaultdict(int)
        for e1, c1 in q.items():
            for e2, c2 in f.items():
                nq[e1 + e2] += c1 * c2
        q = nq
    return (2 ** len(lam)) * (q.get(0, 0) - q.get(1, 0))


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    ok = True
    checked = cross = 0
    for n in range(1, N + 1):
        R = r3(2 * n)
        for lam in partitions(n):
            E = ev(lam)
            lhs = rhs = 0
            for t in E:
                s = (-1) ** len(t)
                for mu in R:
                    c = mn_char(mu, t)
                    f = frob_char(pad3(mu), t)
                    cross += 1
                    if c != f:
                        ok = False
                        print("CROSS-CHECK FAIL", mu, t, c, f)
                    lhs += s * c
                c2 = mn_char((n, n), t)
                f2 = frob_char((n, n, 0), t)
                cross += 1
                if c2 != f2:
                    ok = False
                    print("CROSS-CHECK FAIL", (n, n), t, c2, f2)
                rhs += c2
            ct = ct_reduction(lam)
            checked += 1
            if not (lhs == rhs == ct):
                ok = False
                print("IDENTITY FAIL", n, lam, lhs, rhs, ct)
            if n <= 4 or lam in [(2, 1), (4, 3, 2, 1)]:
                print(f"n={n} lambda={lam}: LHS = RHS = {lhs} (2^r CT = {ct})")
        print(f"n={n}: all {sum(1 for _ in partitions(n))} partitions checked")
    print(f"partitions checked: {checked}; character cross-checks (MN vs Frobenius): {cross}")
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
