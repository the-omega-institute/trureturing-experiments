"""De Castro's p-adic order formula: OEIS A343249-A343253, A345531, and INTEGERS 22 (2022) #A61,
Conjectures 1 and 2.

L_p(k0, n) is the lowest common denominator of {binomial(n, k)/n : 0 < k <= k0, p does not divide k}.
The formula v_p(n) = log_p(n / L) holds iff n / L = p^v_p(n) (L divides n).

Usage: python3 check.py N Q
  N: range of n for A343249-A343253 and the INTEGERS conjectures
  Q: composite q <= Q for A345531
"""
import sys
from math import gcd, comb

DATA = {
    2: [1,1,3,1,5,3,7,1,9,5,11,3,13,7,5,1,17,9,19,5,7,11,23,3,25,13,27,7,29,5,31,1,11,17,7,9,37,19,13,5,41,7,43,11,9,23,47,3,49,25,17,13,53,27,11,7,19,29,59,5,61,31,9,1,13,11,67,17,23,7,71,9,73,37,25,19,11,13,79,5],
    3: [1,2,1,4,5,2,7,8,1,5,11,4,13,7,5,16,17,2,19,5,7,11,23,8,25,13,1,7,29,5,31,32,11,17,7,4,37,19,13,8,41,7,43,11,5,23,47,16,49,25,17,13,53,2,11,8,19,29,59,5,61,31,7,64,13,11,67,17,23,7,71,8,73,37,25,19,11,13,79,16],
    5: [1,2,3,4,1,3,7,8,9,2,11,4,13,7,3,16,17,9,19,4,7,11,23,8,1,13,27,7,29,3,31,32,11,17,7,9,37,19,13,8,41,7,43,11,9,23,47,16,49,2,17,13,53,27,11,8,19,29,59,4,61,31,9,64,13,11,67,17,23,7,71,9,73,37,3,19,11,13,79,16],
    7: [1,2,3,4,5,3,1,8,9,5,11,4,13,2,5,16,17,9,19,5,3,11,23,8,25,13,27,4,29,5,31,32,11,17,5,9,37,19,13,8,41,3,43,11,9,23,47,16,1,25,17,13,53,27,11,8,19,29,59,5,61,31,9,64,13,11,67,17,23,5,71,9,73,37,25,19,11,13,79,16],
    11: [1,2,3,4,5,3,7,8,9,5,1,4,13,7,5,16,17,9,19,5,7,2,23,8,25,13,27,7,29,5,31,32,3,17,7,9,37,19,13,8,41,7,43,4,9,23,47,16,49,25,17,13,53,27,5,8,19,29,59,5,61,31,9,64,13,3,67,17,23,7,71,9,73,37,25,19,7,13,79,16],
}
A345531 = [3,4,7,8,13,16,19,23,25,31,32,41,43,47,49,59,61,64,71,73,79,81,89,97,101,103,107,109,113,121,128,137,139,149,151,157,163,167,169,179,181,191,193,197,199,211,223,227,229,233,239,241,243,256,263,269,271,277]

def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def is_prime(n): return n > 1 and factor(n) == {n: 1}
def is_prime_power(n): return n > 1 and len(factor(n)) == 1

def v(q, n):
    e = 0
    while n % q == 0: n //= q; e += 1
    return e

def den(n, k):  # denominator of binomial(n, k)/n in lowest terms
    return n // gcd(n, comb(n, k))

def lcm(a, b): return a * b // gcd(a, b)

def least_k0(p, n, excluded):
    """least k0 in [1, n] with n / L(k0) = p^v_p(n), where L runs over k <= k0 with not excluded(k); 0 if none.
    L(k0) divides L(k0 + 1), so once L does not divide the target no larger k0 works."""
    target = n // p ** v(p, n)
    L = 1; C = 1
    for k in range(1, n + 1):
        C = C * (n - k + 1) // k
        if not excluded(k):
            d = n // gcd(n, C)
            L = L * d // gcd(L, d)
        if L == target: return k
        if target % L: return 0
    return 0

def gpp(p, n):  # greatest power of a prime different from p dividing n
    return max([1] + [r ** e for r, e in factor(n).items() if r != p])

def spp(m):  # smallest prime power > m
    k = m + 1
    while not is_prime_power(k): k += 1
    return k

def main():
    N, Q = int(sys.argv[1]), int(sys.argv[2])
    ok = True
    # A343249-A343253: data, and a(n) = gpp(p, n)
    for p, data in DATA.items():
        d = [least_k0(p, n, lambda k: k % p == 0) for n in range(1, len(data) + 1)] == data
        g = all(least_k0(p, n, lambda k: k % p == 0) == gpp(p, n) for n in range(1, N + 1))
        print(f"A3432xx p={p}: data {d}; a(n) = greatest prime power (prime != p) dividing n for n <= {N}: {g}")
        ok &= d and g
    # the valuation lemma: for every prime r and k, the r-part of den(n,k) is at most r^min(v_r(n), v_r(k)),
    # and at k = r^e (e = v_r(n)) it is exactly r^e
    lem = True
    for n in range(1, min(N, 600) + 1):
        f = factor(n)
        for k in range(1, n + 1):
            dk = den(n, k)
            for r in factor(dk):
                lem &= v(r, dk) <= min(v(r, n), v(r, k))
        for r, e in f.items():
            lem &= v(r, den(n, r ** e)) == e
    print(f"lemma (r-part of den(n,k) <= r^min(v_r(n), v_r(k)); exactly r^v_r(n) at k = r^v_r(n)) for n <= {min(N, 600)}: {lem}")
    ok &= lem
    # INTEGERS 22 #A61, Conjecture 1: for prime p, L(S_min) = L(S_max) = n / p^v_p(n)
    c1 = True
    for p in [2, 3, 5, 7, 11, 13]:
        for n in range(1, N + 1):
            target = n // p ** v(p, n)
            Lmax = 1
            for k in range(1, n + 1):
                if gcd(k, p) == 1: Lmax = lcm(Lmax, den(n, k))
            Lmin = 1
            for m in range(1, n + 1):
                if is_prime_power(m) or m == 1:
                    if gcd(m, p) == 1: Lmin = lcm(Lmin, den(n, m))
            c1 &= (Lmin == target == Lmax)
    print(f"INTEGERS Conjecture 1: L(S_min) = L(S_max) = n/p^v_p(n) for p <= 13, n <= {N}: {c1}")
    ok &= c1
    # INTEGERS 22 #A61, Conjecture 2: composite q, prime pj | q, n = q^c pj: L(S) = 1 and n/L != q^v_q(n)
    c2 = True; cnt = 0
    for q in range(4, 31):
        if is_prime(q): continue
        for pj in factor(q):
            for c in range(0, 3):
                n = q ** c * pj
                Lmax = 1; C = 1
                for k in range(1, n + 1):
                    C = C * (n - k + 1) // k
                    if gcd(k, q) == 1: Lmax = lcm(Lmax, n // gcd(n, C))
                c2 &= (Lmax == 1) and (n != q ** v(q, n)) and v(q, n) == c
                cnt += 1
    print(f"INTEGERS Conjecture 2: L(S_max) = 1 and the formula fails at n = q^c pj ({cnt} cases, composite q <= 30, c <= 2): {c2}")
    ok &= c2
    # A345531: data, and for composite q the first n with no k0 in [1, n] is q * A345531(i), i = index of minFac(q)
    primes = [m for m in range(2, 400) if is_prime(m)]
    d345 = [spp(p) for p in primes[:len(A345531)]] == A345531
    c3 = True; mism = []
    for q in range(4, Q + 1):
        if is_prime(q): continue
        p = min(factor(q)); i = primes.index(p) + 1; pred = q * spp(primes[i - 1])
        first = next(n for n in range(1, pred + 1) if least_k0(q, n, lambda k: k % q == 0) == 0 or n == pred)
        if not (first == pred and least_k0(q, pred, lambda k: k % q == 0) == 0):
            c3 = False; mism.append(q)
    print(f"A345531: data {d345}; first failure at q * a(i) for every composite q <= {Q}: {c3} {mism[:5]}")
    # the other reading (k0 = n only) fails already at n = q
    def holds_at_n(q, n):
        target = n // q ** v(q, n); L = 1; C = 1
        for k in range(1, n + 1):
            C = C * (n - k + 1) // k
            if k % q: L = lcm(L, n // gcd(n, C))
        return L == target
    other = all(next(n for n in range(1, q + 1) if not holds_at_n(q, n)) == q for q in range(4, Q + 1) if not is_prime(q))
    print(f"A345531, reading k0 = n: the first failure is n = q for every composite q <= {Q}: {other}")
    ok &= other
    ok &= d345 and c3
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
