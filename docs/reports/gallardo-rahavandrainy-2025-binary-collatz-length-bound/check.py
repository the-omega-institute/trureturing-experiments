"""Gallardo-Rahavandrainy, arXiv:2510.07530v1: binary-polynomial Collatz variant. Polynomials over F2 are int bitmasks.
Checks (1) fidelity: the source's Section 3 tables (odd degree sequences and lengths); (2) exhaustively for deg <= 14:
the longest constant-degree block of the odd sequence is floor(n/2) (the source's g(n)), and the source's r_A (Theorem 1.2,
r_A = m + 1 with m least such that A_{2k} = x(x+1), A_{2k+1} = 1 for all k >= m) satisfies r_A <= 2 + floor(n^2/4)."""
def deg(a): return a.bit_length() - 1
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return r
def divmod2(a, b):
    q = 0; db = deg(b)
    while a and deg(a) >= db:
        s = deg(a) - db; q ^= 1 << s; a ^= b << s
    return q, a
X, X1, M = 0b10, 0b11, 0b111
P = mul(X, X1)
def odd_part(a):
    while a & 1 == 0: a >>= 1                       # remove x
    while True:                                     # remove x + 1
        q, r = divmod2(a, X1)
        if r: return a
        a = q
def odd_seq(a, limit=10 ** 6):
    f = odd_part(a); out = [f]
    while f != 1 and len(out) < limit:
        f = odd_part(1 ^ mul(M, f)); out.append(f)
    return out                                     # F_0, ..., F_s with F_s = 1 first
def r_A(a):
    # A_{2k+1} = F_k; A_{2k} = 1 + M F_{k-1} for k >= 1. m = least m >= 1 with A_{2k} = P and A_{2k+1} = 1 for all k >= m.
    F = odd_seq(a); s = len(F) - 1
    # for k > s every term is the fixed pair (P, 1); check k = s and below for an earlier start
    m = s + 1
    for k in range(s, 0, -1):
        if F[k] == 1 and (1 ^ mul(M, F[k - 1])) == P: m = k
        else: break
    return m + 1
def P_of(exps): return sum(1 << e for e in exps)

print("Fidelity, source Section 3.1 (T_n = M^n + 1): odd degree sequences and lengths")
for n in range(9, 17):
    A = 1
    for _ in range(n): A = mul(A, M)
    A ^= 1
    F = odd_seq(A); print(f"  n={n}: {[deg(f) for f in F]} length {len(F)}")
print("Fidelity, source Section 3.2 (x^n + x + 1):")
for n in (31, 32, 33, 34):
    F = odd_seq(P_of([n, 1, 0])); print(f"  n={n}: {[deg(f) for f in F]} length {len(F)}")
print("Fidelity, source Section 3.3 (P1, P2, P3):")
for name, e in (("P1", [14, 10, 9, 8, 3, 2, 0]), ("P2", [14, 12, 9, 6, 4, 3, 0]), ("P3", [14, 13, 9, 8, 6, 5, 3, 2, 0])):
    F = odd_seq(P_of(e)); print(f"  {name}: {[deg(f) for f in F]} length {len(F)}")

print("Exhaustive, every nonzero A with deg A = n:")
ok = True
for n in range(0, 15):
    maxr = 0; maxblock = 0
    for a in range(1 << n, 1 << (n + 1)):
        F = odd_seq(a)
        maxr = max(maxr, r_A(a))
        run = 1
        for i in range(1, len(F)):
            if deg(F[i]) == deg(F[i - 1]) and deg(F[i]) > 0: run += 1
            else: run = 1
            if deg(F[i - 1]) == n and deg(F[i]) == n: maxblock = max(maxblock, run)
        if F and deg(F[0]) == n: maxblock = max(maxblock, 1)
    bound = 2 + n * n // 4
    ok = ok and maxr <= bound
    g = maxblock if n >= 2 else None
    print(f"  n={n:2d}: max r_A={maxr} bound 2+floor(n^2/4)={bound} longest degree-n block={g} floor(n/2)={n // 2}")
print("ALL r_A <= 2 + floor(n^2/4):", ok)
