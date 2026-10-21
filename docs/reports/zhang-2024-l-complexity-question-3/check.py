"""Checks for Question 3 of P. Zhang, "The 2-complexity of even positive integers", arXiv:2411.19364v1 (Integers 24
(2024), #A109). ||n||_l = least number of l's expressing n with + and *; A_l = {n in lZ+ : ||l^m n||_l = m + ||n||_l
for all m >= 0}. Question 3 asks whether every n in lZ+ has some l^{m0} n in A_l.
(a) Reproduce the source: ||54||_2 = 8, ||2^2 * 54||_2 = 9; Theorem 2(1): ||2^m 6^r||_2 = m + 3r (small cases);
    Proposition 1: log_l n <= ||n||_l.
(b) d(m) = ||l^m n||_l - m is non-increasing in m and becomes constant within the computed range."""
import math, sys
import numpy as np

def complexity(l, N):
    INF = 10 ** 9
    c = np.full(N + 1, INF, dtype=np.int64)
    c[l] = 1
    mult = list(range(l, N + 1, l))
    for n in mult:
        if n == l: continue
        best = INF
        a = np.arange(l, n // 2 + 1, l)
        if len(a): best = int(np.min(c[a] + c[n - a]))
        d = l
        while d * d <= n:
            if n % d == 0 and (n // d) % l == 0:
                best = min(best, int(c[d] + c[n // d]))
            d += l
        c[n] = best
    return c

def main():
    bad = 0
    N2 = 1 << 15
    c2 = complexity(2, N2)
    print(f"||54||_2 = {c2[54]} (source 8), ||216||_2 = {c2[216]} (source 9)")
    bad += (c2[54] != 8) + (c2[216] != 9)
    for m in range(0, 6):
        for r in range(1, 5):
            v = 2 ** m * 6 ** r
            if v <= N2:
                ok = c2[v] == m + 3 * r
                bad += not ok
                if not ok: print("Theorem 2 mismatch", m, r, c2[v])
    print("Theorem 2(1) small cases checked")
    for l, N in [(2, N2), (3, 3 ** 9)]:
        c = c2 if l == 2 else complexity(l, N)
        prop = all(math.log(n, l) <= c[n] + 1e-9 for n in range(l, N + 1, l))
        bad += not prop
        incr = 0; stable = 0; tested = 0
        for n in range(l, 400, l):
            seq = []
            m = 0
            while l ** m * n <= N:
                seq.append(int(c[l ** m * n]) - m); m += 1
            if len(seq) < 3: continue
            tested += 1
            if any(seq[i + 1] > seq[i] for i in range(len(seq) - 1)): incr += 1
            if seq[-1] == seq[-2]: stable += 1
        bad += incr
        print(f"l={l}: Proposition 1 lower bound holds {prop}; orbits tested {tested}; d(m) increases in {incr}; last two equal in {stable}")
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
