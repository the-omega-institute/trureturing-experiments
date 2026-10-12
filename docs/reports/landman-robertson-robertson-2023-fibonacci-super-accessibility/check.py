#!/usr/bin/env python3
"""Exact check of the rotation 2-coloring refuting 2-super-accessibility of the Fibonacci set.

Coloring: n >= 1 has color 2 iff frac(alpha*n) < 1/16, alpha = (47 + sqrt 5)/58; otherwise color 1.
Decided with integer arithmetic only (math.isqrt). Reports, for 1 <= n <= N:
  - the longest color-2 F-diffsequence (consecutive gaps in F),
  - the longest color-1 F-arithmetic progression, per gap,
  - the longest run of consecutive color-1 integers (gap of color 2),
  - the Lemma 1 residue cycle of 23 f_n + f_{n+1} mod 29.
Usage: check.py [N]   (default 200000). Exit 0 iff the observed maxima are <= 2 and <= 41.
"""
import math, sys

def color2(n: int) -> bool:
    # alpha*n = (47n + n*sqrt5)/58 ; floor uses floor(n*sqrt5) = isqrt(5 n^2)
    s = math.isqrt(5 * n * n)
    k = (47 * n + s) // 58
    # frac < 1/16  <=>  8 n sqrt5 < 464 k - 376 n + 29 =: R
    R = 464 * k - 376 * n + 29
    return R > 0 and 320 * n * n < R * R

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    fib = [1, 2]
    while fib[-1] + fib[-2] <= N:
        fib.append(fib[-1] + fib[-2])
    c2 = [False] + [color2(n) for n in range(1, N + 1)]
    # longest color-2 F-diffsequence ending at n
    L = [0] * (N + 1)
    best2 = 0
    for n in range(1, N + 1):
        if c2[n]:
            m = 1
            for f in fib:
                if f >= n:
                    break
                if c2[n - f] and L[n - f] + 1 > m:
                    m = L[n - f] + 1
            L[n] = m
            best2 = max(best2, m)
    # longest color-1 F-AP per gap
    per = {}
    for d in fib:
        run = [0] * (N + 1)
        b = 0
        for n in range(1, N + 1):
            if not c2[n]:
                run[n] = (run[n - d] if n - d >= 1 else 0) + 1
                if run[n] > b:
                    b = run[n]
        per[d] = b
    best1 = max(per.values())
    gap, cur = 0, 0
    for n in range(1, N + 1):
        cur = 0 if c2[n] else cur + 1
        gap = max(gap, cur)
    a, b = 1, 24  # 23 f_0 + f_1, 23 f_1 + f_2
    cyc = []
    for _ in range(9):
        cyc.append(a % 29)
        a, b = b, (a + b)
    print(f"N={N} fib_gaps={len(fib)} color2_count={sum(c2)}")
    print(f"max_color2_F_diffsequence={best2}")
    print(f"max_color1_F_AP={best1} per_gap={ {d: per[d] for d in fib[:8]} }")
    print(f"max_consecutive_color1={gap}")
    print(f"residues_23f_n_plus_f_n1_mod29_n0_to_8={cyc}")
    ok = best2 <= 2 and best1 <= 41
    print("bad=" + ("0" if ok else "1"))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
