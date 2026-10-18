"""Giacomelli, arXiv:2604.14391v1, Section 7 (1): is the cone C_inf = {v0 : L^i(a_n) >= 0 for all n, i} of initial
data of an infinitely log-concave sequence convex (d >= 3)? Counterexample with d = 3:
  a_{n+1} = (n+3) a_n - (3n+2) a_{n-1} + 2n a_{n-2},  n >= 2  (p = (1,-3,2), q = (3,-2,0)).
Exact rational arithmetic. Usage: python3 check.py N  (terms a_0..a_N, iterates up to N/2)
"""
import sys
from fractions import Fraction

P, Q = (1, -3, 2), (3, -2, 0)

def seq(v0, N):
    a = [Fraction(x) for x in v0]
    for n in range(2, N):
        a.append(sum((P[k] * n + Q[k]) * a[n - k] for k in range(3)))
    return a

def L(a): return [a[k] ** 2 - a[k + 1] * a[k - 1] for k in range(1, len(a) - 1)]

def inf_lc(a, depth):
    cur = a
    for i in range(depth + 1):
        if any(x < 0 for x in cur): return False, i
        if len(cur) < 3: break
        cur = L(cur)
    return True, None

def main():
    N = int(sys.argv[1]); ok = True
    U, V = (1, 1, 1), (1, 2, 4); M = tuple(Fraction(u + v, 2) for u, v in zip(U, V))
    u, v, m = seq(U, N)[:N], seq(V, N)[:N], seq(M, N)[:N]
    ok &= u == [1] * N and v == [2 ** n for n in range(N)] and m == [Fraction(1 + 2 ** n, 2) for n in range(N)]
    print(f"u_n = 1, v_n = 2^n, m_n = (1+2^n)/2 for n < {N}: {ok}")
    lu, lv = inf_lc(u, N // 2), inf_lc(v, N // 2)
    print(f"u, v: all iterates L^i >= 0 up to i = {N // 2}: {lu[0]}, {lv[0]}; L(u) = 0 and L(v) = 0: {all(x == 0 for x in L(u))}, {all(x == 0 for x in L(v))}")
    lm = L(m)
    print(f"midpoint: L(m)_1 = {lm[0]}, L(m)_2 = {lm[1]}; first failing iterate i = {inf_lc(m, N // 2)[1]}")
    ok &= lu[0] and lv[0] and lm[0] == Fraction(-1, 4) and lm[1] == Fraction(-1, 2)
    print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)

main()
