"""OEIS A338632 (Hanna): G.f. A(x) satisfies 1 = A - x/(A - 3x/(A - 5x/(A - ...))). Conjecture: a(n) = 2 (mod 4) for n > 1.

The program computes A to order N exactly over the rationals, as the fixed point of A -> 1 + x/T_1(A). Here T_j(A) is
the coefficientwise limit of the finite tails T_{j,0} = A, T_{j,r+1} = A - (2j+1)x/T_{j+1,r}. It checks:
- every coefficient is an integer and agrees with the OEIS data;
- the defining relation holds to order N;
- a(n) = 2 (mod 4) for 2 <= n <= N;
- modulo 4, A agrees with (1 + x^2)/(1 - x) (the proof's closed form);
- modulo 3, A satisfies A^2 = A + x (reported, for the separate mod-3 conjectures of the entry).
Standard library only.
"""
from fractions import Fraction as Fr
N = 40
OEIS = [1,1,2,14,166,2714,55866,1377942,39493518,1288115570,47086272754,1906554619166,84711219819062,4098314765667082,214489189682087594,12075596389435432230,727783484288200558110,46755528594469120151010,3189788089674119448202722]

def inv(a):
    b = [Fr(0)] * (N + 1); b[0] = 1 / a[0]
    for n in range(1, N + 1):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1)) / a[0]
    return b

def tail(A, j):
    """T_j(A), truncated deep enough to be exact modulo x^(N+1)."""
    K = j + N + 2
    t = A[:]                      # T_K approximated by A
    for jj in range(K - 1, j - 1, -1):
        q = inv(t)                # 1 / T_{jj+1}
        t = [A[i] - (2 * jj + 1) * (q[i - 1] if i else 0) for i in range(N + 1)]
    return t

A = [Fr(1)] + [Fr(0)] * N
for _ in range(N + 2):
    q = inv(tail(A, 1))
    A = [Fr(1 if i == 0 else 0) + (q[i - 1] if i else 0) for i in range(N + 1)]

ints = all(x.denominator == 1 for x in A)
a = [int(x) for x in A]
ok_data = a[:len(OEIS)] == OEIS
rel = tail(A, 0)
ok_rel = rel[0] == 1 and all(v == 0 for v in rel[1:])
ok_mod4 = all(a[n] % 4 == 2 for n in range(2, N + 1))
B = [1, 1] + [2] * (N - 1)
ok_B = all((a[n] - B[n]) % 4 == 0 for n in range(N + 1))
sq = [sum(a[i] * a[n - i] for i in range(n + 1)) for n in range(N + 1)]
ok_mod3 = all((sq[n] - a[n] - (1 if n == 1 else 0)) % 3 == 0 for n in range(N + 1))
print("integral:", ints, " OEIS data:", ok_data, " relation to order %d:" % N, ok_rel)
print("a(n) = 2 mod 4 for 2..%d:" % N, ok_mod4, " A = (1+x^2)/(1-x) mod 4:", ok_B)
print("A^2 = A + x mod 3 (reported):", ok_mod3)
ok = ints and ok_data and ok_rel and ok_mod4 and ok_B
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
