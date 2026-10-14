#!/usr/bin/env python3
"""Exact check of the proposed answer to Gudder, "Dual Quantum Instruments and Sub-Observables",
arXiv:2302.12243v1, Sec. 3: "we conjecture that U need not be a Sob effect algebra" and Example 3
("If L_x(rho) = a_x o rho is an arbitrary Lueders instrument, we do not know whether U_L is a Sob
effect algebra").

Witness: H = C^2, outcomes {0,1}, K_0 = diag(4/5,3/5), K_1 = diag(3/5,4/5) (positive square roots of
the Lueders effects a_x = K_x^2, a_0 + a_1 = I); dual F(C)_x = K_x C K_x.
E = (1/400)[[151,50],[50,151]]; A = F(E). Then A + A is a sub-observable (positive entries, total
<= I), but its only possible preimage is 2E, which is not an effect. So (S3) fails.
Also the family u^2 + v^2 = 1, u != v, E_r = ((7-r)/16) I + (1/8) X, r = 2uv.
"""
import sys
from fractions import Fraction as F

def mul(A, B): return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def add(A, B): return [[A[i][j]+B[i][j] for j in range(2)] for i in range(2)]
def sc(c, A): return [[c*A[i][j] for j in range(2)] for i in range(2)]
def psd(A):  # real symmetric 2x2: PSD iff diag >= 0 and det >= 0
    return A[0][1] == A[1][0] and A[0][0] >= 0 and A[1][1] >= 0 and A[0][0]*A[1][1] - A[0][1]*A[1][0] >= 0
def pd(A):
    return A[0][1] == A[1][0] and A[0][0] > 0 and A[0][0]*A[1][1] - A[0][1]*A[1][0] > 0
I = [[F(1), F(0)], [F(0), F(1)]]
def sub(A, B): return add(A, sc(F(-1), B))

ok = True
def run(u, v, E, label):
    global ok
    K0 = [[u, F(0)], [F(0), v]]; K1 = [[v, F(0)], [F(0), u]]
    inst = add(mul(K0, K0), mul(K1, K1)) == I            # sum K_x^2 = I (Lueders, K_x PSD)
    effE = pd(E) and pd(sub(I, E))                       # 0 < E < I
    A0 = mul(mul(K0, E), K0); A1 = mul(mul(K1, E), K1)
    sumAA = add(sc(F(2), A0), sc(F(2), A1))
    sob = psd(sc(F(2), A0)) and psd(sc(F(2), A1)) and psd(sub(I, sumAA))
    twoE = sc(F(2), E)
    not_effect = not psd(sub(I, twoE))                   # 2E is not <= I
    # uniqueness of preimage: K_0 C K_0 = 2 K_0 E K_0 with K_0 invertible forces C = 2E
    det0 = u*v
    line = inst and effE and sob and not_effect and det0 != 0
    ok &= line
    print(f"{label}: Lueders instrument={inst}, E effect={effE}, A+A sub-observable={sob}, "
          f"unique preimage 2E not an effect={not_effect} {'OK' if line else 'FAIL'}")

run(F(4,5), F(3,5), [[F(151,400), F(50,400)], [F(50,400), F(151,400)]], "witness u=4/5")
# family members with rational u, v (Pythagorean triples), u != v
for (p, q, h) in [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29)]:
    u, v = F(p, h), F(q, h); r = 2*u*v
    E = add(sc((7 - r)/16, I), sc(F(1, 8), [[F(0), F(1)], [F(1), F(0)]]))
    run(u, v, E, f"family u={u} v={v}")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
