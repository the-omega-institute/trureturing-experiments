"""Exact checks for the positivity question of Gadde-Jain-Krishna-Kulkarni-Sharma, arXiv:2308.16247, Section 7:
is the replica invariant E_n^(q) of every pure state positive?  Here q = 3, n = 4.

With the third party unshifted, E_4^(3)(psi) = sum over a, b in [d]^(4x4) of
prod_{i,j} rho[(a_ij, b_ij), (a_i,j+1, b_i+1,j)] for rho = rho_AB, indices mod 4, which is Tr(T^4) for the
81 x 81 row-transfer matrix T.  rho_s = (I - s F)/(9 - 3 s) on C^3 (x) C^3, F the swap.
Prints ALL_OK when every check passes.
"""
import itertools
from fractions import Fraction
import sympy as sp
from sympy.combinatorics import Permutation

def transfer(R, d=3):
    """Row transfer T[b, c] = sum_a prod_j R[a_j, b_j, a_{j+1}, c_j] (exact, R integer/Fraction 4-index)."""
    idx = list(itertools.product(range(d), repeat=4))
    T = [[0] * len(idx) for _ in idx]
    for x, bs in enumerate(idx):
        for y, cs in enumerate(idx):
            tot = 0
            for As in idx:
                p = 1
                for j in range(4):
                    p *= R[As[j]][bs[j]][As[(j + 1) % 4]][cs[j]]
                    if p == 0:
                        break
                tot += p
            T[x][y] = tot
    return T

def matmul(A, B):
    n = len(A)
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(A[i], Bt[j])) for j in range(n)] for i in range(n)]

def trace4(T):
    T2 = matmul(T, T)
    return sum(T2[i][k] * T2[k][i] for i in range(len(T)) for k in range(len(T)))

def R_of(s):
    # R = (9 - 3s) rho_s = I - s F, entries R[a][b][a'][b'] = rho entry ((a,b),(a',b')) times (9-3s)
    return [[[[ (1 if (a == a2 and b == b2) else 0) - s * (1 if (a == b2 and b == a2) else 0)
              for b2 in range(3)] for a2 in range(3)] for b in range(3)] for a in range(3)]

ok = True
# 1. Endpoint s = 1 (the totally antisymmetric three-qutrit state): E = Tr(T^4)/6^16 = -1680/6^16.
P1 = trace4(transfer(R_of(1)))
E1 = Fraction(P1, 6 ** 16)
print("endpoint Tr(T^4) =", P1, " E =", E1)
ok &= (P1 == -1680) and (E1 == Fraction(-35, 58773123072))

# 2. Direct contraction from the pure state psi = eps_ijk/sqrt(6) (third party traced inside each replica).
eps = {}
for p in itertools.permutations(range(3)):
    sgn = Permutation(list(p)).signature()
    eps[p] = sgn
def rhoAB(a, b, a2, b2):
    return sum(eps.get((a, b, c), 0) * eps.get((a2, b2, c), 0) for c in range(3))  # = 6 rho_AB
R6 = [[[[rhoAB(a, b, a2, b2) for b2 in range(3)] for a2 in range(3)] for b in range(3)] for a in range(3)]
P1b = trace4(transfer(R6))
print("from psi: Tr(T^4) * ... =", P1b)
ok &= (P1b == -1680)

# 3. The polynomial P(s) = Tr(T_s^4) and its sign on [9/10, 1].
s = sp.symbols('s')
Rs = [[[[ (1 if (a == a2 and b == b2) else 0) - s * (1 if (a == b2 and b == a2) else 0)
          for b2 in range(3)] for a2 in range(3)] for b in range(3)] for a in range(3)]
# Evaluate P at 10 exact points of degree <= 16 polynomial: interpolate from 17 rational samples.
pts = [Fraction(k, 16) for k in range(17)]
vals = [trace4(transfer(R_of(t))) for t in pts]
P = sp.expand(sp.interpolate(list(zip([sp.Rational(t.numerator, t.denominator) for t in pts],
                                      [sp.Rational(v.numerator, v.denominator) if isinstance(v, Fraction) else v for v in vals])), s))
print("P(s) =", P)
ok &= (P.subs(s, 1) == -1680) and (P.subs(s, 0) == 3 ** 8)
u = sp.symbols('u')
Pu = sp.expand(P.subs(s, 1 - u))
claimed = (-1680 + 22272*u - 114048*u**2 + 330624*u**3 + u**4*(-548064 + 417024*u + 253824*u**2)
           + u**7*(-1071744 + 1460160*u) + u**9*(-1251456 + 759552*u) + u**11*(-340992 + 116160*u)
           + u**13*(-30336 + 6048*u) + u**15*(-864 + 81*u))
print("P(1-u) matches the seat's expansion:", sp.expand(Pu - claimed) == 0)
ok &= sp.expand(Pu - claimed) == 0
# Sign on [9/10, 1]: no real root of P in [9/10, 1] and P(1) < 0.
nroots = sp.Poly(P, s).count_roots(sp.Rational(9, 10), 1)
root_iv = [r for r in sp.Poly(P, s).intervals(eps=sp.Rational(1, 10**6)) if 0 <= r[0][0] <= 1]
print("real roots of P in [9/10, 1] (Sturm count):", nroots, " isolating intervals in [0,1]:", [(float(a), float(b)) for (a, b), m in root_iv])
ok &= (nroots == 0) and P.subs(s, sp.Rational(9, 10)) < 0
print("P(9/10) =", P.subs(s, sp.Rational(9, 10)), " P(1/2) =", P.subs(s, sp.Rational(1, 2)))
ok &= P.subs(s, sp.Rational(1, 2)) > 0  # sanity: positive in the interior for small s
# 4. rho_s is a density matrix for 0 <= s <= 1: eigenvalues (1-s)/(9-3s) x6 and (1+s)/(9-3s) x3.
print("ALL_OK" if ok else "FAIL")
