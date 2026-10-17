"""da Silva arXiv:2305.05795v1, Section 1 / Section 3: for which dimensions does the tensor product preserve extremality
of UCPT maps? Claimed answer: for all extreme points iff min(n, m) <= 2.
Exact (sympy) checks:
1. The source's qutrit extreme map (E0..E3, eqs. (9)-(12)) and its padded version Phi_d (K0 = E0 + I_{d-3}) for
   d = 3..6: unital and trace preserving; Kraus operators linearly independent; the 16 vectors K_i^dag K_j (+) K_j K_i^dag
   linearly independent (Landau–Streater criterion, the source's Theorem 2.18) => extreme.
2. For (n, m) in {3,4,5} x {3,4,5}: K_0^dag K_1 = 0 and K_2 K_0^dag = 0, so for L_a = K_0 (x) K_0, L_b = K_1 (x) K_2 both
   L_a^dag L_b and L_b L_a^dag vanish; the explicit maps Theta_pm with coefficient matrices I +- H/2 are UCPT, average to
   Phi_n (x) Phi_m and are distinct => not extreme.
3. For n m >= 12 the Choi rank 16 of Phi_n (x) Phi_m satisfies 16 <= sqrt(2) n m, so the source's rank bound
   (Theorem 3.5) does not apply; the zero-product certificate is what proves non-extremality.
"""
import itertools, sys
import sympy as sp

s2, s3 = sp.sqrt(2), sp.sqrt(3)


def E(i, j, d):
    M = sp.zeros(d, d); M[i - 1, j - 1] = 1; return M


def kraus(d):
    K0 = E(1, 1, d) / 2
    for k in range(4, d + 1): K0 += E(k, k, d)
    K1 = (E(2, 1, d) + s2 * E(3, 2, d)) / 2
    K2 = (s2 * E(1, 2, d) + s3 * E(2, 3, d)) / 2
    K3 = (E(1, 3, d) + s2 * E(3, 1, d)) / 2
    return [K0, K1, K2, K3]


def vec(M): return list(M)


def is_ucpt(Ks, d):
    return sp.simplify(sum((K.H * K for K in Ks), sp.zeros(d, d)) - sp.eye(d)) == sp.zeros(d, d) and \
           sp.simplify(sum((K * K.H for K in Ks), sp.zeros(d, d)) - sp.eye(d)) == sp.zeros(d, d)


def rank_of(vectors):
    return sp.Matrix(vectors).rank(simplify=True)


def main():
    fails = []
    for d in range(3, 7):
        Ks = kraus(d)
        u = is_ucpt(Ks, d)
        rk = rank_of([vec(K) for K in Ks])
        joint = [vec(Ki.H * Kj) + vec(Kj * Ki.H) for Ki in Ks for Kj in Ks]
        rj = rank_of(joint)
        print(f"Phi_{d}: UCPT {u}; Kraus rank {rk}/4; joint products rank {rj}/16 -> extreme {u and rk == 4 and rj == 16}")
        if not (u and rk == 4 and rj == 16): fails.append(("factor", d))
    for n, m in itertools.product((3, 4, 5), repeat=2):
        Kn, Km = kraus(n), kraus(m)
        z1 = Kn[0].H * Kn[1] == sp.zeros(n, n); z2 = Km[2] * Km[0].H == sp.zeros(m, m)
        L = [sp.kronecker_product(Kn[i], Km[j]) for i in range(4) for j in range(4)]
        a, b = 0, 1 * 4 + 2
        za = L[a].H * L[b] == sp.zeros(n * m, n * m) and L[b] * L[a].H == sp.zeros(n * m, n * m)
        ok = z1 and z2 and za
        if (n, m) in ((3, 3), (3, 5), (5, 4)):
            N = n * m
            H = sp.zeros(16, 16); H[a, b] = 1; H[b, a] = 1
            maps = []
            for sgn in (1, -1):
                R = sp.eye(16) + sgn * H / 2
                maps.append(R)
                pos = all(ev > 0 for ev in R.eigenvals())
                Kraus_sum_unit = sp.simplify(sum((R[al, be] * L[be].H * L[al] for al in range(16) for be in range(16)), sp.zeros(N, N)) - sp.eye(N)) == sp.zeros(N, N)
                Kraus_sum_tp = sp.simplify(sum((R[al, be] * L[al] * L[be].H for al in range(16) for be in range(16)), sp.zeros(N, N)) - sp.eye(N)) == sp.zeros(N, N)
                ok = ok and pos and Kraus_sum_unit and Kraus_sum_tp
            # distinct: (Theta_+ - Theta_-)(X) for X = |1,1><1,3| (indices (1,1) and (1,3) in n x m)
            X = sp.zeros(N, N); X[0, 2] = 1
            diff = sp.zeros(N, N)
            for al in range(16):
                for be in range(16):
                    c = maps[0][al, be] - maps[1][al, be]
                    if c != 0: diff += c * L[al] * X * L[be].H
            nz = sp.simplify(diff) != sp.zeros(N, N)
            ok = ok and nz
            print(f"({n},{m}): zero products {z1 and z2 and za}; Theta_pm positive and UCPT; Theta_+ != Theta_-: {nz}")
        if not ok: fails.append(("tensor", n, m))
    print("tensor pairs (n, m) in {3,4,5}^2: zero-product certificate holds for all" if not any(f[0] == "tensor" for f in fails) else "tensor FAIL")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
