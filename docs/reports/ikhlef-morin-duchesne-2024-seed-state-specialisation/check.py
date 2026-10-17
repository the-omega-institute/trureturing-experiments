"""Ikhlef–Morin-Duchesne arXiv:2312.14837v2 (SciPost Phys. 17, 132 (2024)), Appendix C, Conjectures C.1 and C.2.

C.1 (p >= 1): X_p = 1, X_{p+1} = a_x + a_y a_w/{p+1},
  X_{k+1} = (a_x + {p} a_y a_w/({k}{k+1})) X_k + ({p}^2-{k}^2)(a_y^2-{k}^2)(a_w^2-{k}^2)/((q-1/q)^2 [2k-1][2k+1] {k}^2) X_{k-1};
  at y = eps q^m:  X_m = (w x)^{p-m} prod_{r=-(m-p-1)/2}^{(m-p-1)/2} (x q^{2r} + eps w)(eps q^{2r} + x w).
C.2 (p = 0): sigma_0 = 1, sigma_1 = a_x + a_y a_w/{1},
  sigma_{k+1} = (a_x + 2 a_y a_w/({k}{k+1})) sigma_k - [k]^2 (a_y^2-{k}^2)(a_w^2-{k}^2)/([2k+1][2k-1]{k}^2) sigma_{k-1};
  at y = eps q^m the same product with p = 0.
Exact rational arithmetic (fractions) at q in {2, 3, 3/2}, random rational x, w, eps = +-1, p <= 5, m <= p + 12.
Control: a perturbed recurrence coefficient breaks the identity.
"""
import random, sys
from fractions import Fraction as F

random.seed(20261017)


def br(q, k): return (q ** k - q ** (-k)) / (q - 1 / q)      # [k]
def cb(q, k): return q ** k + q ** (-k)                        # {k}
def al(z): return z + 1 / z


def X_C1(q, x, y, w, p, m, perturb=0):
    ax, ay, aw = al(x), al(y), al(w)
    if m == p: return F(1)
    prev, cur = F(1), ax + ay * aw / cb(q, p + 1)
    for k in range(p + 1, m):
        A = ax + cb(q, p) * ay * aw / (cb(q, k) * cb(q, k + 1))
        B = (cb(q, p) ** 2 - cb(q, k) ** 2) * (ay ** 2 - cb(q, k) ** 2) * (aw ** 2 - cb(q, k) ** 2) / (
            (q - 1 / q) ** 2 * br(q, 2 * k - 1) * br(q, 2 * k + 1) * cb(q, k) ** 2)
        prev, cur = cur, A * cur + (B + perturb) * prev
    return cur


def X_C2(q, x, y, w, m):
    ax, ay, aw = al(x), al(y), al(w)
    if m == 0: return F(1)
    prev, cur = F(1), ax + ay * aw / cb(q, 1)
    for k in range(1, m):
        A = ax + 2 * ay * aw / (cb(q, k) * cb(q, k + 1))
        B = br(q, k) ** 2 * (ay ** 2 - cb(q, k) ** 2) * (aw ** 2 - cb(q, k) ** 2) / (br(q, 2 * k + 1) * br(q, 2 * k - 1) * cb(q, k) ** 2)
        prev, cur = cur, A * cur - B * prev
    return cur


def product(q, x, w, eps, p, m):
    N = m - p; val = (w * x) ** (p - m)
    for j in range(N):
        s = -N + 1 + 2 * j          # s = 2r
        val *= (x * q ** s + eps * w) * (eps * q ** s + x * w)
    return val


def main():
    fails = []; n1 = n2 = 0
    for q in (F(2), F(3), F(3, 2)):
        for _ in range(4):
            x = F(random.randint(1, 9), random.randint(1, 9)); w = F(random.randint(1, 9), random.randint(1, 9))
            for eps in (1, -1):
                for p in range(1, 6):
                    for m in range(p, p + 13):
                        n1 += 1
                        if X_C1(q, x, eps * q ** m, w, p, m) != product(q, x, w, eps, p, m): fails.append(("C1", q, x, w, eps, p, m))
                for m in range(0, 13):
                    n2 += 1
                    if X_C2(q, x, eps * q ** m, w, m) != product(q, x, w, eps, 0, m): fails.append(("C2", q, x, w, eps, m))
                    # C.2 is the p -> 0 case of the C.1 recurrence ({0} = 2)
                    if X_C1(q, x, eps * q ** m, w, 0, m) != X_C2(q, x, eps * q ** m, w, m): fails.append(("C1p0", m))
    print(f"C.1 cases: {n1}; C.2 cases: {n2}; failures: {len(fails)}")
    # control: perturbing the off-diagonal coefficient breaks the identity for m >= p + 2
    q, x, w = F(2), F(3, 5), F(7, 4)
    broken = sum(1 for p in range(1, 4) for m in range(p + 2, p + 8)
                 if X_C1(q, x, q ** m, w, p, m, perturb=F(1, 10)) != product(q, x, w, 1, p, m))
    print(f"control (perturbed recurrence) breaks {broken}/18 cases (expected 18)")
    if broken != 18: fails.append("control")
    # the specialization matters: a generic y does not factor this way
    gen = X_C1(F(2), F(3, 5), F(5, 3), F(7, 4), 1, 5) != product(F(2), F(3, 5), F(7, 4), 1, 1, 5)
    print(f"control (generic y) differs from the product: {gen}")
    if not gen: fails.append("control-generic")
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:4]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
