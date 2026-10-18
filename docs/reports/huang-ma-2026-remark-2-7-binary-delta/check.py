"""Huang-Ma arXiv:2602.03194v1 Remark 2.7: for fixed n and skew-symmetrizer D, delta(B) = det S(B) mod 4 is binary.

Checks the theorem of the preregistration: delta(B) lies in {0, 2} if n is odd or v_2(prod d_i) is odd,
and in {0, (-1)^(n/2) u mod 4} if n is even and prod d_i = 2^e u with e even, u odd.
 1. det S(B) = det A(B) for the integral companion A (A_ii = 2, A_ij = b_ij for i < j, -b_ij for i > j),
    in floating point with the radicals of S(B) (random matrices).
 2. Lemma 1: det A = (-1)^m det B mod 4 for n = 2m (random matrices).
 3. Exhaustive enumeration for n <= 4 with every pair parameter t in [-2, 2] (d_i b_ij = -d_j b_ji = t lcm(d_i, d_j)),
    for several symmetrizers; random sampling for n = 5, 6. Every residue set has at most two elements and lies in
    the predicted set.
"""
import itertools, math, random
from math import lcm
random.seed(20261018)

def det(M):
    # Bareiss fraction-free determinant over the integers
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]; sign = -sign; break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def build(n, d, ts):
    B = [[0] * n for _ in range(n)]
    for (i, j), t in zip(itertools.combinations(range(n), 2), ts):
        g = lcm(d[i], d[j]); B[i][j] = t * g // d[i]; B[j][i] = -t * g // d[j]
    return B

def companion(B):
    n = len(B); return [[2 if i == j else (B[i][j] if i < j else -B[i][j]) for j in range(n)] for i in range(n)]

def predicted(n, d):
    P = math.prod(d); e = 0
    while P % 2 == 0: P //= 2; e += 1
    if n % 2 == 1 or e % 2 == 1: return {0, 2}
    return {0, ((-1) ** (n // 2) * P) % 4}

def detS_float(B):
    n = len(B)
    S = [[2.0 if i == j else (math.copysign(1, B[i][j]) if B[i][j] else 0.0) * math.sqrt(abs(B[i][j] * B[j][i]))
          if i < j else ((math.copysign(1, B[j][i]) if B[j][i] else 0.0) * math.sqrt(abs(B[i][j] * B[j][i])) if i > j else 2.0)
          for j in range(n)] for i in range(n)]
    # Gaussian elimination with partial pivoting
    M = [r[:] for r in S]; dv = 1.0
    for k in range(n):
        p = max(range(k, n), key=lambda r: abs(M[r][k]))
        if abs(M[p][k]) < 1e-300: return 0.0
        if p != k: M[k], M[p] = M[p], M[k]; dv = -dv
        dv *= M[k][k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n): M[i][j] -= f * M[k][j]
    return dv

ok = True
worst = 0.0
for _ in range(400):
    n = random.randint(2, 6); d = [random.choice([1, 2, 3, 4, 5, 6, 8, 9, 12]) for _ in range(n)]
    B = build(n, d, [random.randint(-3, 3) for _ in range(n * (n - 1) // 2)])
    a = det(companion(B)); s = detS_float(B); worst = max(worst, abs(s - a) / max(1, abs(a)))
print('1. det S = det A: max relative difference', worst); ok &= worst < 1e-9
bad = 0
for _ in range(900):
    n = random.choice([2, 4, 6]); d = [random.choice([1, 2, 3, 4, 5, 6, 8, 12]) for _ in range(n)]
    B = build(n, d, [random.randint(-4, 4) for _ in range(n * (n - 1) // 2)])
    bad += (det(companion(B)) - (-1) ** (n // 2) * det(B)) % 4 != 0
print('2. Lemma 1 failures:', bad, 'of 900'); ok &= bad == 0
syms = {1: [(1,), (2,), (3,)], 2: [(1, 1), (1, 2), (2, 2), (1, 3), (3, 5), (2, 4), (4, 4), (2, 6)],
        3: [(1, 1, 1), (1, 1, 2), (1, 2, 4), (2, 3, 6), (1, 3, 9)],
        4: [(1, 1, 1, 1), (1, 1, 1, 2), (1, 1, 2, 2), (1, 1, 1, 3), (1, 1, 3, 3), (2, 2, 2, 2), (1, 2, 3, 6), (1, 1, 1, 5), (1, 1, 1, 7), (1, 2, 4, 8)]}
for n, ds in syms.items():
    for d in ds:
        vals = set(); count = 0
        for ts in itertools.product(range(-2, 3), repeat=n * (n - 1) // 2):
            vals.add(det(companion(build(n, list(d), ts))) % 4); count += 1
        good = len(vals) <= 2 and vals <= predicted(n, list(d))
        print(f'3. n={n} d={d}: {count} matrices, residues {sorted(vals)}, predicted {sorted(predicted(n, list(d)))}', 'ok' if good else 'FAIL')
        ok &= good
for n in (5, 6):
    for _ in range(6):
        d = [random.choice([1, 2, 3, 4, 6, 9]) for _ in range(n)]; vals = set()
        for _ in range(3000):
            vals.add(det(companion(build(n, d, [random.randint(-3, 3) for _ in range(n * (n - 1) // 2)]))) % 4)
        good = len(vals) <= 2 and vals <= predicted(n, d)
        print(f'3. n={n} d={tuple(d)}: 3000 random, residues {sorted(vals)}, predicted {sorted(predicted(n, d))}', 'ok' if good else 'FAIL')
        ok &= good
print('RESULT', 'PASS' if ok else 'FAIL')
