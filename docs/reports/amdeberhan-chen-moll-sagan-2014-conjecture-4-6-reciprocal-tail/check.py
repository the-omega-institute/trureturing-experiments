"""Amdeberhan-Chen-Moll-Sagan, Ann. Comb. 18 (2014), Conjecture 4.6 (first clause): for integers s > t >= 1,
(s, t) != (2, 1), and n, r >= 1, floor((sum_{k>=n} 1/{rk}_{s,-t})^{-1}) = {rn}_{s,-t} - {r(n-1)}_{s,-t} - 1, where
{0} = 0, {1} = 1, {j} = s{j-1} - t{j-2}."""
import sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 260

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

def U(s, t, N):
    u = [0, 1]
    for j in range(2, N + 1): u.append(s * u[-1] - t * u[-2])
    return u

cases = bad = 0; bounds_bad = 0
for s in range(2, 13):
    for t in range(1, s):
        if (s, t) == (2, 1): continue
        u = U(s, t, 4000)
        if any(x <= 0 for x in u[1:]): bad += 1; continue
        for r in (1, 2, 3, 5):
            for n in (1, 2, 3, 4, 7):
                K = n
                tail = mp.mpf(0); k = n
                while r * k < 4000:
                    term = mp.mpf(1) / u[r * k]
                    tail += term; k += 1
                    if term < mp.mpf(10) ** (-240): break
                lhs = int(mp.floor(1 / tail)); rhs = u[r * n] - u[r * (n - 1)] - 1
                b = u[r * n] - u[r * (n - 1)]
                cases += 1
                if lhs != rhs: bad += 1; print("MISMATCH", s, t, r, n, lhs, rhs)
                if not (mp.mpf(1) / b < tail and (b == 1 or tail <= mp.mpf(1) / (b - 1))): bounds_bad += 1
check(f"floor identity holds in all {cases} cases (2 <= s <= 12, 1 <= t < s, r in {{1,2,3,5}}, n in {{1,2,3,4,7}})", bad == 0)
check(f"the proof's bounds 1/b_n < S_n <= 1/(b_n - 1) hold in all cases (failures {bounds_bad})", bounds_bad == 0)

# Cassini (Lemma 4) and the uniform domination H_k >= c^2 Q^(k-1) + 1 (Lemma 5), exact, on the stride sequences
lem_bad = 0
for s in range(2, 10):
    for t in range(1, s):
        if (s, t) == (2, 1): continue
        u = U(s, t, 400)
        for r in (1, 2, 3, 4):
            a = [u[r * k] for k in range(0, 400 // r)]
            c = a[1]; Q = t ** r
            for k in range(1, len(a) - 1):
                if a[k] ** 2 - a[k - 1] * a[k + 1] != c * c * Q ** (k - 1): lem_bad += 1
                if a[k + 1] - a[k - 1] < c * c * Q ** (k - 1) + 1: lem_bad += 1
check(f"Lemma 4 (Cassini a_k^2 - a_(k-1) a_(k+1) = c^2 Q^(k-1)) and Lemma 5 (H_k >= c^2 Q^(k-1) + 1), exact (failures {lem_bad})", lem_bad == 0)
# the excluded pair diverges: {n}_{2,-1} = n
check("(s, t) = (2, 1) gives {n} = n (harmonic tail, excluded)", U(2, 1, 10) == list(range(11)))
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
