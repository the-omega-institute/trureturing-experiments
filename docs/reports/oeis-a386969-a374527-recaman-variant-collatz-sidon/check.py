"""OEIS A386969 and A374527.

A386969 (J. Beauchamp): a(0) = 0; for n > 0, a(n) = a(n-1) - n if nonnegative and not already in the
sequence, otherwise a(n) = a(n-1) + n + 3. Conjectures:
- (Z. He) g.f. x*(4+2x+11x^3-2x^4-9x^5)/(1-x^2)^2;
- (B. McEachen) a(n) = n + 5 for odd n > 2 and a(n) = 2n + 7 for even n > 3.
Proof checked here: by induction, for n >= 4 the set of earlier values is {0, 2, 4} together with
{m + 5 : m odd, 3 <= m < n} and {2m + 7 : m even, 4 <= m < n}; at even n the candidate a(n-1) - n = 4 is taken,
and at odd n the candidate n + 5 is even and larger than every earlier even value.

A374527 (M. Sigg, H. Pfoertner): odd numbers whose Collatz trajectory is a Sidon set. Conjecture (R. Stephan):
it contains every (4^k - 1)/3 except 5. Proof checked here: for x = (4^k - 1)/3 (k >= 1) the trajectory is
x, 4^k, 2^(2k-1), ..., 2, 1; sums of two powers of two have binary weight <= 2, while x + 2^a and 2x have
weight >= 3 for k >= 3; k = 2 fails since 5 + 5 = 8 + 2.
"""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
K = int(sys.argv[2]) if len(sys.argv) > 2 else 60
seq = [0]; seen = {0}
for n in range(1, N + 1):
    p = seq[-1]
    v = p - n if p - n >= 0 and (p - n) not in seen else p + n + 3
    seq.append(v); seen.add(v)
DATA = [0,4,2,8,15,10,19,12,23,14,27,16,31,18,35,20,39,22,43,24,47,26,51,28,55,30,59,32,63,34,67,36,71,38,75,40,79,42,83,44,87]
ok_data = seq[:len(DATA)] == DATA
mce = all(seq[n] == (n + 5 if n % 2 else 2 * n + 7) for n in range(3, N + 1) if n % 2 or n > 3)
num = [0, 4, 2, 0, 11, -2, -9]
gf = all(sum(c * (seq[n - i] if n - i >= 0 else 0) for i, c in enumerate([1, 0, -2, 0, 1])) == (num[n] if n < len(num) else 0) for n in range(0, N + 1))
inv = True
for n in range(4, N + 1):
    want = {0, 2, 4} | {m + 5 for m in range(3, n, 2)} | {2 * m + 7 for m in range(4, n, 2)}
    if set(seq[:n]) != want: inv = False; break
print(f"A386969: data {ok_data}; McEachen closed forms n<={N}: {mce}; He g.f. (1-x^2)^2 A(x) = numerator, n<={N}: {gf}; set invariant n<={N}: {inv}")
def collatz(x):
    t = [x]
    while x != 1:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        t.append(x)
    return t
def sidon(t):
    s = set()
    for i in range(len(t)):
        for j in range(i, len(t)):
            v = t[i] + t[j]
            if v in s: return False
            s.add(v)
    return True
res = {k: sidon(collatz((4 ** k - 1) // 3)) for k in range(1, K + 1)}
ok_sidon = all(res[k] == (k != 2) for k in res)
traj = all(collatz((4 ** k - 1) // 3) == [(4 ** k - 1) // 3] + [2 ** e for e in range(2 * k, -1, -1)] for k in range(2, K + 1))
early = [x for x in range(1, 400, 2) if sidon(collatz(x))]
print(f"A374527: Sidon for (4^k-1)/3 iff k != 2, k<={K}: {ok_sidon}; trajectory shape: {traj}; odd terms < 400: {early}")
good = ok_data and mce and gf and inv and ok_sidon and traj
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
