"""OEIS A388712 (G. Wiseman): T(n,k) = number of partitions of n into k parts whose
initial intervals {1..y_1}, ..., {1..y_k} are not choosable (no system of distinct
representatives). Conjecture: T(n,k) also counts partitions whose reversed
(weakly increasing) sequence (z_1, ..., z_k) is not superdiagonal (z_i >= i for all i).

The script decides choosability by brute force (search for distinct representatives)
for all partitions of n <= N and checks that it coincides with superdiagonality of the
reversed partition, and compares the row counts with the published triangle and the row
sums with A387118 (partitions of n without choosable initial intervals).
"""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 22
A387118 = [0,0,1,1,2,4,6,8,13,19,28,37,52,70,97,130,172,224,293,378,492,630,806]
DATA = [0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,2,1,1,0,0,0,2,2,1,1,0,0,0,1,3,2,1,1]

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for p in range(min(n, mx), 0, -1):
        for rest in partitions(n - p, p): yield (p,) + rest

def choosable(parts):
    used = set()
    def go(i):
        if i == len(parts): return True
        for x in range(1, parts[i] + 1):
            if x not in used:
                used.add(x)
                if go(i + 1): return True
                used.discard(x)
        return False
    return go(0)

def superdiagonal(z):
    return all(z[i] >= i + 1 for i in range(len(z)))

ok, tri, rowsums = True, [0], [0]  # row n = 0 is [T(0,0)] = [0]: the empty partition is choosable
for n in range(1, N + 1):
    row = [0] * (n + 1)  # columns k = 0..n
    for p in partitions(n):
        c = choosable(p)
        if c != superdiagonal(tuple(reversed(p))): ok = False
        if not c: row[len(p)] += 1
    tri.extend(row)
    rowsums.append(sum(row))
ok_data = tri[:len(DATA)] == DATA and rowsums[:len(A387118)] == A387118
print(f"N={N}: choosable <=> reversed superdiagonal on all partitions: {ok}; triangle (A388712) and row sums (A387118) match data: {ok_data}")
print("RESULT", "PASS" if ok and ok_data else "FAIL")
sys.exit(0 if ok and ok_data else 1)
