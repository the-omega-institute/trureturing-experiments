"""A233172 (n x 6): literal enumeration for small n, transfer sequence, b-file and recurrence check."""
from itertools import product
from fractions import Fraction
import subprocess, json, sys

W = 6
beta = (0, 1, 2, 2, 1, 0)

def scan(s, c):
    x, y = s
    for a in c:
        if a in (2, 3, 4) and not x or a == 3 and not y:
            return None
        x |= a == 1
        y |= a == 2
    return x, y

def legal(p): return all(a != b for a, b in zip(p, p[1:]))
rows = [c for c in product(range(6), repeat=W) if legal(tuple(beta[v] for v in c))]

def compatible(p, q): return not any(p[j] == q[j + 1] or p[j + 1] == q[j] for j in range(W - 1))

# Literal depth-first count for n = 1, 2 straight from the definition.
def literal(n):
    def ok(arr):
        for i in range(n):
            for j in range(W):
                v = arr[i][j]
                for di, dj in ((0, 1), (1, 1), (1, -1)):
                    a, b = i + di, j + dj
                    if 0 <= a < n and 0 <= b < W:
                        u = arr[a][b]
                        if u == v or u == 5 - v: return False
        return True
    cnt = 0
    for arr in product(rows, repeat=n):
        if arr[0][0] != 0 or scan((False, False), [v for r in arr for v in r]) is None: continue
        if ok(arr): cnt += 1
    return cnt

# Transfer on (scan state, label row).
from collections import defaultdict
first = defaultdict(int)
for c in rows:
    if c[0] != 0: continue
    t = scan((False, False), c)
    if t is not None: first[(t, tuple(beta[v] for v in c))] += 1
trans = {}
def step(dist):
    out = defaultdict(int)
    for (s, p), m in dist.items():
        key = (s, p)
        if key not in trans:
            d = defaultdict(int)
            for c in rows:
                q = tuple(beta[v] for v in c)
                if not compatible(p, q): continue
                t = scan(s, c)
                if t is not None: d[(t, q)] += 1
            trans[key] = d
        for k, v in trans[key].items(): out[k] += m * v
    return out
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
D = [dict(first)]
for _ in range(N - 1): D.append(step(D[-1]))
a = [0] + [sum(d.values()) for d in D]
lit = [literal(1), literal(2)]
assert lit == a[1:3], (lit, a[1:3])
bfile = subprocess.check_output(['curl', '-L', '--fail', '--silent', 'https://oeis.org/A233172/b233172.txt'], timeout=60).decode()
b = {int(l.split()[0]): int(l.split()[1]) for l in bfile.splitlines() if l.strip() and not l.startswith('#')}
checked = [n for n in b if n <= N]
assert all(a[n] == b[n] for n in checked)
coeff = [1216, -303104, -6553600, 6056574976, -287762808832, -18622978195456, 1482141674242048,
         -15199648742375424, -522417556774977536, 5764607523034234880]
r = len(coeff)
hold = [n for n in range(r + 1, N + 1) if a[n] == sum(c * a[n - i - 1] for i, c in enumerate(coeff))]
first_fail = [n for n in range(r + 1, N + 1) if n not in hold]
# Hankel rank of the sequence.
H = [[Fraction(a[i + j + 1]) for j in range(r + 2)] for i in range(r + 2)]
rk = 0
for j in range(r + 2):
    piv = next((i for i in range(rk, r + 2) if H[i][j]), None)
    if piv is None: continue
    H[rk], H[piv] = H[piv], H[rk]
    z = H[rk][j]; H[rk] = [x / z for x in H[rk]]
    for i in range(r + 2):
        if i != rk:
            z = H[i][j]; H[i] = [x - z * y for x, y in zip(H[i], H[rk])]
    rk += 1
res = {'literal_n1_n2': lit, 'initial_values': a[1:r + 3], 'states_reached': len(trans),
       'bfile_terms_checked': len(checked), 'recurrence_fail_indices': first_fail,
       'hankel_rank_order_r_plus_2': rk, 'N': N}
print(json.dumps(res, indent=1))
