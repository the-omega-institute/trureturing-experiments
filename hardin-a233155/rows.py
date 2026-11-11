"""A233155 rows: R x n arrays over 0..2 (R = 3..6), horizontal (i,j)-(i,j+1) and antidiagonal (i,j+1)-(i+1,j) bonds
forbid x+y=2. Column transfer: column c (top to bottom) followed by column d is legal when, for every row i,
c[i]+d[i] != 2 (horizontal) and, for every i < R-1, d[i]+c[i+1] != 2 (antidiagonal (i,j+1)-(i+1,j))."""
from itertools import product
from fractions import Fraction
import subprocess, json, sys

def ok(x, y): return x + y != 2

def seq(R, N):
    cols = list(product(range(3), repeat=R))  # no vertical bond: every column is legal
    nxt = {c: [d for d in cols if all(ok(c[i], d[i]) for i in range(R)) and all(ok(d[i], c[i + 1]) for i in range(R - 1))] for c in cols}
    v = {c: 1 for c in cols}
    out = [sum(v.values())]
    for _ in range(N - 1):
        w = dict.fromkeys(cols, 0)
        for c, m in v.items():
            for d in nxt[c]: w[d] += m
        v = w
        out.append(sum(v.values()))
    return [0] + out

def literal(R, n):
    c = 0
    for f in product(range(3), repeat=R * n):
        g = [f[i * n:(i + 1) * n] for i in range(R)]
        if all(ok(g[i][j], g[i][j + 1]) for i in range(R) for j in range(n - 1)) and \
           all(ok(g[i][j + 1], g[i + 1][j]) for i in range(R - 1) for j in range(n - 1)):
            c += 1
    return c

def rank(a, r, s):
    H = [[Fraction(a[s + i + j]) for j in range(r)] for i in range(r)]
    rk = 0
    for j in range(r):
        piv = next((i for i in range(rk, r) if H[i][j]), None)
        if piv is None: continue
        H[rk], H[piv] = H[piv], H[rk]
        z = H[rk][j]; H[rk] = [x / z for x in H[rk]]
        for i in range(r):
            if i != rk:
                z = H[i][j]; H[i] = [x - z * y for x, y in zip(H[i], H[rk])]
        rk += 1
    return rk

ROWS = {3: ('A233156', [5, -2], 4), 4: ('A233157', [9, -15, 6], 7),
        5: ('A233158', [16, -64, 76, 17, -84, 48, -8], 11),
        6: ('A233159', [30, -283, 1038, -573, -6705, 21056, -27636, 17216, -4128], 15)}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 80
res = {}
for R, (aid, coeff, after) in ROWS.items():
    a = seq(R, N)
    lit = [literal(R, n) for n in (1, 2)]
    assert lit == a[1:3], (R, lit, a[1:3])
    b = subprocess.check_output(['curl', '-L', '--fail', '--silent', f'https://oeis.org/{aid}/b{aid[1:]}.txt'], timeout=120).decode()
    b = {int(l.split()[0]): int(l.split()[1]) for l in b.splitlines() if l.strip() and not l.startswith('#')}
    chk = [n for n in b if n <= N]
    assert all(a[n] == b[n] for n in chk), R
    r = len(coeff)
    holds = [n for n in range(r + 1, N + 1) if a[n] == sum(c * a[n - i - 1] for i, c in enumerate(coeff))]
    res[aid] = {'rows': R, 'columns': 3 ** R, 'literal_n1_n2': lit, 'bfile_checked': len(chk), 'source_start': after + 1,
                'holds_from': min(holds), 'fails_at_or_after_source_start': [n for n in range(after + 1, N + 1) if n not in holds],
                'last_failure_before': max([n for n in range(r + 1, N + 1) if n not in holds], default=None),
                'hankel_rank_tail': rank(a, r + 2, after + 1)}
print(json.dumps(res, indent=1))
