"""A233155 family: n x k arrays over 0..2, horizontal (i,j)-(i,j+1) and antidiagonal (i,j+1)-(i+1,j)
bonds forbid x+y=2. Literal small counts, row transfer per width, b-file and recurrence checks."""
from itertools import product
from fractions import Fraction
import subprocess, json, sys

def ok(x, y): return x + y != 2

def rows(k): return [r for r in product(range(3), repeat=k) if all(ok(a, b) for a, b in zip(r, r[1:]))]

def compat(p, q):  # antidiagonal bonds between row p (upper) and row q (lower): (i,j+1)-(i+1,j)
    return all(ok(p[j + 1], q[j]) for j in range(len(p) - 1))

def literal(n, k):
    c = 0
    for f in product(range(3), repeat=n * k):
        g = [f[i * k:(i + 1) * k] for i in range(n)]
        if all(all(ok(a, b) for a, b in zip(r, r[1:])) for r in g) and all(compat(g[i], g[i + 1]) for i in range(n - 1)):
            c += 1
    return c

def seq(k, N):
    R = rows(k)
    d = {r: 1 for r in R}
    out = [sum(d.values())]
    for _ in range(N - 1):
        e = dict.fromkeys(R, 0)
        for p, m in d.items():
            if m:
                for q in R:
                    if compat(p, q): e[q] += m
        d = e
        out.append(sum(d.values()))
    return [0] + out, len(R)

def hrank(a, r):
    H = [[Fraction(a[i + j + 1]) for j in range(r)] for i in range(r)]
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

COLS = {5: ('A233152', [15, -18]), 6: ('A233153', [25, -90, 81]), 7: ('A233154', [42, -351, 972, -810])}
res = {}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
for k, (aid, coeff) in COLS.items():
    a, nrows = seq(k, N)
    lit = [literal(n, k) for n in (1, 2)]
    assert lit == a[1:3], (k, lit, a[1:3])
    num = aid[1:]
    b = subprocess.check_output(['curl', '-L', '--fail', '--silent', f'https://oeis.org/{aid}/b{num}.txt'], timeout=120).decode()
    b = {int(l.split()[0]): int(l.split()[1]) for l in b.splitlines() if l.strip() and not l.startswith('#')}
    chk = [n for n in b if n <= N]
    assert all(a[n] == b[n] for n in chk), k
    r = len(coeff)
    holds = [n for n in range(r + 1, N + 1) if a[n] == sum(c * a[n - i - 1] for i, c in enumerate(coeff))]
    res[aid] = {'width': k, 'legal_rows': nrows, 'literal_n1_n2': lit, 'initial': a[1:r + 3], 'bfile_checked': len(chk),
                'recurrence_holds_from': min(holds), 'fails_in_range': [n for n in range(r + 1, N + 1) if n not in holds],
                'hankel_rank_r_plus_2': hrank(a, r + 2)}
print(json.dumps(res, indent=1))
