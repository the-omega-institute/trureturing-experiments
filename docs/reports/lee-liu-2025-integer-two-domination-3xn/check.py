"""Checks for Lee-Liu, arXiv:2502.00134, Conjecture 2.5: the integer {2}-domination number of the 3 x n grid
P_3 x P_n is floor(3n/2) + 1 for every n >= 1 (functions f: V -> {0,1,2}, closed-neighbourhood sums >= 2).
(1) The potential certificate (certificate.json, the r333q search seat's table) satisfies the three finite lemmas:
    Phi(z,b) <= 2s(b)-3; A(a,b,c) => Phi(b,c) <= Phi(a,b)+2s(c)-3; A(a,b,z) => Phi(a,b) >= 1, Phi = 2s(a)+2s(b)+R.
(2) An exact column-transfer dynamic programme gives the true minimum for 1 <= n <= 40 and it equals floor(3n/2)+1.
(3) The explicit construction attains the bound for 1 <= n <= 40.  Prints ALL_OK."""
import json, itertools, os
here = os.path.dirname(os.path.abspath(__file__))
cert = json.load(open(os.path.join(here, 'certificate.json')))
C = list(itertools.product(range(3), repeat=3))
keys = [''.join(map(str, a)) for a in C]
ok = cert['column_order'] == keys
rows = {a: list(v) for a, v in cert['explicit_rows'].items()}
for value, group in cert['constant_rows'].items():
    for a in group:
        ok &= a not in rows; rows[a] = [int(value)]*27
ok &= set(rows) == set(keys)
R = {(a, b): rows[keys[i]][j] for i, a in enumerate(C) for j, b in enumerate(C)}
z = (0, 0, 0)
def A(a, b, c):
    u = (b[0]+b[1], sum(b), b[1]+b[2])
    return all(a[i]+u[i]+c[i] >= 2 for i in range(3))
phi = lambda a, b: 2*sum(a) + 2*sum(b) + R[a, b]
ok &= all(phi(z, b) <= 2*sum(b) - 3 for b in C)
triples = [(a, b, c) for a, b, c in itertools.product(C, repeat=3) if A(a, b, c)]
ok &= len(triples) == 17138 and all(phi(b, c) <= phi(a, b) + 2*sum(c) - 3 for a, b, c in triples)
terminal = [(a, b) for a, b in itertools.product(C, repeat=2) if A(a, b, z)]
ok &= len(terminal) == 533 and all(phi(a, b) >= 1 for a, b in terminal)
print("certificate lemmas L1-L3 hold:", ok, " admissible triples", len(triples), " terminal pairs", len(terminal))
# exact DP over states (previous column, current column)
INF = 10**9
def gamma(n):
    if n == 1:
        return min(sum(b) for b in C if A(z, b, z))
    dp = {(z, b): sum(b) for b in C}           # columns x0 = z, x1 = b
    for j in range(2, n+1):
        nd = {}
        for (a, b), w in dp.items():
            for c in C:
                if A(a, b, c):
                    k = (b, c); v = w + sum(c)
                    if v < nd.get(k, INF): nd[k] = v
        dp = nd
    return min(w for (a, b), w in dp.items() if A(a, b, z))
for n in range(1, 41):
    g = gamma(n)
    ok &= g == 3*n//2 + 1
print("DP: gamma_{2}(G_{3,n}) = floor(3n/2)+1 for n = 1..40:", ok)
def construction(n):
    if n == 1: return [(0, 2, 0)]
    cols = []
    for j in range(1, n+1):
        if j % 2 == 0: cols.append((1, 1, 1))
        elif j == 1 or j == n: cols.append((0, 1, 0))
        else: cols.append((0, 0, 0))
    return cols
for n in range(1, 41):
    cols = construction(n); ext = [z] + cols + [z]
    ok &= all(A(ext[j-1], ext[j], ext[j+1]) for j in range(1, n+1)) and sum(map(sum, cols)) == 3*n//2 + 1
print("construction attains the bound for n = 1..40:", ok)
print("ALL_OK" if ok else "FAIL")
