"""Checks for Lee-Liu, arXiv:2502.00134, Conjecture 3.7: for n >= 49 the integer {2}-domination number of the
4 x n grid P_4 x P_n is 2n - floor((n-9)/10) (functions f: V -> {0,1,2}, closed-neighbourhood sums >= 2).
(1) Residual automaton: states (last column b, residual demand r on the next column); E_n(s) is the minimum weight of
    the first n columns that satisfy every constraint of columns 1..n-1 and end in s.  Checks 389 states, 15412
    transitions, every coordinate of E_67 finite, E_77 = E_67 + 19 coordinatewise, and gamma(n) = formula for 49..76.
    With min-plus homogeneity F(x + k) = F(x) + k this gives gamma(n + 10) = gamma(n) + 19 for n >= 67, hence all n >= 49.
(2) An independent uncompressed two-column dynamic programme gives the same gamma(n) for 1 <= n <= 120, and gamma(48) = 92
    while the formula gives 93 (the threshold n >= 49 is needed).
(3) The ten-column cyclic pattern of weight 19 satisfies every local constraint (interior density 19/10).  Prints ALL_OK."""
from itertools import product
C = list(product(range(3), repeat=4))
z = (0, 0, 0, 0)
def V(b):
    return tuple((b[i-1] if i else 0) + b[i] + (b[i+1] if i < 3 else 0) for i in range(4))
def R(a, b):
    v = V(b)
    return tuple(max(0, 2-a[i]-v[i]) for i in range(4))
S = sorted({(b, R(a, b)) for a in C for b in C})
ix = {s: i for i, s in enumerate(S)}
edges = [[(ix[(c, R(b, c))], sum(c)) for c in C if all(c[i] >= r[i] for i in range(4))] for b, r in S]
E = [None] * len(S)
for b in C:
    E[ix[(b, R(z, b))]] = sum(b)
terminal = [i for i, (b, r) in enumerate(S) if r == z]
seen, gamma = {}, {}
for n in range(1, 78):
    seen[n] = tuple(E)
    gamma[n] = min(E[i] for i in terminal if E[i] is not None)
    nxt = [None] * len(S)
    for i, old in enumerate(E):
        if old is not None:
            for j, cost in edges[i]:
                value = old + cost
                if nxt[j] is None or value < nxt[j]:
                    nxt[j] = value
    E = nxt
q = lambda n: 2*n - (n-9)//10
ok = len(S) == 389 and sum(map(len, edges)) == 15412
ok &= all(x is not None for x in seen[67]) and all(y == x+19 for x, y in zip(seen[67], seen[77]))
ok &= all(gamma[n] == q(n) for n in range(49, 77))
print("residual automaton: states", len(S), "transitions", sum(map(len, edges)), "E77 = E67 + 19 and gamma(49..76) = formula:", ok)
# (2) uncompressed DP over (previous column, current column)
def A(a, b, c):
    v = V(b)
    return all(a[i] + v[i] + c[i] >= 2 for i in range(4))
need = {}
succ = {}
for a, b in product(C, repeat=2):
    succ[(a, b)] = [c for c in C if A(a, b, c)]
INF = 10**9
dp = {(z, b): sum(b) for b in C}
g2 = {1: min(sum(b) for b in C if A(z, b, z))}
for n in range(2, 121):
    nd = {}
    for (a, b), w in dp.items():
        for c in succ[(a, b)]:
            k = (b, c); v = w + sum(c)
            if v < nd.get(k, INF): nd[k] = v
    dp = nd
    g2[n] = min(w for (a, b), w in dp.items() if A(a, b, z))
ok2 = all(g2[n] == gamma[n] for n in range(1, 78)) and all(g2[n] == q(n) for n in range(49, 121)) and g2[48] == 92 and q(48) == 93
print("uncompressed DP agrees for n <= 77, formula holds for 49..120, gamma(48) = 92 < 93:", ok2)
# (3) ten-column cyclic pattern, digits top to bottom
P = [tuple(int(x) for x in s) for s in "1011 0010 1100 0002 1100 0010 1011 1000 0111 1000".split()]
ok3 = sum(map(sum, P)) == 19 and all(A(P[j-1], P[j], P[(j+1) % 10]) for j in range(10))
print("ten-column cyclic pattern of weight 19 is locally {2}-dominating:", ok3)
print("ALL_OK" if ok and ok2 and ok3 else "FAIL")
