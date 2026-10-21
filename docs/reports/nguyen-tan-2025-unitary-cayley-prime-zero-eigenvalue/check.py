"""Nguyen-Tan, arXiv:2508.10348v1, Question 4.17: for a finite commutative Frobenius ring R and U <= R^x with -1 in U,
Gamma(R,U) connected and anti-connected: is Gamma(R,U) prime iff 0 is not an eigenvalue?
Checks: (1) fidelity, Proposition 4.18 on Z/n with U = R^x; (2) the counterexample R_d = F2[z]/(z^d - 1), U = <z>."""
import itertools, math
import numpy as np

def graph_from(elems, adj):
    n = len(elems); A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and adj(elems[i], elems[j]): A[i, j] = 1
    return A

def connected(A):
    n = len(A); seen = {0}; st = [0]
    while st:
        v = st.pop()
        for w in range(n):
            if A[v, w] and w not in seen: seen.add(w); st.append(w)
    return len(seen) == n

def complement(A):
    n = len(A); return np.ones((n, n), dtype=int) - np.eye(n, dtype=int) - A

def is_prime_bruteforce(A):
    n = len(A)
    # homogeneous set X: every vertex outside is adjacent to all or none of X; non-trivial: 2 <= |X| < n.
    # A non-trivial homogeneous set exists iff one exists that is minimal; check all pairs' closure (module closure).
    for a in range(n):
        for b in range(a + 1, n):
            X = {a, b}
            changed = True
            while changed:
                changed = False
                for v in range(n):
                    if v in X: continue
                    s = sum(A[v, x] for x in X)
                    if 0 < s < len(X):
                        X.add(v); changed = True
            if len(X) < n: return False
    return True

def zero_eig_mult(A):
    ev = np.linalg.eigvalsh(A.astype(float))
    return int(sum(abs(e) < 1e-9 for e in ev))

# (1) Proposition 4.18: R = Z/n, U = R^x, connected and anti-connected => prime iff 0 is not an eigenvalue.
print("Prop 4.18 on Z/n, U = (Z/n)^x:")
for n in range(2, 31):
    units = [u for u in range(n) if math.gcd(u, n) == 1]
    A = graph_from(list(range(n)), lambda a, b: (a - b) % n in units)
    if not (connected(A) and connected(complement(A))): continue
    p = is_prime_bruteforce(A); z = zero_eig_mult(A)
    print(f"  n={n}: prime={p} zero_eigenvalue_mult={z} agree={p == (z == 0)}")

# (2) R_d = F2[z]/(z^d-1): elements are d-bit vectors; multiplication is cyclic convolution mod 2; U = {z^i}.
def ring(d):
    elems = list(itertools.product([0, 1], repeat=d))
    def add(a, b): return tuple((x + y) % 2 for x, y in zip(a, b))
    def mul(a, b):
        c = [0] * d
        for i in range(d):
            if a[i]:
                for j in range(d):
                    if b[j]: c[(i + j) % d] ^= 1
        return tuple(c)
    return elems, add, mul

for d in (3, 4, 5, 6, 8):
    elems, add, mul = ring(d)
    one = tuple([1] + [0] * (d - 1))
    zpow = [tuple(1 if k == i else 0 for k in range(d)) for i in range(d)]
    # U is a subgroup of units: closed under multiplication, contains 1 and -1 = 1.
    closed = all(mul(a, b) in zpow for a in zpow for b in zpow)
    # Frobenius: psi = constant coefficient (F2-linear), and every nonzero a has r with psi(r a) = 1 (no nonzero ideal in ker psi).
    frob = all(any(mul(r, a)[0] == 1 for r in elems) for a in elems if any(a))
    A = graph_from(elems, lambda a, b: add(a, b) in zpow)
    deg = set(int(x) for x in A.sum(axis=1))
    con, acon = connected(A), connected(complement(A))
    p = is_prime_bruteforce(A)
    z = zero_eig_mult(A)
    print(f"R_{d}: |R|={len(elems)} U closed={closed} frobenius={frob} degree={deg} connected={con} anti-connected={acon} prime(bruteforce)={p} zero_eig_mult={z} binom(d,d/2)={math.comb(d, d // 2) if d % 2 == 0 else 0}")
