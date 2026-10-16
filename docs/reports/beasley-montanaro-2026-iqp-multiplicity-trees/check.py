"""Beasley–Montanaro arXiv:2610.19962, Conjecture 6.5 (IQP circuit distinguishability): candidate 11-vertex tree pair.

Checks (all exact or with an explicit tolerance):
1. A and B are trees on 11 vertices, and non-isomorphic (AHU canonical forms over all roots).
2. Direct state-vector simulation of |IQP_G> = prod_{(i,j) in E} exp(i theta X_i X_j) |0^n> (Definition 2.4)
   at 41 angles in [0, pi/2] (endpoints included): the multisets {(|y|, P(y))} coincide for A and B
   (probabilities grouped at 1e-12), i.e. M^k_A(p) = M^k_B(p) for all k, p (Definition 6.2).
3. The exact coefficient tables c_T[j,k] = #{S subset E : |S| = j, |odd(S)| = k} coincide, and
   S -> odd(S) is injective on each tree (so P(y) = sin^{2j} cos^{2(m-j)} with j = |S(y)|).
4. Controls: a non-isomorphic pair that the invariant does separate (path P11 vs. A), and a one-edge
   perturbation of B, non-isomorphic to A and B, whose table and multiplicities differ.
"""
import itertools, math, sys
import numpy as np

A = [(0, 1), (0, 2), (0, 3), (3, 4), (0, 5), (5, 6), (5, 7), (7, 8), (8, 9), (8, 10)]
B = [(0, 1), (0, 2), (0, 3), (3, 4), (3, 5), (0, 6), (6, 7), (7, 8), (7, 9), (9, 10)]
N = 11


def adj(E):
    g = {v: [] for v in range(N)}
    for u, v in E:
        g[u].append(v); g[v].append(u)
    return g


def is_tree(E):
    if len(E) != N - 1: return False
    g = adj(E); seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in g[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == N


def ahu(g, root):
    def enc(x, p):
        return "(" + "".join(sorted(enc(y, x) for y in g[x] if y != p)) + ")"
    return enc(root, -1)


def canon(E):
    g = adj(E)
    return min(ahu(g, r) for r in range(N))


def table(E):
    m = len(E); t = {}
    seen = set()
    injective = True
    for s in range(1 << m):
        bd = 0
        for e in range(m):
            if s >> e & 1:
                u, v = E[e]; bd ^= (1 << u) | (1 << v)
        if bd in seen: injective = False
        seen.add(bd)
        key = (bin(s).count("1"), bin(bd).count("1"))
        t[key] = t.get(key, 0) + 1
    return t, injective


def state(E, th):
    psi = np.zeros(1 << N, complex); psi[0] = 1
    c, s = math.cos(th), math.sin(th)
    idx = np.arange(1 << N)
    for u, v in E:
        flip = idx ^ ((1 << u) | (1 << v))
        psi = c * psi + 1j * s * psi[flip]
    return psi


def multiset(E, th):
    P = np.abs(state(E, th)) ** 2
    w = np.array([bin(y).count("1") for y in range(1 << N)])
    return sorted(zip(w.tolist(), np.round(P, 12).tolist()))


def main():
    fails = []
    ok_tree = is_tree(A) and is_tree(B)
    noniso = canon(A) != canon(B)
    print("both trees:", ok_tree, "| non-isomorphic (AHU over all roots):", noniso)
    if not (ok_tree and noniso): fails.append("trees")
    deg = lambda E: sorted(len(v) for v in adj(E).values())
    print("degree sequences:", deg(A), deg(B))
    tA, injA = table(A); tB, injB = table(B)
    print("boundary map injective:", injA, injB, "| tables equal:", tA == tB)
    if not (injA and injB and tA == tB): fails.append("tables")
    angles = [math.pi / 2 * k / 40 for k in range(41)]
    bad = [th for th in angles if multiset(A, th) != multiset(B, th)]
    print(f"state-vector multiplicities equal at {len(angles) - len(bad)}/{len(angles)} angles")
    if bad: fails.append(("multiplicity", bad[:3]))
    rng = np.random.default_rng(20261017)
    rand = [float(x) for x in rng.uniform(0, math.pi / 2, 20)]
    bad2 = [th for th in rand if multiset(A, th) != multiset(B, th)]
    print(f"random angles equal: {20 - len(bad2)}/20")
    if bad2: fails.append(("multiplicity-random", bad2[:3]))
    # norm check of the simulation
    nrm = max(abs(np.linalg.norm(state(A, th)) - 1) for th in angles)
    print("max |norm-1|:", f"{nrm:.1e}")
    # controls
    P11 = [(i, i + 1) for i in range(10)]
    sep = any(multiset(P11, th) != multiset(A, th) for th in angles)
    print("control P11 vs A separated:", sep)
    if not sep: fails.append("control-path")
    Bp = [e for e in B if e != (0, 1)] + [(1, 2)]
    tBp, _ = table(Bp)
    fresh = is_tree(Bp) and canon(Bp) not in (canon(A), canon(B))
    sepBp = any(multiset(Bp, th) != multiset(A, th) for th in angles)
    print("control B - (0,1) + (1,2): tree, non-isomorphic to A and B:", fresh,
          "| table differs:", tBp != tA, "| separated by multiplicities:", sepBp)
    if not (fresh and tBp != tA and sepBp): fails.append("control-perturb")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
