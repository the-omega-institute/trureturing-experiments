"""Ramane–Kitturmath, Examples and Counterexamples 3 (2023) 100115, Conjecture 3.11:
"There is no regular graph with adjacency energy 2k√2, where k∈Z⁺."

Checks:
1. Every r-regular graph on n <= 8 vertices (all labeled graphs, deduplicated by spectrum) and every circulant
   graph on n <= 26 vertices: E_A(G) differs from 2k√2 for every k >= 1, and more strongly, whenever E_A(G)^2 is
   (numerically) an integer, E_A(G) is an even integer (the Galois argument plus Bapat–Pati).
2. Lemma 2 of the proof route on every spectrum: every submultiset sum lies in [-S, S]; with the eigenvalue r
   removed, every submultiset sum is >= -S, so r - S is a lower bound for every sum containing r.
3. Controls: the path P3 and the stars K_{1,2k^2} (not regular) have energy exactly 2k√2 (k = 1 for P3), so the
   regularity hypothesis is needed; C_8 has energy 4 + 4√2.
"""
import itertools, math, sys
import numpy as np

TOL = 1e-8


def energy(A):
    w = np.linalg.eigvalsh(A)
    return float(np.abs(w).sum()), w


def regular_graphs(n):
    """All labeled graphs on n vertices that are regular, one per spectrum."""
    pairs = list(itertools.combinations(range(n), 2))
    seen = {}
    # backtracking over vertices: choose the neighbourhood of each vertex among later vertices
    for r in range(0, n):
        if (n * r) % 2: continue
        deg = [0] * n; adj = np.zeros((n, n))

        def rec(v):
            if v == n:
                if all(d == r for d in deg):
                    E, w = energy(adj)
                    key = (r, tuple(np.round(np.sort(w), 6)))
                    seen.setdefault(key, adj.copy())
                return
            need = r - deg[v]
            later = [u for u in range(v + 1, n) if deg[u] < r]
            if need < 0 or need > len(later): return
            for S in itertools.combinations(later, need):
                for u in S: adj[v, u] = adj[u, v] = 1; deg[u] += 1
                deg[v] += need
                rec(v + 1)
                deg[v] -= need
                for u in S: adj[v, u] = adj[u, v] = 0; deg[u] -= 1
        rec(0)
    return seen


def circulants(n):
    out = {}
    half = list(range(1, n // 2 + 1))
    for m in range(0, len(half) + 1):
        for C in itertools.combinations(half, m):
            A = np.zeros((n, n))
            for i in range(n):
                for c in C:
                    A[i, (i + c) % n] = A[(i + c) % n, i] = 1
            r = int(A[0].sum())
            E, w = energy(A)
            out.setdefault((r, tuple(np.round(np.sort(w), 6))), A)
    return out


def check_spectrum(r, w):
    E = float(np.abs(w).sum()); S = E / 2
    bad = []
    k = round(E / (2 * math.sqrt(2)))
    if k >= 1 and abs(E - 2 * k * math.sqrt(2)) < TOL: bad.append(("equals 2k sqrt2", E))
    E2 = E * E
    if abs(E2 - round(E2)) < 1e-7 and r >= 1:
        if abs(E - round(E)) > 1e-7 or round(E) % 2: bad.append(("E^2 integer but E not even integer", E))
    # Lemma 2: remove one copy of r, every submultiset sum of the rest is >= -S
    rest = list(w); idx = int(np.argmin(np.abs(np.array(rest) - r))); assert abs(rest[idx] - r) < 1e-7; rest.pop(idx)
    if r >= 1 and sum(x for x in rest if x < 0) < -S - 1e-9: bad.append(("lemma2", E))
    return bad, E


def main():
    fails = []; count = 0; E2int = 0
    for n in range(1, 9):
        gs = regular_graphs(n)
        for (r, _), A in gs.items():
            count += 1
            b, E = check_spectrum(r, np.linalg.eigvalsh(A)); fails += [(n, r) + x for x in b]
            if abs(E * E - round(E * E)) < 1e-7: E2int += 1
        print(f"n={n}: {len(gs)} regular spectra")
    nc = 0
    for n in range(3, 27):
        cs = circulants(n); nc += len(cs)
        for (r, _), A in cs.items():
            b, E = check_spectrum(r, np.linalg.eigvalsh(A)); fails += [("circ", n, r) + x for x in b]
            if abs(E * E - round(E * E)) < 1e-7: E2int += 1
    print(f"exhaustive regular spectra (n<=8): {count}; circulant spectra (n<=26): {nc}; with E^2 integral: {E2int}")
    # controls: regularity is needed. P3 is not regular and has energy 2*sqrt(2).
    P3 = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], float)
    e3, _ = energy(P3)
    print(f"control P3 (not regular): energy = {e3:.12f}, 2*sqrt(2) = {2*math.sqrt(2):.12f}")
    if abs(e3 - 2 * math.sqrt(2)) > 1e-12: fails.append("control-P3")
    # K_{1,2k^2} has energy 2*sqrt(2k^2) = 2k*sqrt(2): k = 1..5
    for k in range(1, 6):
        m = 2 * k * k; A = np.zeros((m + 1, m + 1)); A[0, 1:] = 1; A[1:, 0] = 1
        e, _ = energy(A)
        if abs(e - 2 * k * math.sqrt(2)) > 1e-9: fails.append(("control-star", k))
    print("control stars K_{1,2k^2}, k=1..5: energy 2k*sqrt(2) as expected" if not any(isinstance(f, tuple) and f[0] == "control-star" for f in fails) else "control stars FAIL")
    # C_8: energy 4 + 4 sqrt 2 (irrational square), consistent with the theorem
    C8 = np.zeros((8, 8))
    for i in range(8): C8[i, (i + 1) % 8] = C8[(i + 1) % 8, i] = 1
    print(f"C8 energy = {energy(C8)[0]:.10f} (4 + 4*sqrt2 = {4 + 4*math.sqrt(2):.10f})")
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:5]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
