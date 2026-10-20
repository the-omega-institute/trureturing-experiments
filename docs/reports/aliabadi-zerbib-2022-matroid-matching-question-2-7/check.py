"""Aliabadi-Zerbib, arXiv:2202.07719v3, Question 2.7: a negative answer over G = Z.

For n, s >= 2, E = {1, ..., n + s}, M = U_{n, n+s} and N the rank-n matroid on E whose bases are all n-subsets except
H = {s + 1, ..., s + n}. Checks by brute force for (n, s) in {2, ..., 5}^2: N satisfies the basis axioms, is loopless,
paving, and its dual is paving (so N is sparse paving); 0 is not in E(N); |E(M)| = |E(N)|; the basis (1, ..., n) of M
has no ordered basis (b_1, ..., b_n) of N with a_i + b_i not in E(M) for all i; and the same basis (1, ..., n) of N has
no ordered basis of N with a_i + b_i not in E(N), so N is not matched to itself (against Theorem 2.1 of the source).
"""
import itertools
import sys


def bases_N(n, s):
    E = range(1, n + s + 1)
    H = frozenset(range(s + 1, s + n + 1))
    return {frozenset(B) for B in itertools.combinations(E, n)} - {H}


def basis_exchange(Bs):
    for B in Bs:
        for C in Bs:
            for e in B - C:
                if not any((B - {e}) | {f} in Bs for f in C - B):
                    return False
    return True


def independent(Bs, X):
    return any(X <= B for B in Bs)


def paving(Bs, E, r):
    return all(independent(Bs, frozenset(X)) for X in itertools.combinations(E, r - 1))


ok = True
for n in (2, 3, 4, 5):
    for s in (2, 3, 4, 5):
        E = list(range(1, n + s + 1))
        Bs = bases_N(n, s)
        dual = {frozenset(E) - B for B in Bs}
        loopless = all(any(e in B for B in Bs) for e in E)
        sparse = paving(Bs, E, n) and paving(dual, E, s)
        exch = basis_exchange(Bs)
        a = list(range(1, n + 1))
        matched = any(frozenset(b) in Bs and all(a[i] + b[i] not in E for i in range(n))
                      for b in itertools.permutations(E, n))
        self_matched = any(all(a[i] + b[i] not in E for i in range(n)) and frozenset(b) in Bs
                           for b in itertools.permutations(E, n)) if frozenset(a) in Bs else None
        good = exch and loopless and sparse and 0 not in E and not matched and self_matched is False
        ok = ok and good
        print(f"n={n} s={s}: basis axioms {exch}, loopless {loopless}, sparse paving {sparse}, "
              f"basis {tuple(a)} of M matched to some basis of N: {matched}, of N to some basis of N: {self_matched}")
A = B = {1, 2, 3}
X = {1, 2, 3, 4, 5, 6}
sums = {a + b for a in A for b in B}
lemma_hyp = A <= X and B <= X and sums <= X and 0 not in X
lemma_fails = lemma_hyp and len(X) < len(A) + len(B) + 1
ok = ok and lemma_fails
print(f"Lemma 3.14 with A = B = {{1, 2, 3}}, X = {{1, ..., 6}} in Z: hypotheses {lemma_hyp}, |X| = {len(X)} < "
      f"|A| + |B| + 1 = {len(A) + len(B) + 1}: {lemma_fails}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
