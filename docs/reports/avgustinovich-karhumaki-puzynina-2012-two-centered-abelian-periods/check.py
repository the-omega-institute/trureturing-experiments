"""Finite check of Open problem 6.2 of Avgustinovich-Karhumaki-Puzynina, RAIRO-ITA 46 (2012) 3-15.

Question: does some non-periodic bi-infinite word contain, at every position (cut), a central
abelian 2k-power with period length l1 or l2?  Proposed answer (k = 1 suffices): no; every
bi-infinite word with a centered abelian square of block length p or q at every cut has period
lcm(p, q).
Method: the shift graph on windows of length 2M (M = max(p, q)) whose centre cut satisfies the
hypothesis.  Bi-infinite words satisfying the hypothesis everywhere are exactly bi-infinite
paths.  Prediction registered before the run: for every pair checked, (i) every nontrivial
strongly connected component is a single cycle (no branching), (ii) the periodic word of each
cycle has period lcm(p, q), (iii) no directed path joins two distinct nontrivial components.
Negative control: with every length 1..8 allowed (source Example 3.9: non-periodic (1,8)-ACR
words exist) the same test must fail.
Usage: python3 check.py   (standard library only)
"""
import sys
from math import lcm
from itertools import product
from collections import Counter


def centered_square(win, M, r):
    # win has length 2M; the cut is between win[M-1] and win[M]
    left = win[M - r:M]
    right = win[M:M + r]
    return Counter(left) == Counter(right)


def analyse(alpha, p, q, lengths=None):
    lengths = lengths or (p, q)
    M = max(lengths)
    L = 2 * M
    valid = [w for w in product(range(alpha), repeat=L) if any(centered_square(w, M, r) for r in lengths)]
    vs = set(valid)
    succ = {w: [w[1:] + (a,) for a in range(alpha) if w[1:] + (a,) in vs] for w in valid}
    # Tarjan SCC (iterative)
    index, low, onstack, stack, comps, idx = {}, {}, set(), [], [], [0]
    for root in valid:
        if root in index:
            continue
        work = [(root, iter(succ[root]))]
        index[root] = low[root] = idx[0]; idx[0] += 1; stack.append(root); onstack.add(root)
        while work:
            v, it = work[-1]
            nxt = next(it, None)
            if nxt is None:
                work.pop()
                if work:
                    u = work[-1][0]; low[u] = min(low[u], low[v])
                if low[v] == index[v]:
                    comp = []
                    while True:
                        x = stack.pop(); onstack.discard(x); comp.append(x)
                        if x == v:
                            break
                    comps.append(comp)
            elif nxt not in index:
                index[nxt] = low[nxt] = idx[0]; idx[0] += 1; stack.append(nxt); onstack.add(nxt)
                work.append((nxt, iter(succ[nxt])))
            elif nxt in onstack:
                low[v] = min(low[v], index[nxt])
    comp_of = {}
    nontrivial = []
    for c in comps:
        cs = set(c)
        if len(c) > 1 or any(w in cs for w in succ[c[0]]):
            nontrivial.append(cs)
        for w in c:
            comp_of[w] = len(nontrivial) - 1 if (len(c) > 1 or any(x in cs for x in succ[c[0]])) else None
    m = lcm(*lengths)
    ok = True
    for cs in nontrivial:
        # (i) no branching inside the component
        if any(sum(1 for x in succ[w] if x in cs) != 1 for w in cs):
            ok = False
        # (ii) the cycle word has period m: window w and its m-th successor agree in their first letter shifted
        for w in cs:
            v = w
            for _ in range(m):
                v = next(x for x in succ[v] if x in cs)
            if v != w:
                ok = False
                break
    # (iii) no path from one nontrivial component to another
    comp_index = {}
    for k, cs in enumerate(nontrivial):
        for w in cs:
            comp_index[w] = k
    for k, cs in enumerate(nontrivial):
        seen = set(cs); frontier = list(cs)
        while frontier:
            v = frontier.pop()
            for x in succ[v]:
                if x in seen:
                    continue
                if x in comp_index and comp_index[x] != k:
                    ok = False
                seen.add(x); frontier.append(x)
    return ok, len(valid), len(nontrivial)


def main():
    ok = True
    for alpha, top in ((2, 6), (3, 3)):
        for q in range(2, top + 1):
            for p in range(1, q):
                good, nv, nc = analyse(alpha, p, q)
                ok &= good
                print(f"alphabet {alpha}, (p, q) = ({p}, {q}): {nv} valid windows, {nc} cyclic components, "
                      f"{'ok' if good else 'FAIL'}")
    # negative control (source Example 3.9): all lengths 1..8 admit non-periodic (1,8)-ACR words,
    # so the same test must fail there
    good, nv, nc = analyse(2, 0, 0, tuple(range(1, 9)))
    print(f"control: alphabet 2, lengths 1..8: {nv} valid windows, {nc} cyclic components, test {'passes (unexpected)' if good else 'fails, as required'}")
    ok &= not good
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
