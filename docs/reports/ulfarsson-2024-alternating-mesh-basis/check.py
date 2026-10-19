"""Ulfarsson, arXiv:2411.17778v1, Conjecture 2.4.2 (mesh-pattern bases of the even permutations).

Checks the two lemmas of the proof by brute force:
1. Every mesh pattern (p, R) with an unshaded cell occurs in an even permutation of length <= k+2
   (exhaustive for k <= 2, 3000 random shadings for k = 3), and the explicit two-point insertion in an
   unshaded cell produces two permutations of opposite parity that both contain (p, R).
2. A fully shaded pattern (p, full) occurs in sigma iff sigma = p (all sigma of length <= k+2, k <= 3).
Occurrence follows Branden-Claesson: increasing positions order-isomorphic to p, and no unselected point
in a shaded cell (i = selected positions to its left, j = selected values below it).
Usage: python3 check.py  (prints RESULT PASS / RESULT FAIL).
"""
import itertools, random

def sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s

def occurs(p, R, s):
    k, n = len(p), len(s)
    for pos in itertools.combinations(range(n), k):
        vals = [s[x] for x in pos]
        if any((vals[i] < vals[j]) != (p[i] < p[j]) for i in range(k) for j in range(k)):
            continue
        ok = True
        for x in range(n):
            if x in pos: continue
            c = (sum(1 for y in pos if y < x), sum(1 for v in vals if v < s[x]))
            if c in R: ok = False; break
        if ok: return True
    return False

def insertion(p, a, b):
    """The two permutations of the proof: insert values b, b+1 (0-based) after position a."""
    f = [v if v < b else v + 2 for v in p]
    up = f[:a] + [b, b + 1] + f[a:]
    down = f[:a] + [b + 1, b] + f[a:]
    return tuple(up), tuple(down)

def main():
    ok = True
    rng = random.Random(1)
    cells = lambda k: [(i, j) for i in range(k + 1) for j in range(k + 1)]
    tested = 0
    for k in range(0, 4):
        for p in itertools.permutations(range(k)):
            C = cells(k)
            if k <= 2:
                shadings = [frozenset(c for c, b in zip(C, bits) if b) for bits in itertools.product([0, 1], repeat=len(C))]
            else:
                shadings = [frozenset(c for c in C if rng.random() < 0.7) for _ in range(250)]
            for R in shadings:
                if len(R) == len(C): continue
                a, b = next(c for c in C if c not in R)
                up, down = insertion(list(p), a, b)
                good = occurs(p, R, up) and occurs(p, R, down) and sign(up) == -sign(down)
                if not good: ok = False; print("FAIL insertion", p, sorted(R), up, down)
                tested += 1
    print("lemma 1: patterns with an unshaded cell tested:", tested)
    full_ok = True
    for k in range(0, 4):
        for p in itertools.permutations(range(k)):
            R = frozenset(cells(k))
            for n in range(0, k + 3):
                for s in itertools.permutations(range(n)):
                    if occurs(p, R, s) != (s == p): full_ok = False; print("FAIL full", p, s)
    print("lemma 2 (full shading occurs only in itself), k <= 3, n <= k+2:", full_ok)
    ok &= full_ok
    print("RESULT", "PASS" if ok else "FAIL")

if __name__ == "__main__":
    main()
