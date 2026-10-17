"""Independent check of the length-506 counterexample to Glen–Simpson–Smyth Open Problem 1 (arXiv:1701.00928).
w = G(v), v = ababbababbabbababbabb, G(a) = UUV, G(b) = UUVV, U = aabaabab, V = aabab.
1. w is Lyndon (definition: strictly smaller than every proper suffix) and has exactly 16 distinct Lyndon factors
   (two independent counts: all distinct factors tested by the definition; Lyndon prefixes of every suffix).
2. All 220 lower Christoffel words of length 506 (gcd(p, 506) = 1): the direct count of distinct Lyndon factors equals
   1 + (sum of the partial quotients of p/506), and the minimum is 17.
3. Morphism identities L(A x) = L(x)+1, L(B x) = L(x)+1, L(h x) = L(x)+4 on random binary words containing both letters.
"""
import random, sys
from math import gcd

random.seed(20261017)


def is_lyndon(s):
    return len(s) > 0 and all(s < s[k:] for k in range(1, len(s)))


def lyndon_prefix_lengths(t):
    """Lengths l with t[:l] Lyndon (Duval-style scan)."""
    out = [1]; k = 0
    for j in range(1, len(t)):
        if t[k] < t[j]: k = 0; out.append(j + 1)
        elif t[k] == t[j]: k += 1
        else: break
    return out


def L_fast(s):
    return len({s[i:i + l] for i in range(len(s)) for l in lyndon_prefix_lengths(s[i:])})


def L_brute(s):
    facs = {s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)}
    return len([t for t in facs if is_lyndon(t)]), sorted((t for t in facs if is_lyndon(t)), key=len)


def christoffel(p, n):
    return "".join("b" if ((i + 1) * p) // n - (i * p) // n else "a" for i in range(n))


def qsum(a, b):
    s = 0
    while b: s += a // b; a, b = b, a % b
    return s


A = {"a": "a", "b": "ab"}; B = {"a": "ab", "b": "b"}; H = {"a": "aab", "b": "aabb"}
app = lambda m, x: "".join(m[c] for c in x)


def main():
    fails = []
    U, V = "aabaabab", "aabab"; v = "ababbababbabbababbabb"
    G = {"a": U + U + V, "b": U + U + V + V}
    w = app(G, v)
    # the morphism description agrees: v = BABAB(ab), G = (AB)^2 h
    v2 = app(B, app(A, app(B, app(A, app(B, "ab")))))
    w2 = app(A, app(B, app(A, app(B, app(H, v)))))
    print(f"|w| = {len(w)}, letters a/b = {w.count('a')}/{w.count('b')}, v = BABAB(ab): {v2 == v}, w = (AB)^2 h (v): {w2 == w}")
    if not (len(w) == 506 and v2 == v and w2 == w): fails.append("construction")
    lyn = is_lyndon(w)
    nb, facs = L_brute(w); nf = L_fast(w)
    print(f"w Lyndon: {lyn}; distinct Lyndon factors: brute {nb}, prefix-scan {nf}; lengths {[len(t) for t in facs]}")
    if not (lyn and nb == 16 and nf == 16): fails.append("w")
    cnt = {}
    for p in range(1, 506):
        if gcd(p, 506) != 1: continue
        c = christoffel(p, 506)
        if not is_lyndon(c): fails.append(("christoffel-not-lyndon", p))
        cnt[p] = L_fast(c)
        if cnt[p] != qsum(506, p) + 1 - 0 and cnt[p] != qsum(p, 506 - p) + 0:  # report mismatch to the quotient formula
            pass
    qs = {p: qsum(p, 506) for p in cnt}  # p/506 = [0; q1, ...]: qsum(p, 506) = 0 + q1 + ... (first quotient of p//506 = 0)
    mism = [p for p in cnt if cnt[p] != qs[p] + 1]
    print(f"Christoffel words of length 506: {len(cnt)}; min direct count = {min(cnt.values())}; "
          f"argmin count = {sum(1 for x in cnt.values() if x == min(cnt.values()))}; direct = quotient sum + 1 for all: {not mism}")
    if len(cnt) != 220 or min(cnt.values()) != 17 or mism: fails.append(("christoffel", min(cnt.values()), mism[:5]))
    # spot brute-force check of the direct count on a few Christoffel words
    for p in (97, 1, 253 if gcd(253, 506) == 1 else 251):
        b, _ = L_brute(christoffel(p, 506))
        if b != cnt[p]: fails.append(("brute-christoffel", p))
    print("brute-force spot check of Christoffel counts (p = 97, 1, 251): ok" if not any(isinstance(f, tuple) and f[0] == "brute-christoffel" for f in fails) else "spot check FAIL")
    # morphism identities
    bad = 0
    for _ in range(300):
        x = "".join(random.choice("ab") for _ in range(random.randint(2, 14)))
        if "a" not in x or "b" not in x: continue
        Lx = L_fast(x)
        if L_fast(app(A, x)) != Lx + 1 or L_fast(app(B, x)) != Lx + 1 or L_fast(app(H, x)) != Lx + 4: bad += 1
    print(f"morphism identities on random words: {bad} failures")
    if bad: fails.append("morphisms")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
