"""Checks for Conjecture 1 of H. Manabe, "An Amalgamation Nim with Restriction", arXiv:2411.14488v1.
Moves (Definition 5): lower one heap strictly; merge two heaps that are both >= 2, the sum going to either of the two
places and the other becoming 0. g = mex of the option values.
(a) Theorem 2: the zero set of g equals the set P of Definition 6 (as transcribed).
(b) Conjecture 1: floor(g/2) = floor((x xor y xor z)/2).
(c) The search seat's exact formula g = 2h + (p xor delta)."""
import sys
from functools import lru_cache
sys.setrecursionlimit(100000)

def moves(x, y, z):
    for u in range(x): yield (u, y, z)
    for v in range(y): yield (x, v, z)
    for w in range(z): yield (x, y, w)
    if x >= 2 and y >= 2: yield (x + y, 0, z); yield (0, x + y, z)
    if y >= 2 and z >= 2: yield (x, y + z, 0); yield (x, 0, y + z)
    if x >= 2 and z >= 2: yield (x + z, y, 0); yield (0, y, x + z)

@lru_cache(maxsize=None)
def g(x, y, z):
    s = {g(*m) for m in moves(x, y, z)}
    k = 0
    while k in s: k += 1
    return k

def inP(x, y, z):
    def P01(a, b, c): return a + b == c and a <= c and b <= c and a ^ b ^ c == 0 and min(a, b) < 2
    def P02(a, b, c): return a + b > c + 2 and a <= c and b <= c and a ^ b ^ c == 0
    def N01(a, b, c): return a + b == c and a <= c and b <= c and a ^ b ^ c == 0 and a >= 2 and b >= 2
    def N02(a, b, c): return a + b == c + 2 and a <= c and b <= c and a ^ b ^ c == 0
    def P11(a, b, c):
        return ((a + b) % 2 == 0 and c >= 1 and N01(a, b, c - 1)) or ((a + b) % 2 == 1 and N01(a, b, c + 1))
    def P12(a, b, c):
        return ((a + b) % 2 == 0 and c >= 1 and N02(a, b, c - 1)) or ((a + b) % 2 == 1 and N02(a, b, c + 1))
    def anyperm(f):
        return f(x, y, z) or f(x, z, y) or f(y, z, x)
    return anyperm(lambda a, b, c: P01(a, b, c) or P02(a, b, c)) or anyperm(lambda a, b, c: P11(a, b, c) or P12(a, b, c))

def formula(x, y, z):
    a = [x // 2, y // 2, z // 2]; e = [x % 2, y % 2, z % 2]
    h = a[0] ^ a[1] ^ a[2]; p = e[0] ^ e[1] ^ e[2]
    pairs = [(0, 1), (0, 2), (1, 2)]
    B = any(a[i] > 0 and a[j] > 0 and (a[i] & a[j]) == 0 for i, j in pairs)
    E = all(v > 0 for v in a) and all((a[i] & a[j]) == 0 for i, j in pairs)
    delta = int(B and not E)
    return 2 * h + (p ^ delta)

def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    bad_thm2 = bad_conj = bad_form = 0
    n = 0
    for x in range(K + 1):
        for y in range(K + 1):
            for z in range(K + 1):
                v = g(x, y, z); n += 1
                if (v == 0) != inP(x, y, z): bad_thm2 += 1
                if v // 2 != (x ^ y ^ z) // 2: bad_conj += 1
                if v != formula(x, y, z): bad_form += 1
    print(f"positions with x,y,z <= {K}: {n}")
    print(f"(a) Theorem 2 (zero set = P of Definition 6): mismatches {bad_thm2}")
    print(f"(b) Conjecture 1: mismatches {bad_conj}")
    print(f"(c) exact formula 2h + (p xor delta): mismatches {bad_form}")
    print(f"examples: g(0,2,4)={g(0,2,4)} (xor 6); g(2,2,0)={g(2,2,0)}; g(2,4,6)={g(2,4,6)} (xor 0)")
    bad = bad_conj + bad_form
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
