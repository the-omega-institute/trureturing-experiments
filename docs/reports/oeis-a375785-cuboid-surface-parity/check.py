"""OEIS A375785 (Felix Huber): a(n) = number of unordered positive solutions of xy + yz + zx = 3n^2. Conjecture: all terms are odd.
1. Brute force for n <= 300 (x <= y <= z, z = (N - xy)/(x + y)): all a(n) odd; the first 69 values agree with the OEIS data.
2. The reduction a(n) = H(3n^2) - d(3n^2)/2, H = number of reduced positive definite forms [a, b, c] of discriminant -12 n^2
   (|b| <= a <= c, b >= 0), via the bijection (x, y, z) -> [x + y, 2x, x + z] on nonnegative solutions, checked for n <= 120.
3. Orbit-type split: o1 = 1 (only (n, n, n)), and o3 + o6 even, listed for n <= 30.
"""
import json, os, sys
from math import isqrt

def a_brute(n):
    N = 3 * n * n; cnt = 0; o3 = 0; o6 = 0; o1 = 0
    x = 1
    while 3 * x * x <= N:
        y = x
        while x * y + y * y + x * y <= N:          # z >= y  =>  xy + 2yz... use xy + y*z + x*z >= xy + y^2 + x*y
            num = N - x * y; den = x + y
            if num > 0 and num % den == 0:
                z = num // den
                if z >= y:
                    cnt += 1
                    if x == y == z: o1 += 1
                    elif x == y or y == z: o3 += 1
                    else: o6 += 1
            y += 1
        x += 1
    return cnt, o1, o3, o6

def reduced_forms(D):  # D = -disc > 0; count reduced forms with b >= 0, |b| <= a <= c
    cnt = 0
    a = 1
    while 3 * a * a <= D:
        for b in range(0, a + 1):
            if (b * b + D) % (4 * a) == 0:
                c = (b * b + D) // (4 * a)
                if c >= a: cnt += 1
        a += 1
    return cnt

def ndiv(m):
    c = 0; i = 1
    while i * i <= m:
        if m % i == 0: c += 1 if i * i == m else 2
        i += 1
    return c

def main():
    fails = []
    here = os.path.dirname(os.path.abspath(__file__))
    vals = {}
    for n in range(1, 301):
        cnt, o1, o3, o6 = a_brute(n); vals[n] = (cnt, o1, o3, o6)
        if cnt % 2 == 0: fails.append(("even", n))
        if o1 != 1 or (o3 + o6) % 2: fails.append(("split", n))
    data = None
    if os.path.exists(os.path.join(here, 'a375785.json')):
        d = json.load(open(os.path.join(here, 'a375785.json'))); r = d[0] if isinstance(d, list) else d['results'][0]
        data = [int(v) for v in r['data'].split(',')]
        if [vals[n][0] for n in range(1, len(data) + 1)] != data: fails.append("oeis")
    print(f"a(n) odd for all n <= 300: {not any(f[0] == 'even' for f in fails if isinstance(f, tuple))}; matches OEIS data ({len(data) if data else 0} terms): {'oeis' not in fails}")
    bad = [n for n in range(1, 121) if vals[n][0] != reduced_forms(12 * n * n) - ndiv(3 * n * n) // 2]
    print(f"a(n) = H(3n^2) - d(3n^2)/2 for n <= 120: mismatches {bad}")
    if bad: fails.append("reduction")
    print("n: (a, o1, o3, o6) for n <= 30:", {n: vals[n] for n in range(1, 31)})
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:5]}")
    sys.exit(0 if not fails else 1)

if __name__ == "__main__":
    main()
