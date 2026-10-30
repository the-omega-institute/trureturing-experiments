#!/usr/bin/env python3
"""Exact finite checks for OEIS A048585, ceiling Pisot L(6,7).

Source: https://oeis.org/A048585/internal, revision 18.
OEIS data are attributed to the OEIS Foundation and David W. Wilson,
under CC BY-SA 4.0 (https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The finite computations here are not an all-index proof or a Lean certificate.
"""
from fractions import Fraction
from hashlib import sha256
import argparse

SOURCE_PREFIX = (6,7,9,12,16,22,31,44,63,91,132,192,280,409,598,875,1281,1876,2748,4026,5899,8644,12667,18563,27204,39868,58428,85629,125494,183919,269545,395036,578952,848494,1243527,1822476,2670967,3914491,5736964,8407928,12322416,18059377,)


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def ceiling(q):
    require(isinstance(q, Fraction), "rational ceiling input")
    return -((-q.numerator) // q.denominator)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terms", type=int, default=10002)
    args = parser.parse_args()
    require(args.terms >= len(SOURCE_PREFIX) and args.terms >= 15, "sufficient terms")
    a = [6, 7]
    quotients = []
    for n in range(args.terms - 2):
        require(a[n] > 0, f"positive denominator at {n}")
        q = Fraction(a[n + 1] ** 2, a[n])
        quotients.append(q)
        term = ceiling(q)
        require(term == (a[n + 1] ** 2 + a[n] - 1) // a[n], f"ceiling cross-check {n}")
        a.append(term)
    require(tuple(a[:len(SOURCE_PREFIX)]) == SOURCE_PREFIX, "verbatim OEIS displayed data")
    require(quotients[2] == 16 and a[4] == 16, "integer quotient unchanged")
    for n in range(4, len(a)):
        require(a[n] == 2*a[n-1] - a[n-2] + a[n-3] - a[n-4], f"published identity {n}")

    v = [3, 4, 6]
    for n in range(args.terms - 3):
        v.append(v[n + 2] + v[n])
    u = [x + 3 for x in v]
    require(a == u, "independent linear recursion agrees with direct ceiling")
    s = [v[n+2] - 2*v[n+1] + v[n] for n in range(len(v)-2)]
    d = [v[n]*v[n+2] - v[n+1]**2 for n in range(len(v)-2)]
    e = [u[n]*u[n+2] - u[n+1]**2 for n in range(len(u)-2)]
    for n in range(len(s)):
        require(s[n] > 0, f"second difference positive {n}")
        require(e[n] == d[n]+3*s[n], f"translation identity {n}")
        require(0 <= e[n] < u[n], f"ceiling error window {n}")
        require(Fraction(u[n+2]-1) < quotients[n] <= u[n+2], f"ceiling characterization {n}")
        if n+1 < len(s):
            require(v[n] < v[n+1] and s[n] <= s[n+1], f"monotonicity {n}")
        if n+3 < len(s):
            require(s[n+3] == s[n+2]+s[n], f"second difference recursion {n}")
            require(d[n+3] == d[n]-d[n+1], f"determinant recursion {n}")
            require(v[n+3]-4*s[n+3] == (v[n+2]-4*s[n+2])+(v[n]-4*s[n]), f"margin recursion {n}")
        if n >= 7:
            require(abs(d[n]) <= s[n] and v[n] > 4*s[n], f"induction bounds {n}")
            require(0 < 2*s[n] <= e[n] <= 4*s[n] < v[n] < u[n], f"final bounds {n}")
    bases = [(v[n], s[n], d[n], v[n]-4*s[n]) for n in (7,8,9)]
    require(bases == [(41,9,8,5),(60,13,-4,8),(88,19,-9,12)], "three induction bases")
    initial = list(zip(e[:7],u[:7]))
    require(initial == [(5,6),(3,7),(0,9),(8,12),(12,16),(7,22),(17,31)], "seven initial error windows")
    for n in range(len(u)-3):
        require(u[n+3] == u[n+2]+u[n]-3, f"affine recursion {n}")
    digest = sha256("\n".join(map(str,a[:1002])).encode()).hexdigest()
    print(f"FIDELITY: all {len(SOURCE_PREFIX)} displayed OEIS A048585 terms agree")
    print("ROUNDING: 12^2/9 = 16; ceiling leaves 16 unchanged")
    print(f"DIRECT CEILING VS LINEAR: {len(a)} terms, indices 0..{len(a)-1}, agree")
    print(f"PUBLISHED ORDER-FOUR IDENTITY: indices 4..{len(a)-1}, all agree")
    print(f"DETERMINANTS AND ERROR BOUNDS: indices 0..{len(e)-1}, all agree")
    print(f"BASES (v,S,D,v-4S) at 7,8,9: {bases}")
    print(f"INITIAL (E,U) at 0..6: {initial}")
    print(f"FIRST_1002_TERMS_SHA256 (decimal lines, no final newline): {digest}")
    print("SCOPE: exact finite checks only; the unbounded outline is not Lean-verified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
