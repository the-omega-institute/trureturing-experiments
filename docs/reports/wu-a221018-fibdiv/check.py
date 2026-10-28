#!/usr/bin/env python3
"""Exact, finite checks for the first FORMULA conjecture in OEIS A221018.

Python integers only. This is numerical evidence, not an unbounded proof.
The adjacent-index hypotheses are tested directly, without residue filters.
"""
from collections import Counter
from hashlib import sha256
from math import lcm
from pathlib import Path


LIMIT = 500_000
DIRECT_LIMIT = 20_000
SOURCE_URL = "https://oeis.org/A221018/b221018.txt"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def fib_pair_mod(n, modulus):
    """Return (F_n,F_{n+1}) mod modulus by binary fast doubling."""
    require(n >= 0 and modulus >= 1, "invalid Fibonacci arguments")
    a, b = 0, 1 % modulus
    for bit in bin(n)[2:]:
        c = a * (2 * b - a) % modulus
        d = (a * a + b * b) % modulus
        a, b = (c, d) if bit == "0" else (d, (c + d) % modulus)
    return a, b


def matrix_fib_mod(n, modulus):
    """Independent binary 2x2 matrix powering; no doubling formulas."""
    def mul(x, y):
        a, b, c, d = x
        e, f, g, h = y
        return ((a * e + b * g) % modulus,
                (a * f + b * h) % modulus,
                (c * e + d * g) % modulus,
                (c * f + d * h) % modulus)
    acc, base = (1, 0, 0, 1), (1, 1, 1, 0)
    while n:
        if n % 2:
            acc = mul(acc, base)
        base = mul(base, base)
        n //= 2
    return acc[1]


def belongs(k):
    require(k >= 2, "membership domain starts at 2")
    # lcm, rather than the product, is equivalent to the THREE conditions
    # for odd as well as even k. No parity assumption enters enumeration.
    return fib_pair_mod(k, lcm(k - 1, k, k + 1))[0] == 0


def main():
    source = Path(__file__).with_name("b221018.txt").read_bytes()
    rows = [tuple(map(int, line.split())) for line in source.decode().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    require(len(rows) == 10_000, "source must contain exactly 10000 terms")
    require([n for n, _ in rows] == list(range(1, 10_001)), "b-file indices")
    terms = [k for _, k in rows]
    require(all(a < b for a, b in zip(terms, terms[1:])), "strict ordering")
    print(f"source: {SOURCE_URL}")
    print(f"b-file sha256: {sha256(source).hexdigest()}")

    # Source fidelity: independently recompute the sequence, including gaps.
    found = [k for k in range(2, LIMIT + 1) if belongs(k)]
    require(found == [k for k in terms if k <= LIMIT], "b-file prefix/gaps mismatch")
    published_prefix = [108, 2688, 3528, 4968, 6480, 15288, 18048, 20808,
                        21000, 25308, 35448, 47520, 62928, 68208, 81648,
                        82008, 97608]
    require([k for k in found if k < 100_000] == published_prefix,
            "OEIS Python enumeration / displayed terms mismatch")
    print(f"fidelity: all integers 2 <= k <= {LIMIT}; {len(found)} members; exact b-file prefix")
    print(f"source Python range(2, 100000): {published_prefix}")

    # A direct full-integer recurrence checks every small hypothesis and the
    # fast-doubling implementation on every small case, independently.
    a, b = 0, 1
    direct_found = []
    for n in range(DIRECT_LIMIT + 1):
        if n >= 2:
            m = lcm(n - 1, n, n + 1)
            require(fib_pair_mod(n, m) == (a % m, b % m), "doubling/recurrence mismatch")
            member = all(a % d == 0 for d in (n - 1, n, n + 1))
            require(member == (a % m == 0), "lcm membership mismatch")
            if member:
                direct_found.append(n)
                require(n % 10 in (0, 2, 8), "claim fails on direct recurrence")
        a, b = b, a + b
    require(direct_found == [k for k in found if k <= DIRECT_LIMIT], "direct membership mismatch")
    print(f"independent full-integer recurrence: n=0..{DIRECT_LIMIT}; all pairs and memberships agree")

    # Source comment fidelity, plus checking all supplied listed examples.
    for k in terms:
        require(belongs(k), f"listed nonmember {k}")
        require(k % 10 in (0, 2, 8), f"claim fails at listed k={k}")
    histogram = dict(sorted(Counter(k % 10 for k in terms).items()))
    non8 = sum(k % 10 != 8 for k in terms[:7777])
    require(non8 == 81, "source's first-7777 count differs from b-file")
    print(f"fidelity: first 7777 terms have exactly {non8} terms not ending in 8")
    print(f"all 10000 listed members verified; final-digit histogram: {histogram}")
    for k in found:
        require(k % 10 in (0, 2, 8), f"small-case counterexample {k}")
        require(k % 12 == 0 and k % 5 not in (1, 4), "outline consequence fails")
    print(f"expected settlement: no counterexample for ANY integer 2 <= k <= {LIMIT}")

    # Finite-state ingredients of the proposed proof, including all zero
    # positions in a period. The unbounded induction belongs to the outline.
    for modulus, period, rank in ((2, 3, 3), (3, 8, 4), (5, 20, 5)):
        a, b = 0, 1
        values = []
        for _ in range(period):
            values.append(a)
            a, b = b, (a + b) % modulus
        require((a, b) == (0, 1), "state does not return")
        zeros = [i for i, value in enumerate(values) if value == 0]
        require(zeros == list(range(0, period, rank)), "wrong zero positions")
        for n in range(601):
            require((fib_pair_mod(n, modulus)[0] == 0) == (n % rank == 0),
                    "rank ingredient failed")
        print(f"rank ingredient: modulus={modulus}, period={period}, residues={values}, zero positions={zeros}")

    # A225219's stronger comment must not be treated as a prior theorem.
    k = 14_930_352
    m = k * (k - 1) * (k + 1)
    require(terms[125] == k, "A221018 term 126 mismatch")
    pair = fib_pair_mod(k, m)
    require(pair == (0, 1), "large witness pair")
    require(matrix_fib_mod(k, m) == pair[0], "large matrix check")
    a, b = 0, 1
    for _ in range(k):
        a, b = b, (a + b) % m
    require((a, b) == pair, "large direct recurrence disagrees")
    require(k % 30 == 12 and k % 10 == 2, "large witness residues")
    print(f"A225219 stronger-comment counterexample: k={k}, k%30={k % 30}, k%10={k % 10}")
    print(f"modulus k(k-1)(k+1)={m}; pair={pair}; doubling, matrix, and {k} direct steps agree")
    for k in terms[:20] + [terms[125], terms[-1]]:
        m = lcm(k - 1, k, k + 1)
        require(matrix_fib_mod(k, m) == 0, "sampled matrix check")
    print("scope: finite evidence and source fidelity only; the universal claim is not Lean-proved")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
