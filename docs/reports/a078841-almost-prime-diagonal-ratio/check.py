#!/usr/bin/env python3
"""Exact finite checks of Wilson's two diagonal almost-prime ratio candidates.

Inputs: the official, SHA-256-pinned OEIS b-files credited in README.md.
Arithmetic: Python integers and fractions.Fraction only. No limit is proved.
"""
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from urllib.request import Request, urlopen

TABLES = {
    "A078841": (
        "https://oeis.org/A078841/b078841.txt",
        "dbd55c9151f0b2840b36efe4ae1d93e7cccd6249f73ba636f7c6aa5d860b5c55",
    ),
    "A101695": (
        "https://oeis.org/A101695/b101695.txt",
        "62aa01f436c5e62b6625b729b95e0b5bc400801ced1583617a538f874dc2b762",
    ),
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def official_table(name):
    url, expected_hash = TABLES[name]
    with urlopen(Request(url, headers={"User-Agent": "Wilson-OEIS-exact-check/1"}), timeout=60) as response:
        raw = response.read()
    digest = sha256(raw).hexdigest()
    require(digest == expected_hash, f"{name}: source hash changed: {digest}")
    table = {}
    for line in raw.decode("utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            n, a = map(int, line.split())
            require(n not in table, f"{name}: duplicate index")
            table[n] = a
    print(f"SOURCE {name} terms={len(table)} sha256={digest}")
    return table


def primes_through(limit):
    prime = bytearray(b"\x01") * (limit + 1)
    prime[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            for m in range(p * p, limit + 1, p):
                prime[m] = 0
    return [p for p in range(2, limit + 1) if prime[p]]


def monoid_through(bound):
    """Enumerate every nondecreasing odd-prime word whose product of p/2 <= bound.

    All factors >= 3/2, so pruning is exact. Every used p satisfies p <= 2*bound.
    Reduced denominators retain the word length, independently checked below.
    """
    primes = [p for p in primes_through((2 * bound).numerator // (2 * bound).denominator) if p != 2]
    lengths = {}

    def visit(value, first, length):
        require(value not in lengths, f"nonunique odd-prime word: {value}")
        require(value.denominator == 1 << length, f"length lost in reduction: {value}")
        lengths[value] = length
        for i in range(first, len(primes)):
            child = value * Fraction(primes[i], 2)
            if child > bound:
                break
            visit(child, i, length + 1)

    visit(Fraction(1), 0, 0)
    return sorted(lengths), lengths, primes[-1]


def omega_sieve(limit):
    """Independent integer sieve: add one for each prime power dividing m."""
    omega = bytearray(limit + 1)
    for p in range(2, limit + 1):
        if omega[p] == 0:
            power = p
            while power <= limit:
                for m in range(power, limit + 1, power):
                    omega[m] += 1
                power *= p
    return omega


def main():
    a = official_table("A078841")
    d = official_table("A101695")
    require(set(a) == set(range(1001)), "A078841 index range")
    require(set(d) == set(range(1, 1001)), "A101695 index range")
    require([a[n] for n in range(6)] == [1, 3, 9, 20, 54, 112], "A078841 published examples")
    require([d[n] for n in range(1, 6)] == [2, 6, 18, 40, 108], "A101695 published examples")

    bound = Fraction(4437, 16)
    b, lengths, largest_prime = monoid_through(bound)
    require(len(b) >= 1001, "monoid bound does not reach rank 1000")
    max_length = 0
    for n in range(1001):
        max_length = max(max_length, lengths[b[n]])
        require(max_length <= n, f"length filtration at rank {n}")
        require(b[n] <= Fraction(3, 2) ** n, f"power upper bound at rank {n}")
        require(Fraction(a[n], 1 << n) == b[n], f"main rank identity at n={n}")
    print(f"MONOID bound={bound} elements={len(b)} largest_odd_prime={largest_prime} b_1000={b[1000]}")
    print("FIDELITY all 1001 official A078841 terms match independent Fraction enumeration")
    for n in range(1000):
        require(Fraction(d[n + 1], 1 << (n + 1)) == b[n], f"companion rank at n={n}")
        require(d[n + 1] == 2 * a[n], f"companion identity at n={n}")
    print("FIDELITY all 1000 official A101695 terms match Fraction enumeration")
    print("IDENTITY A101695(n+1)=2*A078841(n) checked for every n=0..999")

    # Recompute source ranks without the monoid identification.
    last = 16
    limit = max(a[last], d[last + 1])
    omega = omega_sieve(limit)
    rows = [[] for _ in range(last + 2)]
    for m in range(1, limit + 1):
        k = omega[m]
        if k < len(rows) and len(rows[k]) < last + 2:
            rows[k].append(m)
    for n in range(last + 1):
        require(len(rows[n]) > n, f"missing main brute rank at {n}")
        require(rows[n][n] == a[n] == b[n] * (1 << n), f"independent main rank at {n}")
        require(len(rows[n + 1]) > n, f"missing companion brute rank at {n}")
        require(rows[n + 1][n] == d[n + 1], f"independent companion rank at {n}")
    require(rows[0] == [1], "Omega=0 positive-number edge case")
    print(f"INDEPENDENT_INTEGER_SIEVE limit={limit} main_n=0..{last} companion_n=1..{last+1}")
    print("EDGE_CASE Omega(1)=0; zero excluded; A078841(0)=1; A101695(1)=2")
    print(f"DIAGNOSTIC exact_main_ratio_999={Fraction(a[1000], a[999])}")
    print("SCOPE finite fidelity and rank checks; neither asymptotic limit is certified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
