"""Exact good ideals and sextic symbols in Z[omega], omega**2+omega+1=0.

Ideals use their unique generator congruent to 1 modulo 3; good means
coprime to 2 and 3. Symbols return exponents of zeta_6=1+omega, with -1
for a natural character zero. Only integer arithmetic is used here.
"""
from __future__ import annotations
from math import isqrt

Pair = tuple[int, int]


def norm(x: Pair) -> int:
    a, b = x
    return a*a - a*b + b*b


def mul(x: Pair, y: Pair) -> Pair:
    a, b = x
    c, d = y
    return a*c - b*d, a*d + b*c - b*d


UNITS = ((1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))


def primary(x: Pair) -> Pair:
    """Canonical associate, defined exactly when the norm is prime to 3."""
    if norm(x) % 3 == 0:
        raise ValueError("primary requires norm prime to 3")
    return next(y for u in UNITS if (y := mul(u, x))[0] % 3 == 1
                and y[1] % 3 == 0)


def div_exact(x: Pair, y: Pair) -> Pair | None:
    n = norm(y)
    if not n:
        raise ZeroDivisionError("Eisenstein division by zero")
    a, b = mul(x, (y[0] - y[1], -y[1]))
    return (a // n, b // n) if a % n == b % n == 0 else None


def _rational_factors(n: int) -> dict[int, int]:
    result = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        result[n] = 1
    return result


def _ideals(limit: int):
    # (2a-b)^2+3b^2=4N makes the ellipse enumeration exact.
    for b in range(-isqrt(4*limit // 3), isqrt(4*limit // 3) + 1):
        if b % 3:
            continue
        r = isqrt(4*limit - 3*b*b)
        for a in range(-((r-b) // 2), (b+r) // 2 + 1):
            if a % 3 == 1 and norm((a, b)) % 6 == 1:
                yield a, b


def _split_primes(p: int) -> tuple[Pair, ...]:
    result = set()
    for b in range(-isqrt(4*p // 3), isqrt(4*p // 3) + 1):
        if b % 3:
            continue
        r = isqrt(4*p - 3*b*b)
        if r*r != 4*p - 3*b*b:
            continue
        for c in (b-r, b+r):
            if c % 2 == 0 and (c // 2) % 3 == 1:
                result.add((c // 2, b))
    if len(result) != 2:
        raise ArithmeticError("split prime representation failed")
    return tuple(sorted(result))


def _modmul(x: Pair, y: Pair, p: int) -> Pair:
    a, b = mul(x, y)
    return a % p, b % p


def _modpow(x: Pair, k: int, p: int) -> Pair:
    result = (1, 0)
    while k:
        if k & 1:
            result = _modmul(result, x, p)
        x = _modmul(x, x, p)
        k //= 2
    return result


class Engine:
    """Enumerate all good ideals up to limit; factor extra inputs on demand.

    Additional factors are cached in .factors but never added to .ideals.
    Non-primary associates are accepted by the public ideal methods.
    """
    def __init__(self, limit: int):
        if not isinstance(limit, int) or limit < 1:
            raise ValueError("limit must be a positive integer")
        self.limit = limit
        self.ideals = sorted(_ideals(limit), key=lambda x: (norm(x), x))
        self.factors: dict[Pair, dict[Pair, int]] = {}
        self._split = {}
        self._fields = {}
        self._symbols = {}
        for x in self.ideals:
            self._factor(x)

    def _factor(self, x: Pair) -> dict[Pair, int]:
        x = primary(x)
        if norm(x) % 6 != 1:
            raise ValueError("ideal must be prime to 2 and 3")
        if x not in self.factors:
            left, factors = x, {}
            for q in _rational_factors(norm(x)):
                if q % 3 == 2:
                    candidates = (primary((q, 0)),)
                else:
                    if q not in self._split:
                        self._split[q] = _split_primes(q)
                    candidates = self._split[q]
                for p in candidates:
                    while (quotient := div_exact(left, p)) is not None:
                        factors[p] = factors.get(p, 0) + 1
                        left = quotient
            if left != (1, 0):
                raise ArithmeticError("incomplete ideal factorization")
            self.factors[x] = factors
        return self.factors[x]

    def mu(self, x: Pair) -> int:
        values = self._factor(x).values()
        return 0 if any(e > 1 for e in values) else (-1)**len(values)

    def liouville(self, x: Pair) -> int:
        return (-1)**sum(self._factor(x).values())

    def is_sixth_free(self, x: Pair) -> bool:
        return all(e < 6 for e in self._factor(x).values())

    def _field(self, prime: Pair):
        prime = primary(prime)
        if prime not in self._fields:
            if self._factor(prime) != {prime: 1}:
                raise ValueError("expected a good prime ideal")
            n = norm(prime)
            p = isqrt(n) if isqrt(n)**2 == n else n
            root = None if p*p == n else -prime[0]*pow(prime[1], -1, p) % p
            roots = {self._residue_raw(u, p, root): k for k, u in enumerate(UNITS)}
            if len(roots) != 6:
                raise ArithmeticError("sixth roots are not distinct")
            self._fields[prime] = p, root, roots
        return self._fields[prime]

    @staticmethod
    def _residue_raw(x: Pair, p: int, root: int | None) -> Pair:
        return ((x[0] + x[1]*root) % p, 0) if root is not None else (x[0] % p, x[1] % p)

    def prime_symbol_exponent(self, prime: Pair, row: Pair) -> int:
        prime = primary(prime)
        p, root, roots = self._field(prime)
        residue = self._residue_raw(row, p, root)
        if residue == (0, 0):
            return -1
        key = prime, residue
        if key not in self._symbols:
            power = _modpow(residue, (norm(prime)-1) // 6, p)
            self._symbols[key] = roots[power]
        return self._symbols[key]

    def symbol_exponent(self, denom: Pair, row: Pair) -> int:
        result = 0
        for p, e in self._factor(denom).items():
            k = self.prime_symbol_exponent(p, row)
            if k == -1:
                return -1
            result += e*k
        return result % 6

    def unit_weight(self, denom: Pair) -> int:
        """Exponent of chi_d(zeta_6), preserving each prime valuation."""
        return sum(e*((norm(p)-1) // 6) for p, e in self._factor(denom).items()) % 6
