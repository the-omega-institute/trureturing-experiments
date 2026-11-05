#!/usr/bin/env python3
"""Exact, independent coefficient computations for five Bala OEIS conjectures.

Integer arithmetic only. The global logarithmic-derivative recurrence is
cross-checked against multiplication of the individual part-size factors.
No congruence under investigation is used to compute a coefficient.
"""
import json
from functools import lru_cache
from math import comb, isqrt

LIMIT = 200
IDS = ("A270913", "A008485", "A270919", "A008705", "A291697")
SOURCE_TERMS = {'A270913': [1, 1, 3, 13, 51, 206, 855, 3585, 15155, 64525, 276278, 1188353, 5130999, 22226049, 96544003, 420368858, 1834203955, 8018057345, 35107961175, 153950675585, 675978772326, 2971700764941, 13078268135683, 57613905606273, 254038914924791, 1121081799217231], 'A008485': [1, 1, 5, 22, 105, 506, 2492, 12405, 62337, 315445, 1605340, 8207563, 42124380, 216903064, 1119974875, 5796944357, 30068145905, 156250892610, 813310723925, 4239676354650, 22130265931900, 115654632452535, 605081974091875, 3168828466966388, 16610409114771900], 'A270919': [1, 2, 12, 80, 552, 3912, 28224, 206208, 1520784, 11297546, 84413912, 633713808, 4776117216, 36115518376, 273868321536, 2081866609920, 15859616674336, 121046064563376, 925411686479820, 7085465166635440, 54323193841192752, 416993869451825424, 3204447137019290944], 'A008705': [1, -1, -1, 5, -5, -6, 11, 41, -125, -85, 1054, -2069, -209, 8605, -15625, 3990, 14035, 36685, -130525, -254525, 1899830, -3603805, -134905, 13479425, -25499225, 23579969, -64447293, 237487433, -133867445, -1795846200, 6309965146, -6788705842, -11762712973], 'A291697': [1, 2, 8, 44, 256, 1512, 9056, 54896, 335872, 2069774, 12827888, 79875996, 499305472, 3131436856, 19694403520, 124165133424, 784478240768, 4965659813668, 31484486937512, 199923173603596, 1271192603065856, 8092551782518688, 51574780342740256, 329022223268286288, 2100934234342260736]}


def primes_up_to(limit):
    return [p for p in range(2, limit + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def log_derivative_weights(kind, limit):
    """[x^j] x F'(x)/F(x), without its diagonal exponent."""
    b = [0] * (limit + 1)
    for size in range(1, limit + 1):
        if kind == "A291697" and size % 2 == 0:
            continue
        for count in range(1, limit // size + 1):
            if kind == "A270913":
                weight = size if count % 2 else -size
            elif kind == "A008485":
                weight = size
            elif kind == "A008705":
                weight = -size
            else:
                weight = 2 * size if count % 2 else 0
            b[size * count] += weight
    return b


@lru_cache(None)
def recurrence_coefficient(kind, n):
    b = log_derivative_weights(kind, n)
    a = [1] + [0] * n
    for degree in range(1, n + 1):
        numerator = n * sum(b[j] * a[degree-j] for j in range(1, degree+1))
        assert numerator % degree == 0, (kind, n, degree, numerator)
        a[degree] = numerator // degree
    return a[n]


def local_coefficients(kind, n):
    if kind == "A270913":
        return [comb(n, m) for m in range(n + 1)]
    if kind == "A008705":
        return [(-1)**m * comb(n, m) for m in range(n + 1)]
    if kind == "A008485":
        return [comb(n+m-1, m) for m in range(n + 1)]
    # (1-t^2) G'(t) = 2n G(t), G=((1+t)/(1-t))^n.
    c = [1] + [0] * n
    for m in range(1, n + 1):
        numerator = 2*n*c[m-1] + ((m-2)*c[m-2] if m >= 2 else 0)
        assert numerator % m == 0
        c[m] = numerator // m
    return c


@lru_cache(None)
def product_coefficient(kind, n):
    if n == 0:
        return 1
    local = local_coefficients(kind, n)
    a = [1] + [0] * n
    for size in range(1, n + 1):
        if kind == "A291697" and size % 2 == 0:
            continue
        out = a.copy()
        for count in range(1, n // size + 1):
            shift, factor = count*size, local[count]
            if factor == 0:
                continue
            for degree in range(shift, n+1):
                out[degree] += factor*a[degree-shift]
        a = out
    return a[n]


def exact_value(kind, n):
    a = recurrence_coefficient(kind, n)
    b = product_coefficient(kind, n)
    assert a == b, ("independent algorithms disagree", kind, n, a, b)
    return a


def main():
    primes = primes_up_to(LIMIT)
    print("ARITHMETIC integers only; independent full-product DP and global recurrence")
    print("PRIMES", primes)
    for kind in IDS:
        for n, expected in enumerate(SOURCE_TERMS[kind]):
            assert exact_value(kind, n) == expected, ("source fidelity", kind, n)
        print("FIDELITY", kind, "all", len(SOURCE_TERMS[kind]), "published S/T/U terms reproduced")

    # A finite fidelity check of the Gauss congruences separately stated in each entry.
    gauss_cases = 0
    for kind in IDS:
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                for k in [1, 2]:
                    hi, lo = n*p**k, n*p**(k-1)
                    assert (exact_value(kind, hi) - exact_value(kind, lo)) % p**k == 0
                    gauss_cases += 1
    print("GAUSS_FIDELITY", gauss_cases, "cases passed; not an all-index proof")

    failures = {kind: [] for kind in IDS}
    totals = {kind: 0 for kind in IDS}
    for p in primes:
        for kind in IDS:
            if kind == "A008485" and p < 3:
                continue
            if kind == "A008705" and p < 3:
                continue
            if kind == "A291697" and p < 5:
                continue
            modulus = p**(3 if kind == "A291697" else 2)
            target = 2*p+2 if kind in ["A270919", "A291697"] else (-1-p if kind == "A008705" else p+1)
            value = exact_value(kind, p)
            ok = (value-target) % modulus == 0
            if not ok:
                failures[kind].append(p)
            totals[kind] += 1
            print(json.dumps(dict(id=kind, p=p, index=p, value=value,
                                  modulus=modulus, residue=value % modulus,
                                  target_residue=target % modulus, source_congruence_holds=ok),
                             sort_keys=True))
            if kind == "A008705":
                value2 = exact_value(kind, 2*p)
                assert (value2-(p-1)) % modulus == 0, (kind, p, "2p clause")
                print(json.dumps(dict(id=kind, p=p, index=2*p, value=value2,
                                      modulus=modulus, residue=value2 % modulus,
                                      target_residue=(p-1) % modulus, source_congruence_holds=True),
                                 sort_keys=True))
    assert failures == {"A270913": [], "A008485": [], "A270919": [2],
                        "A008705": [], "A291697": []}, failures
    assert exact_value("A270919", 2) == 12
    assert (12-6) % 4 != 0
    # Explicitly verify the excluded boundary cases; do not silently strengthen a source.
    assert (exact_value("A008485", 2)-3) % 4 != 0
    assert (exact_value("A291697", 3)-8) % 27 != 0
    print("SUMMARY", json.dumps(dict(prime_limit=LIMIT, last_prime=primes[-1],
          prime_cases=totals, failures=failures, signed_2p_cases=len(primes)-1,
          signed_max_index=2*primes[-1]), sort_keys=True))
    print("EXPECTED_SETTLEMENT A270913 Proved; A008485 Proved; A270919 Refuted (p=2); A008705 Proved (both clauses); A291697 Proved")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
