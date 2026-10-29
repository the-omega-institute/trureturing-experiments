#!/usr/bin/env python3
"""Exact EGF fidelity and bounded checks for Lang's OEIS A130187 conjecture.

Python 3 standard library only. This is finite arithmetic evidence, not a proof.
The EGF route uses the generic formal identity E'=g'E for E=exp(g),
with g=X^2/4. The candidate route independently multiplies odd integers.
"""
from fractions import Fraction
from hashlib import sha256
from math import factorial, gcd
from pathlib import Path

MAX_N = 1000
SOURCE_R = ['1', '3/4', '5/4', '105/32', '189/16', '3465/64', '19305/64', '2027025/1024', '3828825/256', '130945815/1024', '1249937325/1024', '105411381075/8192', '608142583125/4096', '30494006668125/16384', '412685556908625/16384', '191898783962510625/524288', '372509404162520625/65536', '24627010608522196875/262144', '431620764875678503125/262144', '63966197354575554163125/2097152', '624431926556570885878125/1048576', '51260184516416682722540625/4194304', '1103208318940272084680765625/4194304', '397522730924811374513302546875/67108864', '2337433657837890882138218975625/16777216', '229248301057177759594325322609375/67108864', '5850077015866499124462598047328125/67108864', '1241052052651678742832422585754609375/536870912', '17075164448552407530694365921244453125/268435456', '1947707091431544619001204006083283953125/1073741824']
SOURCE_Z = ['0', '-1/2', '0', '-3/8', '0', '-5/8', '0', '-105/64', '0', '-189/32', '0', '-3465/128', '0', '-19305/128', '0', '-2027025/2048', '0', '-3828825/512', '0', '-130945815/2048', '0', '-1249937325/2048', '0', '-105411381075/16384', '0', '-608142583125/8192', '0', '-30494006668125/32768', '0', '-412685556908625/32768', '0', '-191898783962510625/1048576', '0', '-372509404162520625/131072', '0', '-24627010608522196875/524288', '0', '-431620764875678503125/524288', '0', '-63966197354575554163125/4194304']
BFILE_SHA256 = '513e2f71231546b181fd885e9619ac665c6e02c2053e66c9bb477b546bf19a36'


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def oddpart(m):
    require(m > 0, "oddpart domain")
    while m % 2 == 0:
        m //= 2
    return m


def main():
    # m*[X^m]exp(g) = sum_k k*[X^k]g * [X^(m-k)]exp(g).
    degree = 2 * MAX_N + 2
    g = {2: Fraction(1, 4)}
    e = [Fraction(1)]
    for m in range(1, degree + 1):
        e.append(sum((k * c * e[m-k] for k, c in g.items() if k <= m), Fraction(0)) / m)
    z = [factorial(j) * (-2 * e[j+1]) for j in range(degree)]
    r = [-2 * z[2*n+1] for n in range(MAX_N+1)]
    require(r[:len(SOURCE_R)] == list(map(Fraction, SOURCE_R)), "source r table")
    require(z[:len(SOURCE_Z)] == list(map(Fraction, SOURCE_Z)), "source z table")
    require(all(z[j] == 0 for j in range(0, degree, 2)), "even z entries")
    print(f"Fidelity: all {len(SOURCE_R)} displayed r values and {len(SOURCE_Z)} displayed z values agree with the formal EGF")
    print("Source label z(n), n=0..40 contains 40 displayed entries (indices 0..39); EGF also gives z(40)=0")

    numerators = []
    double_factorial = 1
    for n in range(MAX_N+1):
        double_factorial *= 2*n+1
        u = oddpart(n+1)
        t = 0
        m = n+1
        while m % 2 == 0:
            m //= 2
            t += 1
        require(double_factorial % u == 0, f"odd factor divisibility at {n}")
        expected_num = double_factorial // u
        expected_den = 2 ** (n+t)
        require(expected_num % 2 == 1, f"odd quotient at {n}")
        require(gcd(expected_num, expected_den) == 1, f"coprimality at {n}")
        require(r[n].numerator == expected_num, f"numerator at {n}")
        require(r[n].denominator == expected_den, f"denominator at {n}")
        require(r[n] == Fraction(factorial(2*n+1), 4**n * factorial(n+1)), f"coefficient extraction at {n}")
        require(r[n] == Fraction(double_factorial, 2**n*(n+1)), f"factorial separation at {n}")
        require(expected_num == double_factorial // gcd(n+1, double_factorial), f"disclosed LODA formula at {n}")
        numerators.append(expected_num)
    print(f"Settlement checks: all {MAX_N+1} indices 0..{MAX_N}; numerator, denominator, coefficient bridge, cancellation, and LODA gcd formula; zero discrepancies")
    path = Path(__file__).with_name("b130187.txt")
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == BFILE_SHA256, "b-file snapshot hash")
    rows = []
    for line in raw.decode().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            n, value = map(int, line.split())
            rows.append((n, value))
    require([n for n, _ in rows] == list(range(404)), "published b-file index range")
    require(all(numerators[n] == value for n, value in rows), "published b-file values")
    print(f"Published b-file: all {len(rows)} rows 0..403 agree; sha256={BFILE_SHA256}")
    print("First five rationals: " + ", ".join(map(str, r[:5])))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
