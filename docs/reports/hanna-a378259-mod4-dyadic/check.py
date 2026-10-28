#!/usr/bin/env python3
"""Exact finite checks for Paul D. Hanna's OEIS A378259 conjecture.

The original rational equation is reconstructed before any residue tests.
This program gives finite evidence; it does not prove the unbounded claim.
"""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

N = 500
SCALED_DEGREE = 60
PUBLISHED_PREFIX = (
    1, 2, 14, 32, 286, 480, 3532, 7520, 75126, 138560, 897876,
    1997408, 17039340, 29800896, 233106072, 497063872, 5068814502,
    10555899520, 67510782340, 148762640480, 961465207492,
    1640373169344, 14064469534248, 31207376374080, 312901302489020,
    640840797612416, 4625855789637960, 9196764841428416,
)


def original_coefficients(nterms):
    """Solve A=sum_m Q(x^m), A=x U, Q=x U(x^2)/U, over Q.

    q[n] + sum_{j=1}^{n-1} u[j]q[n-j] = [x^(n-1)]U(x^2).
    Also u[n-1]=sum_{d|n}q[d], with q[1]=u[0]=1.
    Removing q[n] from these equations leaves 2*u[n-1] known.
    """
    u = [Fraction(1)]
    q = [Fraction(0), Fraction(1)]
    for n in range(2, nterms + 1):
        numerator = u[(n - 1) // 2] if n % 2 else Fraction(0)
        cross = sum((u[j] * q[n - j] for j in range(1, n - 1)), Fraction(0))
        proper = sum((q[d] for d in range(1, n) if n % d == 0), Fraction(0))
        u.append((numerator - cross + proper) / 2)
        q.append(u[-1] - proper)
    return u, q


def product(left, right, degree):
    out = [0] * (degree + 1)
    for i, a in enumerate(left[:degree + 1]):
        if a:
            for j, b in enumerate(right[:degree + 1 - i]):
                if b:
                    out[i + j] += a * b
    return out


def substitute(series, scalar, power, degree):
    out = [0] * (degree + 1)
    for j in range(min(len(series) - 1, degree // power) + 1):
        out[power * j] = series[j] * scalar ** j
    return out


def inverse_unit(series, degree):
    assert series[0] == 1
    out = [1] + [0] * degree
    for n in range(1, degree + 1):
        out[n] = -sum(series[i] * out[n - i] for i in range(1, n + 1))
    return out


def dyadic(n):
    return n > 0 and (n & (n - 1)) == 0


def main():
    fixture = Path(__file__).with_name("source-b378259.txt")
    source_data = fixture.read_bytes()
    table = {}
    for line in source_data.decode().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            n, value = map(int, line.split())
            assert n not in table
            table[n] = value
    assert set(table) == set(range(1, N + 1))

    u, q = original_coefficients(N)
    assert q[0] == 0 and q[1] == u[0] == 1
    for n in range(1, N + 1):
        # Original NAME equation, summed coefficientwise, not numerically.
        assert u[n - 1] == sum(q[d] for d in range(1, n + 1) if n % d == 0)
        # Independent product check of the cancelled quotient Q=A(x^2)/A(x).
        assert sum(u[j] * q[n - j] for j in range(n)) == (
            u[(n - 1) // 2] if n % 2 else 0
        )
    terms = []
    for n in range(1, N + 1):
        scaled = 4 ** (n - 1) * u[n - 1]
        assert scaled.denominator == 1, (n, scaled)
        terms.append(scaled.numerator)
    assert tuple(terms[:len(PUBLISHED_PREFIX)]) == PUBLISHED_PREFIX
    assert all(terms[n - 1] == table[n] for n in range(1, N + 1))
    assert next(n for n, value in enumerate(terms, 1) if value < 0) == 224
    assert terms[0] == 1
    assert all(value % 2 == 0 for value in terms[1:])

    # b(x)=A(4x)/(4x); b[r]=a(r+1).
    b = terms
    degree = SCALED_DEGREE
    tail = [0] * (degree + 1)
    for m in range(2, degree + 2):
        residual_degree = degree - (m - 1)
        numerator = substitute(b, 4 ** (2 * m - 1), 2 * m, residual_degree)
        denominator = substitute(b, 4 ** (m - 1), m, residual_degree)
        ratio = product(numerator, inverse_unit(denominator, residual_degree), residual_degree)
        for j, value in enumerate(ratio):
            tail[j + m - 1] += 4 ** (m - 1) * value
    rhs = substitute(b, 4, 2, degree)
    btail = product(b, tail, degree)
    assert product(b, b, degree) == [a + c for a, c in zip(rhs, btail)]
    assert [value % 8 for value in tail] == [0, 4] + [0] * (degree - 1)

    square = product(b, b, N - 1)
    assert [value % 8 for value in square] == [1, 4] + [0] * (N - 2)
    c = [0] + [value // 2 for value in b[1:]]
    csquare = product(c, c, N - 1)
    assert [(c[r] + csquare[r]) % 2 for r in range(N)] == [0, 1] + [0] * (N - 2)
    assert all(c[r] % 2 == int(dyadic(r)) for r in range(1, N))
    support = []
    for n in range(2, N + 1):
        expected = 2 if dyadic(n - 1) else 0
        assert terms[n - 1] % 4 == expected, (n, terms[n - 1], expected)
        if expected == 2:
            support.append(n)

    print("arithmetic: Python integers and fractions.Fraction only")
    print(f"source table: {fixture.name}; sha256={sha256(source_data).hexdigest()}")
    print(f"fidelity: all {N} independently reconstructed terms match the OEIS table")
    print("first ten: " + ",".join(map(str, terms[:10])))
    print(f"source NAME equation and quotient product: degrees 1..{N}, exact")
    print("source first negative term: a(224)=" + str(terms[223]))
    print(f"scaled fixed-point square equation: degrees 0..{degree}, exact")
    print(f"integrality: n=1..{N}; evenness: n=2..{N}")
    print(f"b^2=1+4x modulo 8: degrees 0..{N - 1}")
    print(f"c+c^2=x modulo 2 and dyadic support: degrees 0..{N - 1}")
    print(f"conjecture: all {N - 1} cases n=2..{N}; residue-2 support={support}")
    print("scope: finite exact checks; the unbounded conjecture is not proved here")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
