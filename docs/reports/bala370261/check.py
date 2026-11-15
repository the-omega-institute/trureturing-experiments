#!/usr/bin/env python3
"""Exact finite fidelity and identity checks for OEIS A370261; no universal proof."""
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def cheb_t(m, x):
    """T_0(x)=1, T_1(x)=x, T_(m+2)(x)=2*x*T_(m+1)(x)-T_m(x)."""
    require(m >= 0, "T index")
    previous, current = 1, x
    if m == 0:
        return previous
    for _ in range(1, m):
        previous, current = current, 2 * x * current - previous
    return current


def cheb_u(m, x):
    """U_0(x)=1, U_1(x)=2*x, with the same recurrence."""
    require(m >= 0, "U index")
    previous, current = 1, 2 * x
    if m == 0:
        return previous
    for _ in range(1, m):
        previous, current = current, 2 * x * current - previous
    return current


def matrix_mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def cheb_t_matrix(m, x):
    """Independent binary matrix-power evaluation of the first-kind recurrence."""
    if m == 0:
        return 1
    accumulator = ((1, 0), (0, 1))
    power = ((2 * x, -1), (1, 0))
    exponent = m - 1
    while exponent:
        if exponent % 2:
            accumulator = matrix_mul(accumulator, power)
        power = matrix_mul(power, power)
        exponent //= 2
    return accumulator[0][0] * x + accumulator[0][1]


def read_bfile(name):
    rows = {}
    for line in (ROOT / name).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            n, a = map(int, line.split())
            require(n not in rows, f"duplicate index in {name}")
            rows[n] = a
    return rows


def main():
    require(sys.version_info >= (3, 12), "Python 3.12 or newer required")
    sources = json.loads((ROOT / "source-fields.json").read_text())
    target = read_bfile("b370261.txt")
    parent = read_bfile("b370259.txt")
    require(sources["A370261"]["offset"] == "1,2", "source offset")
    require(sources["A370261"]["revision"] == 15, "source revision")
    require(sources["A370261"]["comment"] ==
            ["The sequence is conjectured to be integral."], "source conjecture")
    require(sorted(target) == list(range(1, 301)), "target b-file range")
    require(sorted(parent) == list(range(1, 351)), "parent b-file range")
    listed = [int(v) for v in sources["A370261"]["data"].split(",")]
    require(listed == [target[n] for n in range(1, len(listed) + 1)],
            "published target prefix")
    for n, a in parent.items():
        require(Fraction(cheb_t(n, n + 1) - 1, n ** 3) == a,
                f"parent fidelity n={n}")
    print("SOURCE_FIDELITY A370259: all 350 published b-file rows agree")
    print(f"SOURCE_FIDELITY A370261: offset 1,2; revision 15; {len(listed)} listed terms agree")

    witnesses = []
    for n in range(1, 301):
        x = 2 * n + 1
        tn, t2n, u = cheb_t(n, x), cheb_t(2 * n, x), cheb_u(n - 1, x)
        require(cheb_t_matrix(2 * n, x) == t2n, f"matrix cross-check n={n}")
        require(t2n - 1 == 2 * (tn * tn - 1), f"doubling n={n}")
        require(tn * tn - (x * x - 1) * u * u == 1, f"Pell identity n={n}")
        require(cheb_u(n - 1, 1) == n, f"U at one n={n}")
        require((u - n) % (2 * n) == 0, f"polynomial congruence n={n}")
        require(u % n == 0, f"integer quotient n={n}")
        q = u // n
        require(q % 2 == 1, f"odd witness n={n}")
        denominator = (n + 1) * (2 * n) ** 3
        require(denominator > 0, f"positive denominator n={n}")
        radicand = Fraction(t2n - 1, denominator)
        name_radicand = Fraction(cheb_t(2 * n, 2 * n + 1) - 1, (2 * n) ** 3) / (n + 1)
        require(radicand == name_radicand == q * q, f"normalized square n={n}")
        require(radicand.denominator == 1 and radicand.numerator >= 0,
                f"real square-root domain n={n}")
        r = isqrt(radicand.numerator)
        require(r * r == radicand.numerator and r == abs(q) == target[n],
                f"source b-file agreement n={n}")
        witnesses.append(q)
    print("TARGET n=1..300: exact T matrix cross-check, doubling, Pell, congruence, square and b-file agreement")
    print("FIRST_WITNESSES", ",".join(map(str, witnesses[:5])))
    n, x = 4, 9
    print("EXAMPLE n=4: U_3(9)=" + str(cheb_u(3, x)) +
          "; T_8(9)-1=" + str(cheb_t(8, x) - 1) + "; denominator=2560; q=1449")

    count = 0
    for n in range(1, 101):
        for c in range(-7, 8):
            x = 1 + c * n
            u = cheb_u(n - 1, x)
            require(u % n == 0, f"general by-product n={n}, c={c}")
            if c:
                require((u - n) % abs(c * n) == 0,
                        f"stronger polynomial congruence n={n}, c={c}")
            count += 1
    print(f"BYPRODUCT {count} cases: n=1..100, c=-7..7; includes odd, zero and negative c")

    count = 0
    for x in range(-4, 5):
        # v_m=U_(m-1)(x), m>=1, with v_0=0 and v_1=1.
        v, next_v = 0, 1
        for m in range(21):
            require(next_v ** 2 - 2 * x * next_v * v + v ** 2 == 1,
                    f"invariant x={x}, m={m}")
            require(cheb_t(m, x) == next_v - x * v, f"T bridge x={x}, m={m}")
            v, next_v = next_v, 2 * x * next_v - v
            count += 1
    print(f"INVARIANT {count} cases: x=-4..4, m=0..20")
    for name in ["source-fields.json", "b370261.txt", "b370259.txt"]:
        print("INPUT_SHA256", name, sha256((ROOT / name).read_bytes()).hexdigest())
    print("Finite checks support the preregistered argument; they are not a universal or Lean proof.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
