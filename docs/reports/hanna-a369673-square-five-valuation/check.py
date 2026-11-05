#!/usr/bin/env python3
"""Exact finite checks of Hanna's OEIS A369673 conjectures C.1--C.4.

Only Python's integer arithmetic is used. These checks do not prove the
universal conjectures. The published terms are fidelity checks independent
of the proposed pairing and factorial identities.
"""

from math import factorial, isqrt, prod


def valuation(p: int, value: int) -> int:
    if p < 2 or value <= 0:
        raise ValueError("valuation requires p >= 2 and a positive integer")
    exponent = 0
    while value % p == 0:
        value //= p
        exponent += 1
    return exponent


def a(index: int) -> int:
    return prod(2**k + 2 ** (index - k) for k in range(index + 1))


def b(index: int) -> int:
    return valuation(5, factorial(10 * index)) - valuation(5, factorial(5 * index))


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def main() -> None:
    # OEIS A369673 revision #21, %S/%T/%U (all 15 displayed terms).
    published_a = [
        2, 9, 100, 2916, 231200, 50808384, 31258240000,
        54112148361216, 264265663201280000, 3645603832850650497024,
        142153785549232537600000000, 15673043740102659990892604030976,
        4886752115388739132874502963200000000,
        4309225323078788454199311474023086952546304,
        10747393363422494556085100202291563069440000000000,
    ]
    # OEIS A127428 revision #21, the first 31 displayed terms.
    published_b = [
        0, 1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 17, 18, 19,
        20, 21, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 35, 36, 37,
    ]
    require([a(n) for n in range(15)] == published_a, "A369673 source fidelity")
    require([b(n) for n in range(31)] == published_b, "A127428 source fidelity")
    print("FIDELITY A369673 N=0..14: 15/15 published terms matched")
    print("FIDELITY A127428 n=0..30: 31/31 published terms matched")

    values = [a(N) for N in range(401)]
    for N, value in enumerate(values):
        require(value > 0, f"positivity at N={N}")
        root = isqrt(value)
        require((root * root == value) == (N % 4 != 0), f"C.1 at N={N}")
        if N % 2:
            m = N // 2
            half = prod(2**k + 2 ** (N - k) for k in range(m + 1))
            require(value == half * half, f"odd pairing at N={N}")
            require(value % 5 != 0, f"C.2 at n={m}")
        else:
            m = N // 2
            Q = prod(1 + 4**j for j in range(1, m + 1))
            require(Q % 2 == 1, f"odd Q at m={m}")
            require(value == 2 ** (m * m + 1) * Q * Q, f"even pairing at N={N}")
            require(valuation(2, value) == m * m + 1, f"two-adic order at N={N}")
            if N % 4 == 0:
                require(valuation(5, value) == 2 * b(N // 4), f"C.3 at n={N // 4}")
            else:
                require(valuation(5, value) == 2 * b(N // 4 + 1), f"C.4 at n={N // 4}")
    print("C.1 N=0..400: 401/401 exact integer-square checks passed")
    print("C.2 n=0..199: 200/200 nondivisibility checks passed")
    print("C.3 n=0..100: 101/101 direct factorial-valuation checks passed")
    print("C.4 n=0..99: 100/100 direct factorial-valuation checks passed")
    print("PAIRING even m=0..200 and odd m=0..199: 401/401 identities passed")

    for j in range(1, 5001):
        expected = 0 if j % 2 == 0 else 1 + valuation(5, j)
        require(valuation(5, 4**j + 1) == expected, f"lifting at j={j}")
    print("LIFTING j=1..5000: 5000/5000 exact valuations passed")

    for r in range(202):
        odd_sum = sum(1 + valuation(5, j) for j in range(1, 2 * r + 1, 2))
        factorial_sum = r + valuation(5, factorial(2 * r)) - valuation(5, factorial(r))
        require(b(r) == odd_sum == factorial_sum, f"factorial bridge at r={r}")
        require(valuation(5, factorial(5 * r)) == r + valuation(5, factorial(r)),
                f"factorial step at r={r}")
    print("FACTORIAL BRIDGE r=0..201: 202/202 identities and factorial steps passed")
    require(a(0) == 2 and b(0) == 0, "zero-index conventions")
    print("EDGE a(0)=2; b(0)=0; 0!=1; empty Q product=1")
    for N in (0, 1, 2, 4, 6):
        print(f"SAMPLE N={N} a(N)={values[N]} v5={valuation(5, values[N])}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
