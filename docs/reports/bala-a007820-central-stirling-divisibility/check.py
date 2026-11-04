#!/usr/bin/env python3
"""Exact finite checks for OEIS A007820 and Bala's odd-index conjecture."""

from math import comb, factorial


MAX_M = 201
SOURCE_DATA = (
    1, 1, 7, 90, 1701, 42525, 1323652, 49329280, 2141764053,
    106175395755, 5917584964655, 366282500870286,
    24930204590758260, 1850568574253550060, 148782988064375309400,
    12879868072770626040000, 1194461517469807833782085,
    118144018577011378596484455,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def stirling_table(max_row, max_column):
    """S(0,0)=1; S(r+1,k)=k*S(r,k)+S(r,k-1), with zero boundaries."""
    rows = [[0] * (max_column + 1) for _ in range(max_row + 1)]
    rows[0][0] = 1
    for r in range(1, max_row + 1):
        for k in range(1, min(r, max_column) + 1):
            rows[r][k] = k * rows[r - 1][k] + rows[r - 1][k - 1]
    return rows


def alternating_power_sum(r, m):
    """Python integer powers use the inclusion-exclusion convention 0**0=1."""
    return sum((-1) ** (m - i) * comb(m, i) * i ** r
               for i in range(m + 1))


def main():
    s = stirling_table(2 * MAX_M, MAX_M)

    # Fidelity: the complete displayed DATA block of A007820 revision 89.
    for m, expected in enumerate(SOURCE_DATA):
        require(s[2 * m][m] == expected, f"source DATA mismatch at m={m}")
    print(f"FIDELITY_DATA_OK: m=0..{len(SOURCE_DATA) - 1}; {len(SOURCE_DATA)} values")

    # A second, independent formula checks the recurrence table, including
    # r<m (zero) and r=0. All operations are on exact Python integers.
    ie_count = 0
    for m in range(26):
        for r in range(2 * m + 1):
            require(alternating_power_sum(r, m) == factorial(m) * s[r][m],
                    f"inclusion-exclusion mismatch at (r,m)=({r},{m})")
            ie_count += 1
    print(f"INCLUSION_EXCLUSION_OK: m=0..25, r=0..2*m; {ie_count} pairs")

    # Fidelity: Kruchinin's FORMULA (2013), already present in the source.
    for m in range(MAX_M + 1):
        reflected = sum(
            (-1) ** (m - j) * m ** j * comb(2 * m, j) * s[2 * m - j][m]
            for j in range(m + 1)
        )
        require(s[2 * m][m] == reflected, f"reflection mismatch at m={m}")
    print(f"FIDELITY_REFLECTION_OK: m=0..{MAX_M}; {MAX_M + 1} values")

    odd_count = 0
    for m in range(1, MAX_M + 1, 2):
        require(s[2 * m][m] % (m * m) == 0,
                f"Bala square-divisibility failure at m={m}")
        # Check the modular extraction independently of the signed source sum.
        tail = sum(
            (-1) ** (m + j) * comb(2 * m, j) * m ** j * s[2 * m - j][m]
            for j in range(2, m + 1)
        )
        require(2 * s[2 * m][m] == 2 * m * m * s[2 * m - 1][m] + tail,
                f"j=0,1 extraction mismatch at m={m}")
        require(tail % (m * m) == 0, f"tail divisibility failure at m={m}")
        odd_count += 1
    print(f"ODD_SQUARE_AND_EXTRACTION_OK: m=1..{MAX_M}, odd; {odd_count} values")

    for m in (1, 3, 5, 7, 9):
        print(f"EXAMPLE: m={m}, S(2*m,m)={s[2*m][m]}, quotient={s[2*m][m] // (m*m)}")

    require(s[4][2] == 7 and s[4][2] % 4 != 0, "even-index boundary mismatch")
    print("EVEN_BOUNDARY_OK: m=2, S(4,2)=7 is not divisible by 4")
    print("SCOPE: finite exact checks only; the universal assertion is not proved here")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
