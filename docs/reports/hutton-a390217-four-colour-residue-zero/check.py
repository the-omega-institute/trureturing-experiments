#!/usr/bin/env python3
"""Exact independent checks of Thomas Hutton's OEIS A390217 conjectures.

Source: https://oeis.org/A390217, revision 29, and Sean A. Irvine's
https://oeis.org/A390217/b390217.txt (n = 0..56).
The embedded source data are attributed to OEIS contributors under CC BY-SA 4.0.
This finite check is supplemental evidence, not a proof for every n.
"""

import json
from collections import Counter

LIMIT = 128
SOURCE_TERMS = (1, 1, 4, 10, 28, 63, 146, 310, 651, 1295, 2536, 4802, 8941, 16240, 29050, 50996, 88303, 150587, 253720, 422100, 694753, 1131190, 1824376, 2914730, 4617358, 7253925, 11309296, 17501344, 26897135, 41061345, 62291644, 93926230, 140814226, 209939415, 311348100, 459390240, 674525928, 985753748, 1434084330, 2077219340, 2996150106, 4304048329, 6158657696, 8778971080, 12468244683, 17644910187, 24884869250, 34978155576, 49005883236, 68443028815, 95297089736, 132292699854, 183119206909, 252759207575, 347926751840, 477646862112, 654024681269)


def residue_dp(limit):
    """Unbounded multiplicities for independent (positive size, color) modes."""
    counts = [[0] * 4 for _ in range(limit + 1)]
    counts[0][0] = 1
    for size in range(1, limit + 1):
        for color in range(4):
            for energy in range(size, limit + 1):
                for residue in range(4):
                    counts[energy][residue] += counts[energy - size][(residue - color) % 4]
    return counts


def ordinary_partitions(limit):
    """Euler's pentagonal recurrence, independent of the residue DP."""
    p = [1] + [0] * limit
    for n in range(1, limit + 1):
        k = 1
        while k * (3 * k - 1) // 2 <= n:
            sign = 1 if k % 2 else -1
            for degree in (k * (3 * k - 1) // 2, k * (3 * k + 1) // 2):
                if degree <= n:
                    p[n] += sign * p[n - degree]
            k += 1
    return p


def convolution(a, b, limit):
    return [sum(a[k] * b[n - k] for k in range(n + 1)) for n in range(limit + 1)]


def literal_multisets(n):
    """Yield unordered colored multisets by choosing each cell's multiplicity."""
    cells = [(size, color) for size in range(1, n + 1) for color in range(4)]

    def visit(index, remaining, chosen):
        if remaining == 0:
            yield tuple(chosen)
            return
        if index == len(cells):
            return
        size, color = cells[index]
        for copies in range(remaining // size + 1):
            yield from visit(index + 1, remaining - copies * size, chosen + [(size, color)] * copies)

    yield from visit(0, n, [])


def main():
    residues = residue_dp(LIMIT)
    p = ordinary_partitions(LIMIT)
    p2 = convolution(p, p, LIMIT)
    p4 = convolution(p2, p2, LIMIT)

    for n, value in enumerate(SOURCE_TERMS):
        assert residues[n][0] == value, ("OEIS b-file mismatch", n, residues[n][0], value)
    assert len(SOURCE_TERMS) == 57
    assert p[:9] == [1, 1, 2, 3, 5, 7, 11, 15, 22]

    expected_n2 = {
        ((2, 0),),
        ((1, 0), (1, 0)),
        ((1, 1), (1, 3)),
        ((1, 2), (1, 2)),
    }
    actual_n2 = {m for m in literal_multisets(2) if sum(c for _, c in m) % 4 == 0}
    assert actual_n2 == expected_n2, ("source example", actual_n2)
    for n in range(9):
        count = Counter(sum(c for _, c in m) % 4 for m in literal_multisets(n))
        assert [count[r] for r in range(4)] == residues[n], ("literal enumeration", n)

    rows = []
    for n in range(LIMIT + 1):
        a0, a1, a2, a3 = residues[n]
        even = p2[n // 2] if n % 2 == 0 else 0
        fourth = p[n // 4] if n % 4 == 0 else 0
        assert a0 + a1 + a2 + a3 == p4[n], ("unfiltered", n)
        assert a0 - a1 + a2 - a3 == even, ("half-turn", n)
        assert a1 - a3 == 0, ("quarter-turn imaginary part", n)
        assert a0 - a2 == fourth, ("quarter-turn real part", n)
        assert 4 * a0 == p4[n] + even + 2 * fourth, ("claimed identity", n)
        if n <= 8 or n in (16, 32, 56, 64, 100, 128):
            rows.append({"n": n, "residues": residues[n], "p4": p4[n], "even_p2": even, "fourth_p": fourth})

    print(json.dumps({
        "source": "OEIS A390217 revision 29; Irvine b-file n=0..56",
        "source_terms_matched": len(SOURCE_TERMS),
        "displayed_terms_matched": 35,
        "literal_enumeration_interval": [0, 8],
        "identity_interval": [0, LIMIT],
        "arithmetic": "Python integers only",
        "source_n2_multisets": [list(m) for m in sorted(actual_n2)],
        "rows": rows,
        "all_n_proof": False,
    }, indent=2, sort_keys=True))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
