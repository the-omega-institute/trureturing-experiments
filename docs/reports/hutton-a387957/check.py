#!/usr/bin/env python3
"""Exact finite checks of Thomas Hutton's OEIS A387957, revision 20.

Each (size, color) has an unbounded nonnegative multiplicity. Equal size/color
parts are indistinguishable. Supplemental checks do not prove the all-n claim.
"""

LIMIT = 128
# Sean A. Irvine's OEIS A387957 b-file, n=0..56.
# https://oeis.org/A387957/b387957.txt; OEIS attribution: CC BY-SA 4.0.
OEIS_TERMS = (0, 1, 3, 10, 25, 63, 141, 310, 640, 1295, 2518, 4802, 8907, 16240, 28995, 50996, 88208, 150587, 253570, 422100, 694509, 1131190, 1824000, 2914730, 4616770, 7253925, 11308411, 17501344, 26895795, 41061345, 62289666, 93926230, 140811304, 209939415, 311343865, 459390240, 674519798, 985753748, 1434075585, 2077219340, 2996137664, 4304048329, 6158640195, 8778971080, 12468220150, 17644910187, 24884835175, 34978155576, 49005836080, 68443028815, 95297024980, 132292699854, 183119118315, 252759207575, 347926631420, 477646862112, 654024518194)


def residue_counts(limit):
    """Unbounded-part DP with four color-sum residues, not the product formula."""
    counts = [[0] * 4 for _ in range(limit + 1)]
    counts[0][0] = 1
    for size in range(1, limit + 1):
        for color in range(4):
            for total in range(size, limit + 1):
                for residue in range(4):
                    counts[total][(residue + color) % 4] += counts[total - size][residue]
    return counts


def partitions_by_pentagonal_recurrence(limit):
    """Euler's exact recurrence; independent of the colored-multiplicity DP."""
    p = [0] * (limit + 1)
    p[0] = 1
    for n in range(1, limit + 1):
        k = 1
        while k * (3 * k - 1) // 2 <= n:
            sign = 1 if k % 2 else -1
            for offset in (k * (3 * k - 1) // 2, k * (3 * k + 1) // 2):
                if offset <= n:
                    p[n] += sign * p[n - offset]
            k += 1
    return p


def convolution(a, b, limit):
    return [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(limit + 1)]


def colored_multisets(n, colors=4):
    """Literal recursive enumeration, used only for n<=8 and source fidelity."""
    modes = [(size, color) for size in range(1, n + 1) for color in range(colors)]

    def enumerate_at(index, left, prefix):
        if index == len(modes):
            if left == 0:
                yield prefix
            return
        size, color = modes[index]
        for multiplicity in range(left // size + 1):
            yield from enumerate_at(index + 1, left - multiplicity * size,
                                    prefix + ((size, color),) * multiplicity)

    return enumerate_at(0, n, ())


def main():
    assert LIMIT >= len(OEIS_TERMS) - 1
    counts = residue_counts(LIMIT)
    p = partitions_by_pentagonal_recurrence(LIMIT)
    p2 = convolution(p, p, LIMIT)
    p4 = convolution(p2, p2, LIMIT)

    expected_n2 = {((2, 1),), ((1, 0), (1, 1)), ((1, 2), (1, 3))}
    actual_n2 = {m for m in colored_multisets(2) if sum(c for _, c in m) % 4 == 1}
    assert actual_n2 == expected_n2, (actual_n2, expected_n2)
    print('FIDELITY: all three source n=2 residue-one multisets reproduced exactly')
    print('n=2 objects:', sorted(actual_n2))

    for n in range(9):
        literal = [0] * 4
        for m in colored_multisets(n):
            literal[sum(c for _, c in m) % 4] += 1
        assert literal == counts[n], (n, literal, counts[n])
        plain = sum(1 for _ in colored_multisets(n, colors=1))
        assert plain == p[n], (n, plain, p[n])
    print('INDEPENDENCE: literal multiset enumeration matches all four DP residues for n=0..8')
    print('INDEPENDENCE: literal uncolored enumeration matches Euler recurrence for n=0..8')

    for n, expected in enumerate(OEIS_TERMS):
        assert counts[n][1] == expected, (n, counts[n][1], expected)
    print(f'FIDELITY: all {len(OEIS_TERMS)} OEIS b-file values n=0..{len(OEIS_TERMS)-1} matched')

    print('n,a(n),p4(n),even_correction')
    for n in range(LIMIT + 1):
        correction = p2[n // 2] if n % 2 == 0 else 0
        assert counts[n][1] == counts[n][3], (n, counts[n])
        assert sum(counts[n]) == p4[n], (n, counts[n], p4[n])
        assert counts[n][0] - counts[n][1] + counts[n][2] - counts[n][3] == correction
        assert 4 * counts[n][1] == p4[n] - correction, (n, counts[n], p4[n], correction)
        print(f'{n},{counts[n][1]},{p4[n]},{correction}')
    print(f'IDENTITY: 4*a(n)=p4(n)-even_correction for every n=0..{LIMIT}; integer arithmetic only')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
