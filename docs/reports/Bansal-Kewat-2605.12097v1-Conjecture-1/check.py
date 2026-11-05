#!/usr/bin/env python3
"""Exact GF(2) checks for Bansal–Kewat, arXiv:2605.12097v1, Conjecture 1.

Bit i of an integer is the coefficient of x^i. Addition is XOR;
multiplication is carry-less. No floating-point arithmetic is used.
"""


def mul(p, q):
    result = 0
    while q:
        if q & 1:
            result ^= p
        p <<= 1
        q >>= 1
    return result


def power(p, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, p)
        p = mul(p, p)
        exponent >>= 1
    return result


def remainder(p, q):
    assert q
    while p.bit_length() >= q.bit_length():
        p ^= q << (p.bit_length() - q.bit_length())
    return p


def dot(u, v):
    return (u & v).bit_count() & 1


def rref(rows, width):
    rows = list(rows)
    pivots = []
    r = 0
    for column in range(width):
        found = next((i for i in range(r, len(rows)) if (rows[i] >> column) & 1), None)
        if found is None:
            continue
        rows[r], rows[found] = rows[found], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> column) & 1):
                rows[i] ^= rows[r]
        pivots.append(column)
        r += 1
    return rows[:r], pivots


def nullspace(rows, width):
    rows, pivots = rref(rows, width)
    free = [i for i in range(width) if i not in pivots]
    basis = []
    for f in free:
        v = 1 << f
        for row, p in zip(rows, pivots):
            if (row >> f) & 1:
                v ^= 1 << p
        assert all(dot(row, v) == 0 for row in rows)
        basis.append(v)
    return basis


def span_words(basis):
    words = [0]
    for v in basis:
        words += [w ^ v for w in words]
    return words


def generator(p, L, j):
    assert 0 <= j <= L
    n = (p.bit_length() - 1) * L
    k = (p.bit_length() - 1) * (L - j)
    g = power(p, j)
    rows = [g << i for i in range(k)]
    assert all(row.bit_length() <= n for row in rows)
    return n, k, rows


def gram_rank(rows):
    gram = [sum(dot(u, v) << i for i, v in enumerate(rows)) for u in rows]
    return len(rref(gram, len(rows))[1])


def check_lcd(a, T, j):
    p = 1 ^ (1 << a) ^ (1 << (2 * a))
    L = 1 << T
    n, k, rows = generator(p, L, j)
    assert len(rref(rows, n)[1]) == k, ('generator rank', a, T, j)
    assert gram_rank(rows) == k, ('nontrivial hull', a, T, j)


def fidelity_example_6():
    # Source Example 6, Section V: <x^3+x+1> modulo (x^3+x+1)^8.
    p, L, j = 0b1011, 8, 1
    n, k, rows = generator(p, L, j)
    rank = len(rref(rows, n)[1])
    assert (n, k, rank, gram_rank(rows)) == (24, 21, 21, 21)
    dual = nullspace(rows, n)
    assert len(dual) == 3
    dual_words = span_words(dual)
    assert len(set(dual_words)) == 8
    assert all(all(dot(c, row) == 0 for row in rows) for c in dual_words)
    dual_distance = min(c.bit_count() for c in dual_words if c)
    assert dual_distance == 13
    # A singleton cannot be divisible by p; a weight-two word is divisible by p.
    assert all(remainder(1 << i, p) != 0 for i in range(n))
    weight_two = (1 << 7) ^ 1
    assert remainder(weight_two, p) == 0
    assert weight_two.bit_count() == 2 and weight_two.bit_length() <= n
    print('SOURCE_FIDELITY Example 6: primal [24,21,2], dual [24,3,13], Gram rank 21')


def main():
    fidelity_example_6()
    limits = {1: 7, 2: 4, 3: 4, 9: 2}
    nontrivial = 0
    cases = set()
    for a, bound in limits.items():
        for T in range(1, bound + 1):
            for j in range((1 << T) + 1):
                check_lcd(a, T, j)
                cases.add((a, T, j))
                nontrivial += int(0 < j < (1 << T))
    assert nontrivial == 303
    print('LCD_CHECK source-search ranges: 303 nontrivial codes and 34 endpoint codes; all full Gram rank')
    # Additional reducible/non-source spacings and the T=0 extension.
    for a in range(1, 7):
        for T in range(4):
            for j in range((1 << T) + 1):
                if (a, T, j) not in cases:
                    check_lcd(a, T, j)
                    cases.add((a, T, j))
    print('LCD_CHECK total distinct (a,T,j):', len(cases), '; a in {1,2,3,4,5,6,9}; includes all j and endpoints in each stated grid')
    # Independently exercise the proposed gap lemma for every small multiplier.
    gap_tests = 0
    for a, Tmax in [(1, 3), (2, 2), (3, 2)]:
        p = 1 ^ (1 << a) ^ (1 << (2 * a))
        for T in range(Tmax + 1):
            L = 1 << T
            assert power(p, L) == (1 ^ (1 << (a * L)) ^ (1 << (2 * a * L)))
            for j in range(1, L):
                g = power(p, j)
                bound = a * (L - j)
                mask = ((1 << (a * L - a * j)) - 1) << (a * j)
                for h in range(1 << bound):
                    assert (mul(g, h) & mask == 0) == (h == 0), ('gap', a, T, j, h)
                    gap_tests += 1
    print('GAP_CHECK all multipliers in grids (a=1,T<=3), (a=2 or 3,T<=2):', gap_tests)
    # Negative controls test that the LCD checker does not always answer yes.
    for p, L, j in [(0b111, 3, 1), (0b1011, 2, 1)]:
        n, k, rows = generator(p, L, j)
        assert gram_rank(rows) < k
    print('NEGATIVE_CONTROLS two non-LCD codes detected')
    print('SCOPE exact finite corroboration only; no universal proof or Lean-kernel verification')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
