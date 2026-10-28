#!/usr/bin/env python3
"""Exact finite checks for Freeman's printed backward-ratio limsup conjecture.

This program checks finite instances and construction identities. It does not
certify the infinite limit or eventual nonperiodicity. See README.md for scope.
"""
from fractions import Fraction
from itertools import product

IDENTITY = (1, 0, 0, 1)


def multiply(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return (a * e + b * g, a * f + b * h,
            c * e + d * g, c * f + d * h)


def matrix(word):
    result = IDENTITY
    for digit in word:
        assert isinstance(digit, int) and digit > 0
        result = multiply(result, (digit, 1, 1, 0))
    return result


def cf(word):
    """Evaluate from the right, independently of transfer matrices."""
    assert word and all(isinstance(a, int) and a > 0 for a in word)
    result = Fraction(word[-1])
    for digit in reversed(word[:-1]):
        result = digit + 1 / result
    return result


def is_m_palindrome(word, m):
    return bool(word) and cf(word) == m * cf(word[::-1])


def transpose(x):
    a, b, c, d = x
    return a, c, b, d


def theta(x, m):
    a, b, c, d = x
    return a, m * c, Fraction(b, m), d


def encode(word, m, r, s):
    return tuple(d for letter in word
                 for t in [r if letter == 'a' else s] for d in (m * t, t))


def all_prefix_lengths(digits, m):
    """Scan EVERY positive prefix length, not just designated block ends."""
    result = []
    p, q, u, v = IDENTITY
    for length, digit in enumerate(digits, 1):
        p, q, u, v = p * digit + q, p, u * digit + v, u
        if q == m * u:
            result.append(length)
    return result


def stages(count):
    result = ['a']
    for k in range(count):
        w = result[-1]
        result.append(w * (k + 2) + 'b' * (k + 1) + w * (k + 2))
    return result


def fidelity():
    # Section 2's explicit example of distinct words with the same value.
    assert cf((1, 1, 1)) == cf((1, 2)) == Fraction(3, 2)
    # Source labels X:short_pal and L:periodic_m_pals: two-digit words and their squares.
    for m, n in product(range(1, 7), repeat=2):
        word = (m * n, n)
        assert is_m_palindrome(word, m)
        assert is_m_palindrome(word * 2, m)
    # Section 7: the first Fibonacci words and all S_n^* used here.
    fib = [0, 1]
    for _ in range(12):
        fib.append(fib[-1] + fib[-2])
    words = ['a', 'ab']
    for _ in range(8):
        words.append(words[-1] + words[-2])
    parameters = [(m, r, s) for m in range(1, 5)
                  for r in range(1, 4) for s in range(1, 4) if r != s]
    count = 0
    for n in range(2, 10):
        central = words[n][:-2]
        assert central == central[::-1]
        for m, r, s in parameters:
            digits = encode(central, m, r, s)
            assert len(digits) == 2 * fib[n + 2] - 4
            assert is_m_palindrome(digits, m)
            count += 1
    print('FIDELITY: [1,1,1]=[1,2]=3/2; 36 two-digit examples and squares;'
          f' {count} Fibonacci-derived source prefixes (n=2..9)')


def independent_arithmetic():
    count = 0
    for length in range(1, 8):
        for word in product((1, 2, 3), repeat=length):
            p, q, u, v = x = matrix(word)
            assert p > 0 and q > 0 and u > 0
            assert cf(word) == Fraction(p, u)
            assert cf(word[::-1]) == Fraction(p, q)
            assert matrix(word[::-1]) == transpose(x)
            assert p * v - q * u == (-1) ** length
            for m in range(1, 6):
                assert is_m_palindrome(word, m) == (q == m * u)
            count += 1
    # Verify the reversal of multiplication in the weighted transpose separately.
    samples = [(a, b, c, d) for a, b, c, d in product((-1, 0, 2), repeat=4)]
    for m in range(2, 6):
        for x, y in zip(samples, reversed(samples)):
            assert theta(multiply(x, y), m) == multiply(theta(y, m), theta(x, m))
    print(f'ARITHMETIC: {count} positive words of length 1..7; rational CF,'
          ' transpose, determinant, and m=1..5 criterion agree')
    print('WEIGHTED TRANSPOSE: 324 matrix-pair cases, m=2..5')


def counterfamily():
    ws = stages(4)
    assert [len(w) for w in ws] == [1, 5, 32, 259, 2594]
    for k, (w, nxt) in enumerate(zip(ws, ws[1:])):
        assert w and w[0] == w[-1] == 'a' and w == w[::-1]
        assert nxt.startswith(w)
        start = (k + 2) * len(w)
        assert nxt[start:start + k + 1] == 'b' * (k + 1)
        assert nxt[start + k + 1] == 'a'
    parameters = [(m, r, s) for m in range(2, 6)
                  for r in range(1, 4) for s in range(1, 4) if r != s]
    scans, pairs = 0, 0
    for m, r, s in parameters:
        digits = encode(ws[-1], m, r, s)
        assert set(digits) <= {m * r, r, m * s, s}
        assert digits[1::2] == tuple(r if x == 'a' else s for x in ws[-1])
        lengths = all_prefix_lengths(digits, m)
        assert lengths and all(a < b for a, b in zip(lengths, lengths[1:]))
        scans += len(digits)
        for k, w in enumerate(ws[:-1]):
            t = 2 * (k + 1) * len(w)
            tprime = 2 * (k + 2) * len(w)
            assert is_m_palindrome(digits[:t], m)
            assert is_m_palindrome(digits[:tprime], m)
            idx = lengths.index(t)
            successor = lengths[idx + 1]
            actual = Fraction(t, successor)
            lower = Fraction(k + 1, k + 2)
            assert successor <= tprime
            assert lower <= actual < 1
            # For k>=1 the EXACT rational lower bound exceeds 1/phi.
            # On x>=0, x>1/phi iff x*x+x>1 (phi=(1+sqrt(5))/2).
            if k >= 1:
                assert lower * lower + lower > 1
            pairs += 1
        if (m, r, s) == (2, 2, 3):
            print('CONCRETE: m=2 r=2 s=3; alphabet=2,3,4,6;'
                  f' scanned_digits={len(digits)} qualifying_prefixes={len(lengths)}')
            for k, w in enumerate(ws[:-1]):
                t = 2 * (k + 1) * len(w)
                nxt = lengths[lengths.index(t) + 1]
                print(f'  k={k} W_length={len(w)} t={t} actual_next={nxt}'
                      f' ratio={Fraction(t, nxt)} lower={Fraction(k + 1, k + 2)}')
        # Finite falsification of all tail/period pairs in this stated window.
        for start in range(33):
            for period in range(1, 33):
                assert any(digits[i] != digits[i + period]
                           for i in range(start, len(digits) - period))
    # Length recurrence and prescribed ratios can be checked without huge words.
    length = 1
    for k in range(100):
        t, tprime = 2 * (k + 1) * length, 2 * (k + 2) * length
        assert Fraction(t, tprime) == Fraction(k + 1, k + 2)
        assert t >= 2 * (k + 1)
        length = 2 * (k + 2) * length + k + 1
    print(f'COUNTERFAMILY: {len(parameters)} parameter triples; {scans} prefix positions;'
          f' {pairs} designated pairs checked against the COMPLETE finite prefix lists')
    print('APERIODICITY FINITE SCREEN: all tail starts 0..32 and periods 1..32 rejected'
          ' for each parameter triple; not an infinite aperiodicity certificate')
    print('LENGTH/RATIO: exact recurrence and (k+1)/(k+2), k=0..99')
    print('SCOPE: finite evidence only; the infinite construction, aperiodicity and'
          ' limsup=1 require the mathematical argument in the preregistration')


def main():
    fidelity()
    independent_arithmetic()
    counterfamily()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
