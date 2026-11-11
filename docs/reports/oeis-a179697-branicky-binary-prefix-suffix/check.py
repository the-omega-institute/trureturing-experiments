#!/usr/bin/env python3
"""Exact diagnostics for OEIS A179697; no claim of a universal computational proof."""
from hashlib import sha256
from pathlib import Path

LIMIT = 10**7


def reverse(n: int) -> int:
    assert n >= 0
    return int(format(n, 'b')[::-1], 2)


def reverse_shift(n: int) -> int:
    out = 0
    while n:
        out = 2*out + (n & 1)
        n >>= 1
    return out


def source_term(k: int) -> bool:
    r = reverse(k)
    d = abs(k-r)
    return k > 0 and d + reverse(d) == k and abs(k+r-reverse(k+r)) == k


def normal_form(k: int) -> bool:
    if k <= 0 or k % 2:
        return False
    p = k//2
    return p % 2 == 1 and reverse(p) == p and reverse(3*p) == 5*p


def target(k: int) -> bool:
    w = format(k, 'b')
    numeric = len(w) >= 5 and k//(2**(len(w)-4)) == 11 and k % 32 == 26
    literal = w.startswith('1011') and w.endswith('11010')
    assert numeric == literal
    return literal


def main() -> None:
    table_path = Path(__file__).with_name('b179697.txt')
    table_bytes = table_path.read_bytes()
    rows = [tuple(map(int, line.split())) for line in table_bytes.decode().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    assert [n for n, _ in rows] == list(range(1, 499))
    values = [k for _, k in rows]
    assert values == sorted(set(values))
    assert values[0] == 90
    assert reverse(100) == 19  # A030101's explicit example.
    assert [reverse(n) for n in range(16)] == [0,1,1,3,1,5,3,7,1,9,5,13,3,11,7,15]
    k = 90
    r, s = reverse(k), k+reverse(k)
    assert (r, s, reverse(s), abs(s-reverse(s)), abs(k-r)+reverse(abs(k-r))) == (45,135,225,90,90)
    assert source_term(k) and normal_form(k) and target(k)
    assert reverse(0) == 0 and abs(0+reverse(0)-reverse(0+reverse(0))) == 0
    assert abs(0-reverse(0))+reverse(abs(0-reverse(0))) == 0
    assert not source_term(0) and not target(0)
    print('Fidelity: A030101(100)=19; A030101 first 16 values match; A179697 example k=90 verified.')
    print('Boundary: k=0 satisfies both raw equations and is explicitly excluded from positive membership.')
    for n in range(1, 65536):
        r = reverse(n)
        odd_part = n//(n & -n)
        assert reverse_shift(n) == r
        assert r > 0 and r % 2 == 1 and r < 2*n
        assert reverse(r) == odd_part
        if n % 2 == 0:
            assert r < n
        else:
            assert r.bit_length() == n.bit_length()
    print('Reversal diagnostics: all 1<=n<=65535; shift implementation, odd-part involution, parity, and bounds pass.')
    for k in values:
        assert reverse_shift(k) == reverse(k)
        assert source_term(k) and normal_form(k) and target(k), k
    print(f'Published table: all {len(values)} entries satisfy both equations, normal form, and target; maximum={max(values)}.')
    print('Table sha256: '+sha256(table_bytes).hexdigest())
    found = []
    for k in range(1, LIMIT+1):
        actual = source_term(k)
        normal = normal_form(k)
        assert actual == normal, (k, actual, normal)
        if actual:
            assert target(k), k
            found.append(k)
    assert found == [k for k in values if k <= LIMIT]
    print(f'Exhaustive exact scan: 1<=k<={LIMIT}; {len(found)} solutions; matches the entire published table prefix in this interval.')
    print('Solutions: '+','.join(map(str, found)))
    print('Normal-form equivalence and prefix/suffix hold for every scanned k; this is finite evidence only.')
    print('Small odd palindromes (p,R(3p),5p): '+repr([(p,reverse(3*p),5*p) for p in [1,3,5,7]]))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
