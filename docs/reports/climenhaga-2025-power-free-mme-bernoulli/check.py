#!/usr/bin/env python3
"""Exact finite checks for the aligned-product Bernoulli research proposal.

No floating point, entropy approximation, or claimed test of measure isomorphism.
Run with Python 3; sympy is used only for an exact polynomial expansion.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from sympy import Poly, symbols, expand


def words(d, n):
    return product(range(d), repeat=n)


def fourth(w):
    return bool(w) and len(w) % 4 == 0 and w == w[:len(w)//4] * 4


@lru_cache(None)
def good(w):
    return all(not fourth(w[:4*k]) and not fourth(w[-4*k:])
               for k in range(1, len(w)//4 + 1))


@lru_cache(None)
def locally_admissible(w):
    # Avoid exponents strictly greater than 12. This suffices for either
    # convention at every beta > 12. A prohibited factor has a prefix
    # of length 12*p+1 with period p, and conversely.
    n = len(w)
    for p in range(1, (n-1)//12 + 1):
        m = 12*p + 1
        for start in range(n-m+1):
            v = w[start:start+m]
            if all(v[i] == v[i % p] for i in range(m)):
                return False
    return True


def b_block(w):
    if not w or len(w) % 4:
        return False
    return any(fourth(tuple(a[c] for a in w)) for c in (0, 1))


def bstar(w):
    reached = {0}
    for j in range(4, len(w)+1, 4):
        if any(i in reached and b_block(w[i:j]) for i in range(0, j, 4)):
            reached.add(j)
    return len(w) in reached


def synchronized_cut(w):
    n = len(w)
    left = max(i for i in range(n+1) if bstar(w[:i]))
    right = max(j for j in range(n-left+1) if bstar(w[n-j:]))
    middle = w[left:n-right]
    assert left + right <= n
    assert w == w[:left] + middle + w[n-right:]
    assert all(good(tuple(a[c] for a in middle)) for c in (0, 1))
    return left, right


def main():
    # Source fidelity: Climenhaga v2 (proof of entropy gap), G4 and G8.
    for d in range(2, 7):
        n4 = sum(good(w) for w in words(d, 4))
        assert n4 == d**4-d
        print(f'source G4: d={d}, count={n4}=d^4-d')
    n8 = sum(good(w) and locally_admissible(w) for w in words(2, 8))
    assert n8 >= 2**8-2*2**5-2**2 == 188
    print(f'source binary G8: exact={n8}, source lower bound=188')

    # New helper: count good locally admissible binary length-16 words.
    n16 = sum(good(w) and locally_admissible(w) for w in words(2, 16))
    lower = 2**16-2*sum(2**(16-3*k) for k in range(1, 5))-4*2**4
    assert lower == 46752 and lower > 2**15 and n16 == 49100 and n16 >= lower
    print(f'binary G16: exact={n16}, union lower={lower}>2^15={2**15}')

    # Universal d>=3 positivity reduced to nonnegative coefficients at d=3+x.
    x = symbols('x')
    poly = Poly(expand((x+3)**8-2*(x+3)**5-(x+3)**2-(x+3)**7), x)
    assert all(c > 0 for c in poly.all_coeffs())
    assert 3**9 < 3**11
    assert (3*2**5)**3 < 2**20
    print(f'ternary-and-larger count surplus coefficients={poly.all_coeffs()}')
    print(f'product entropy gap, binary: 96^3={(3*2**5)**3}<2^20={2**20}')
    print(f'product entropy gap, d>=3: 3^9={3**9}<3^11={3**11}<=d^11')

    # Small tests of the LOCAL equal-length gluing statement used before
    # compactness/extendability. They do not test infinite extendability.
    glue_cases = 0
    for d, max_n, tau in ((2, 6, 2), (3, 4, 1)):
        bridges = list(words(d, tau))
        for n in range(1, max_n+1):
            g = [w for w in words(d, n) if good(w) and locally_admissible(w)]
            for v in g:
                for w in g:
                    assert any(good(v+u+w) and locally_admissible(v+u+w) for u in bridges)
                    glue_cases += 1
    print(f'local equal-length gluing: {glue_cases} ordered pairs checked')

    # All paired binary words through length 8, including empty words and
    # overlapping potential prefix/suffix strips. At these lengths every
    # coordinate is locally admissible. The cut uses one SHARED interval.
    paired_alphabet = tuple(product(range(2), repeat=2))
    cut_cases = 0
    for n in range(9):
        count_bstar = 0
        for w in product(paired_alphabet, repeat=n):
            synchronized_cut(w)
            count_bstar += bstar(w)
            cut_cases += 1
        if n and n % 4 == 0:
            k = n//4
            bound = 2*3**(k-1)*2**(5*k)
            assert count_bstar <= bound
            print(f'paired B-star length={n}: exact={count_bstar}, bound={bound}')
    print(f'synchronized product cuts: {cut_cases} paired words checked, lengths 0..8')

    # All composition lengths up to 12, exact integer recurrence.
    for d in range(2, 7):
        a = [1]
        for k in range(1, 13):
            a.append(sum(2*d**(5*j)*a[k-j] for j in range(1, k+1)))
            closed = 2*3**(k-1)*d**(5*k)
            composition_sum = d**(5*k)*sum(2**m*comb(k-1, m-1) for m in range(1, k+1))
            assert a[k] == closed == composition_sum
    print('composition recurrence: d=2..6, k=1..12, all exact')

    # Cyclic factor calculation: diagonal union invariant, mass exactly 1/p.
    for p in range(1, 65):
        diagonal = {(i, i) for i in range(p)}
        image = {((i+1) % p, (j+1) % p) for i, j in diagonal}
        mass = sum((Fraction(1, p*p) for _ in diagonal), Fraction(0))
        assert image == diagonal and mass == Fraction(1, p)
        assert (mass in (0, 1)) == (p == 1)
    print('cyclic-period diagonal invariant: p=1..64, mass=1/p exactly')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
