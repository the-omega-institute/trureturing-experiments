#!/usr/bin/env python3
"""Exact finite checks for Mabilat's golden quiddity classification.

Pairs (a,b) mean a+b*phi, with phi**2=phi+1. All evidence uses integers.
The exhaustive irreducibility decision solves the endpoints of every possible
splice factor; it does not assume the proposed zero/triple reduction theorem.
"""
from collections import Counter
from itertools import product

ZERO, ONE, PHI = (0, 0), (1, 0), (0, 1)
I = (ONE, ZERO, ZERO, ONE)
MINUS_I = ((-1, 0), ZERO, ZERO, (-1, 0))


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    a, b = x
    c, d = y
    return a*c + b*d, a*d + b*c + b*d


def mmul(a, b):
    return (add(mul(a[0], b[0]), mul(a[1], b[2])),
            add(mul(a[0], b[1]), mul(a[1], b[3])),
            add(mul(a[2], b[0]), mul(a[3], b[2])),
            add(mul(a[2], b[1]), mul(a[3], b[3])))


def eta(x):
    return x, (-1, 0), ONE, ZERO


def monodromy(entries):
    m = I
    for x in entries:
        m = mmul(eta(x), m)
    return m


def golden_entries(c):
    return tuple((0, k) for k in c)


def is_quiddity(c):
    return len(c) >= 1 and monodromy(golden_entries(c)) in (I, MINUS_I)


def splice(a, b):
    assert len(a) >= 2 and len(b) >= 2
    return (a[0]+b[-1],) + a[1:-1] + (a[-1]+b[0],) + b[1:-1]


def dihedral(c):
    for d in (c, c[::-1]):
        for j in range(len(d)):
            yield d[j:] + d[:j]


def displayed(c):
    if len(c) == 5:
        return c in ((1,)*5, (-1,)*5)
    return len(c) == 4 and ((c[0] == c[2] == 0 and c[1] == -c[3]) or
                            (c[1] == c[3] == 0 and c[0] == -c[2]))


def decomposition(c):
    """Complete cut search for an input quiddity, with unbounded endpoints.

    Write B=(x, middle, y), P=M(middle)=[[p,q],[r,s]]. If M(B)=eI,
    then p=-e, x=-q/p, y=r/p. Conversely these suffice since det(P)=1.
    Thus all possible B are recovered by scanning the cyclic/reversed middle,
    testing p=+-1, and testing x,y have zero rational coordinate. No coefficient
    cutoff is imposed on either factor. A is then uniquely forced by splice.
    """
    for d in dihedral(c):
        for r in range(1, len(c)-2):
            p, q, t, _ = monodromy(golden_entries(d[-r:]))
            if p not in (ONE, (-1, 0)):
                continue
            x, y = neg(mul(q, p)), mul(t, p)  # p^{-1}=p
            if x[0] != 0 or y[0] != 0:
                continue
            b = (x[1],) + d[-r:] + (y[1],)
            a = (d[0]-y[1],) + d[1:-r-1] + (d[-r-1]-x[1],)
            assert len(a) >= 3 and len(b) >= 3
            assert splice(a, b) == d
            assert is_quiddity(b)
            assert is_quiddity(a)
            return a, b
    return None


def reduction(c):
    """Check the specific zero/pentagon splice claimed by the proof outline."""
    assert len(c) >= 6
    n = len(c)
    if 0 in c:
        j = c.index(0)
        # d=(a,0,b,c4,...,cn)
        d = c[j-1:] + c[:j-1] if j else c[-1:] + c[:-1]
        a, _, b = d[:3]
        factor_a = d[3:] + (a+b,)
        factor_b = (-b, 0, b, 0)
        target = d[3:] + d[:3]
    else:
        j = next(j for j in range(n)
                 if c[j] == c[(j+1) % n] == c[(j+2) % n] and abs(c[j]) == 1)
        start = (j+3) % n
        d = c[start:] + c[:start]
        p = c[j]
        factor_a = (d[0]-p,) + d[1:-4] + (d[-4]-p,)
        factor_b = (p,)*5
        target = d
    assert len(factor_a) >= 3 and len(factor_b) >= 3
    assert len(factor_a) < n and len(factor_b) < n
    assert splice(factor_a, factor_b) == target
    assert is_quiddity(factor_a) and is_quiddity(factor_b)


def regular_representation(x):
    a, b = x
    return ((a, b), (b, a+b))


def integer_block_matrix(m):
    return tuple(tuple(regular_representation(m[(i//2)*2+j//2])[i % 2][j % 2]
                       for j in range(4)) for i in range(4))


def int_matrix_mul(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4))
                       for j in range(4)) for i in range(4))


def independent_monodromy(c):
    # Independent 4x4 integer matrix product (regular representation of Q(phi)).
    m = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    for k in c:
        m = int_matrix_mul(integer_block_matrix(eta((0, k))), m)
    return m


def source_fidelity():
    examples = [((-1, 2, 4), (3, 0, 1), (0, 2, 7, 0)),
                ((-2, 1, 3, 1), (2, 3, 1), (-1, 1, 3, 3, 3)),
                ((2, 1, 0, 2), (1, -3, 2, 5, 1), (3, 1, 0, 3, -3, 2, 5))]
    for a, b, expected in examples:
        assert splice(a, b) == expected
    assert monodromy(golden_entries((1,)*5)) == MINUS_I
    assert monodromy(golden_entries((-1,)*5)) == I
    # Proposition 3.1, reproduced on all integer tuples in [-3,3] of sizes 1..4.
    for n in range(1, 5):
        for c in product(range(-3, 4), repeat=n):
            scalar = monodromy(tuple((k, 0) for k in c)) in (I, MINUS_I)
            expected = ((n == 2 and c == (0, 0)) or
                        (n == 3 and c in ((1, 1, 1), (-1, -1, -1))))
            if n == 4:
                a, b, x, y = c
                expected = ((a == -x and b == -y and a*b == 0) or
                            (a == x and b == y and a*b == 2))
            assert scalar == expected, ('source-fidelity', c)
    print('FIDELITY_OK: Definition 2.2 examples; signed golden pentagons; '
          'Proposition 3.1 on integer tuples [-3,3], lengths 1..4')


def identity_checks():
    u = eta(PHI)
    expected = [u, (PHI, neg(PHI), PHI, (-1, 0)),
                (ONE, neg(PHI), PHI, neg(PHI)),
                (ZERO, (-1, 0), ONE, neg(PHI)), MINUS_I]
    power = I
    for e in expected:
        power = mmul(power, u)
        assert power == e
    for a, b in product(range(-5, 6), repeat=2):
        lhs = mmul(mmul(eta((0, a)), eta(ZERO)), eta((0, b)))
        rhs = tuple(neg(x) for x in eta((0, a+b)))
        assert lhs == rhs
    for k in list(range(-20, 21)) + [10**80+7, -10**80-7]:
        for c in [(0, k, 0, -k), (k, 0, -k, 0)]:
            assert is_quiddity(c) and decomposition(c) is None
            assert independent_monodromy(c) == integer_block_matrix(monodromy(golden_entries(c)))
    print('IDENTITIES_OK: U^1..U^5; zero contraction on [-5,5]^2; '
          '43 parameters in each zero-quadruple family, including 81-digit values')


def exhaustive():
    max_n = 9
    alphabet = range(-2, 3)
    visited = Counter()
    solutions = Counter()
    irreducibles = Counter()
    cross_checks = 0

    def visit(c, m):
        nonlocal cross_checks
        n = len(c)
        if n:
            visited[n] += 1
            if m in (I, MINUS_I):
                solutions[n] += 1
                assert monodromy(golden_entries(c)) == m
                assert independent_monodromy(c) == integer_block_matrix(m)
                cross_checks += 1
                if n >= 3:
                    irreducible = decomposition(c) is None
                    assert irreducible == displayed(c), ('classification', c)
                    irreducibles[n] += irreducible
                    if 0 not in c:
                        assert any(c[j] == c[(j+1) % n] == c[(j+2) % n]
                                   and abs(c[j]) == 1 for j in range(n)), ('triple', c)
                    if n >= 6:
                        reduction(c)
        if n == max_n:
            return
        # Left multiplication by eta(k*phi); this is the source's product order.
        (a, b), (d, e), f, g = m
        for k in alphabet:
            nxt = ((k*b-f[0], k*(a+b)-f[1]),
                   (k*e-g[0], k*(d+e)-g[1]), (a, b), (d, e))
            visit(c+(k,), nxt)

    visit((), I)
    for n in range(1, max_n+1):
        assert visited[n] == 5**n
        print(f'EXHAUSTIVE n={n} tuples={visited[n]} quiddities={solutions[n]} '
              f'irreducibles={irreducibles[n]}')
    assert irreducibles == Counter({4: 9, 5: 2})
    print(f'EXHAUSTIVE_OK: {sum(visited.values())} tuples; '
          f'{cross_checks} scalar words independently checked by 4x4 integer matrices; '
          'all irreducibles match the displayed list; all long solutions admit '
          'the proposed strictly shorter splice')


if __name__ == '__main__':
    source_fidelity()
    identity_checks()
    exhaustive()
    print('ALL_CHECKS_OK')
