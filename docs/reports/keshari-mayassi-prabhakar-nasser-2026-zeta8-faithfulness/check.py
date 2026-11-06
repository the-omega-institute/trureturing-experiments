#!/usr/bin/env python3
"""Exact checks for arXiv:2605.13090v1, Question 4.3; requires SymPy.
Rational functions over Q(a,b,y0,y1,y2), with b and y parameters nonzero.
Normal forms are for free products of involutions, not the source group.
"""
from collections import Counter
from itertools import product
import sympy as sp


def s(i):
    return ('s', i, -1)


def r(i, alpha):
    return ('r', i, alpha)


def relations(k, n):
    """All source relations (11)-(19), with zero-based generator indices."""
    for i in range(n - 1):
        yield 11, (s(i), s(i)), ()
    for i, j in product(range(n - 1), repeat=2):
        if abs(i - j) >= 2:
            yield 12, (s(i), s(j)), (s(j), s(i))
    for i, alpha in product(range(n - 1), range(k)):
        yield 13, (r(i, alpha), r(i, alpha)), ()
    for i, j, alpha, beta in product(range(n - 1), range(n - 1), range(k), range(k)):
        if abs(i - j) >= 2:
            yield 14, (r(i, alpha), r(j, beta)), (r(j, beta), r(i, alpha))
    for i, j, alpha in product(range(n - 1), range(n - 1), range(k)):
        if abs(i - j) >= 2:
            yield 15, (r(i, alpha), s(j)), (s(j), r(i, alpha))
    for i, alpha in product(range(n - 2), range(k)):
        yield 16, (r(i, alpha), r(i + 1, alpha), r(i, alpha)), (r(i + 1, alpha), r(i, alpha), r(i + 1, alpha))
    for i, alpha, beta in product(range(n - 2), range(k), range(k)):
        if alpha < beta:
            yield 17, (r(i, alpha), r(i + 1, beta), r(i, beta)), (r(i + 1, beta), r(i, beta), r(i + 1, alpha))
            yield 18, (r(i, alpha), r(i + 1, alpha), r(i, beta)), (r(i + 1, beta), r(i, alpha), r(i + 1, alpha))
    for i, alpha in product(range(n - 2), range(k)):
        yield 19, (r(i, alpha), r(i + 1, alpha), s(i)), (s(i + 1), r(i, alpha), r(i + 1, alpha))


def reduce_involutions(word):
    out = []
    for letter in word:
        if out and out[-1] == letter:
            out.pop()
        else:
            out.append(letter)
    return tuple(out)


def quotient_word(word):
    return reduce_involutions(g[2] for g in word if g[0] == 'r')


def equal_matrix(left, right):
    return all(sp.cancel(e) == 0 for e in left - right)


def local(block, i, n):
    out = sp.eye(n)
    out[i:i + 2, i:i + 2] = block
    return out


a, b, y0, y1, y2 = sp.symbols('a b y0 y1 y2')
ys = (y0, y1, y2)
S = sp.Matrix([[-a, -(a*a - 1)/b], [b, a]])


def R(y):
    return sp.Matrix([[0, 1/y], [y, 0]])


def image(word, n, avalue=a, yvalues=ys):
    out = sp.eye(n)
    for kind, i, alpha in word:
        block = S.subs(a, avalue) if kind == 's' else R(yvalues[alpha])
        out = out * local(block, i, n)
    return out.applyfunc(sp.cancel)


def check_source_matrices(k, n):
    counts = Counter()
    for number, left, right in relations(k, n):
        assert equal_matrix(image(left, n), image(right, n)), (number, left, right)
        counts[number] += 1
    print(f'SOURCE_RELATIONS k={k} n={n}: ' + ', '.join(f'{number}:{counts[number]}' for number in range(11, 20)))


def check_quotient(k):
    counts = Counter()
    for number, left, right in relations(k, 3):
        lred, rred = quotient_word(left), quotient_word(right)
        assert lred == rred, (number, left, right, lred, rred)
        print(f'QUOTIENT k={k} relation={number} lhs={lred} rhs={rred}')
        counts[number] += 1
    print(f'QUOTIENT_ALL_RELATIONS k={k} n=3: ' + ', '.join(f'{number}:{counts[number]}' for number in range(11, 20)))


def main():
    assert not any(e.has(sp.Float) for e in list(S) + list(R(y0)))
    print('Exact arithmetic: rational functions over Q(a,b,y0,y1,y2); b,y0,y1,y2 nonzero.')
    print('Source: arXiv:2605.13090v1, relations (11)-(19), Theorem 4.2(b), Question 4.3.')
    check_source_matrices(2, 3)
    check_source_matrices(3, 3)
    check_source_matrices(3, 4)  # exercises far commutations
    fidelity = [(-1, (s(0), r(1, 0))*4), (0, (s(0), r(1, 0))*3), (1, (s(1), r(0, 0))*4)]
    for avalue, word in fidelity:
        assert equal_matrix(image(word, 3, avalue=avalue), sp.eye(3))
        print(f'FIDELITY Theorem 4.2(b) a={avalue}: product is I_3 (symbolic b,y0).')
    check_quotient(2)
    old_word = (r(0, 0), r(0, 1))*2
    assert quotient_word(old_word) == (0, 1, 0, 1)
    assert equal_matrix(image(old_word, 3, yvalues=(y0, -y0, y2)), sp.eye(3))
    print('HANDOFF k=2 n=3 y1=-y0: w=(rho_1^0 rho_1^1)^2 maps to I_3; quotient normal form=(0,1,0,1).')
    print('HANDOFF unrestricted y: image(w)=' + str(image(old_word, 3)))
    check_quotient(3)
    first = (r(0, 0), r(0, 1))
    second = (r(0, 0), r(0, 2))
    word = first + second + first[::-1] + second[::-1]
    normal_form = quotient_word(word)
    assert normal_form == (0, 1, 0, 2, 1, 0, 2, 0)
    assert normal_form and all(normal_form[i] != normal_form[i+1] for i in range(len(normal_form)-1))
    assert equal_matrix(image(first, 3), sp.diag(y1/y0, y0/y1, 1))
    assert equal_matrix(image(second, 3), sp.diag(y2/y0, y0/y2, 1))
    assert equal_matrix(image(word, 3), sp.eye(3))
    print('GENERIC k=3 n=3: P=diag(y1/y0,y0/y1,1), Q=diag(y2/y0,y0/y2,1).')
    print('GENERIC quotient reduced word=' + str(normal_form) + '; length=8; nonempty.')
    print('GENERIC commutator image=' + str(image(word, 3)))
    # A separate finite action witnesses nonidentity without relying on the
    # normal-form implementation: u=(01), v=(12), t=(23) in S4.
    involutions = ((1, 0, 2, 3), (0, 2, 1, 3), (0, 1, 3, 2))
    def permutation_image(source_word):
        perm = tuple(range(4))
        for kind, _, alpha in source_word:
            gen = tuple(range(4)) if kind == 's' else involutions[alpha]
            perm = tuple(perm[gen[j]] for j in range(4))
        return perm
    for number, left, right in relations(3, 3):
        assert permutation_image(left) == permutation_image(right), number
    action = permutation_image(word)
    assert action != tuple(range(4))
    print('INDEPENDENT S4 action: u=(01), v=(12), t=(23); all 20 source relations hold.')
    print('INDEPENDENT S4 image of kernel word=' + str(action) + '; nonidentity.')
    print('GENERIC conclusion: one fixed nontrivial kernel word for every a, b!=0, y0*y1*y2!=0; no parameter collision or root-of-unity assumption.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
