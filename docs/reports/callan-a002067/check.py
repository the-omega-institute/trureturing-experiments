#!/usr/bin/env python3
"""Exact finite checks for David Callan's two OEIS A002067 conjectures.

Source: https://oeis.org/A002067, revision 116, COMMENTS dated July 21 2017.
Python standard library only. No floating point and no Lean proof.
"""
from collections import Counter
from fractions import Fraction
from functools import cache
from math import comb, factorial, prod
import json

SOURCE_TERMS = (
    1, 1, 7, 127, 4369, 243649, 20036983, 2280356863,
    343141433761, 65967241200001, 15773461423793767,
)
MAX_N = 10


def source_sequence(limit):
    a = [1]
    for n in range(1, limit + 1):
        a.append(sum(comb(2*n, 2*k)*a[k]*a[n-k-1] for k in range(n)))
    return a


@cache
def hafnian(labels):
    """Pair the least remaining label with every possible partner, once each."""
    if not labels:
        return 1
    i, *tail = labels
    return sum(i * hafnian(tuple(tail[:j] + tail[j+1:]))
               for j in range(len(tail)))


def dyck_paths(n):
    def visit(word, ups, downs):
        if ups == downs == n:
            yield word
        if ups < n:
            yield from visit(word + 'U', ups + 1, downs)
        if downs < ups:
            yield from visit(word + 'D', ups, downs + 1)
    yield from visit('', 0, 0)


def path_statistics(word):
    height = 0
    up_positions, upper_heights, down_start_heights = [], [], []
    for j, step in enumerate(word, 1):
        if step == 'U':
            height += 1
            up_positions.append(j)
            upper_heights.append(height)
            assert j + height == 2*len(up_positions)
        else:
            assert step == 'D' and height > 0
            down_start_heights.append(height)
            height -= 1
    assert height == 0
    assert Counter(upper_heights) == Counter(down_start_heights)
    weight = prod(up_positions)*prod(upper_heights)
    return up_positions, upper_heights, down_start_heights, weight


def matching_list(labels):
    """Explicit enumeration, separate from the memoized hafnian sum."""
    if not labels:
        yield ()
    else:
        i = labels[0]
        for k in range(1, len(labels)):
            j = labels[k]
            rest = labels[1:k] + labels[k+1:]
            for matching in matching_list(rest):
                yield ((i, j),) + matching


def matching_fibres(n):
    counts, weights = Counter(), Counter()
    total = 0
    for matching in matching_list(tuple(range(1, 2*n+1))):
        openings = {i for i, _ in matching}
        word = ''.join('U' if i in openings else 'D'
                       for i in range(1, 2*n+1))
        counts[word] += 1
        weights[word] += prod(i for i, _ in matching)
        total += 1
    expected_count = prod(range(1, 2*n, 2))
    assert total == expected_count
    assert set(counts) == set(dyck_paths(n))
    for word, count in counts.items():
        positions, heights, _, weight = path_statistics(word)
        assert count == prod(heights)
        assert weights[word] == weight
    return total


def polynomials(limit):
    """Coefficient arrays for P_(r+1) = P_r' + (r+1) X P_r."""
    result = [[1]]
    for r in range(limit):
        p = result[-1]
        q = [0]*(len(p)+1)
        for h, coefficient in enumerate(p):
            q[h+1] += (r+1)*coefficient
            if h:
                q[h-1] += h*coefficient
        result.append(q)
    return result


def enumerate_meanders(length):
    """Enumerate all nonnegative walks, using their literal step weights."""
    coefficients = [0]*(length+1)
    def visit(position, height, weight):
        if position == length:
            coefficients[height] += weight
            return
        visit(position+1, height+1, weight*(position+1))
        if height:
            visit(position+1, height-1, weight*height)
    visit(0, 0, 1)
    return coefficients


def multiply(f, g, degree):
    return [sum((f[k]*g[n-k] for k in range(n+1)), Fraction(0))
            for n in range(degree+1)]


def exponential(f, degree):
    assert f[0] == 0
    term = [Fraction(1)] + [Fraction(0)]*degree
    answer = term.copy()
    for k in range(1, degree+1):
        term = [c/k for c in multiply(term, f, degree)]
        answer = [x+y for x, y in zip(answer, term)]
    return answer


def formal_series_checks(a):
    degree = 2*MAX_N
    y = [Fraction(0)]*(degree+2)
    for n in range(MAX_N+1):
        y[2*n+1] = Fraction(a[n], factorial(2*n+1))
    quadratic = [v/2 for v in multiply(y, y, degree)]
    positive = exponential(quadratic, degree)
    negative = exponential([-v for v in quadratic], degree)
    A = [(n+1)*y[n+1] for n in range(degree+1)]
    assert A == [Fraction(a[n//2], factorial(n)) if n % 2 == 0 else Fraction(0)
                 for n in range(degree+1)]
    assert positive == A
    assert multiply(A, negative, degree) == [Fraction(1)] + [Fraction(0)]*degree
    C = [Fraction(1)-negative[0]] + [-v for v in negative[1:]]
    for m in range(degree-1):
        assert (m+1)*(m+2)*C[m+2] == A[m]
    assert C[0] == C[1] == 0
    return {'exp_y_squared_over_2_degree': degree,
            'A_times_exp_negative_degree': degree,
            'C_second_derivative_degree': degree-2}


def main():
    a = source_sequence(MAX_N)
    assert tuple(a) == SOURCE_TERMS
    example = path_statistics('UUDUUDDD')
    assert example[0] == [1, 2, 4, 5]
    assert example[1] == [1, 2, 2, 3]
    assert example[3] == 480
    print('FIDELITY ' + json.dumps({'path': 'UUDUUDDD', 'positions': example[0],
                                   'upper_heights': example[1], 'weight': example[3]}))
    P = polynomials(2*MAX_N)
    assert P[1] == [0, 1]
    assert P[2] == [1, 0, 2]
    assert P[3] == [0, 7, 0, 6]
    for r in range(13):
        assert enumerate_meanders(r) == P[r]
    for n in range(MAX_N+1):
        paths = list(dyck_paths(n))
        W = sum(path_statistics(p)[3] for p in paths)
        H = hafnian(tuple(range(1, 2*n+1)))
        assert len(paths) == comb(2*n, n)//(n+1)
        assert a[n] == H == W == P[2*n][0]
        print('IDENTITY ' + json.dumps({'n': n, 'recurrence': a[n], 'hafnian': H,
                                       'weighted_dyck': W, 'polynomial_constant': P[2*n][0],
                                       'dyck_count': len(paths)}))
    fibre_counts = [matching_fibres(n) for n in range(7)]
    print('MATCHING_FIBRES ' + json.dumps({'n_range': [0, 6],
                                         'explicit_matching_counts': fibre_counts}))
    print('MEANDER_POLYNOMIALS ' + json.dumps({'length_range': [0, 12]}))
    print('FORMAL_SERIES ' + json.dumps(formal_series_checks(a)))
    print('All arithmetic is exact. Finite checks do not prove the universal assertions.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
