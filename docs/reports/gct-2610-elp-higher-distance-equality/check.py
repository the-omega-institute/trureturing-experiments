#!/usr/bin/env python3
"""Exact finite checks for Ge--Cai--Tang arXiv:2610.15252v1, VI-C.

All evidence is integer arithmetic or rational coefficients of logarithms.
NumPy is used only for bounded int64 arrays. No floating-point operations occur.
The checks do not prove the unbounded parameter statement.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
import numpy as np


def subsets(n):
    for size in range(n + 1):
        yield from combinations(range(n), size)


def words(b, k):
    """Enumerate all b**k words, once, as bounded exact integers."""
    out = np.empty((b**k, k), dtype=np.int64)
    labels = np.arange(b**k, dtype=np.int64)
    for j in range(k):
        out[:, j] = labels % b
        labels //= b
    assert not labels.any()
    return out


def prime_factors(n):
    out = Counter()
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] += 1
            n //= p
        p += 1
    if n > 1:
        out[n] += 1
    return out


def entropy_signature(counts):
    """H(full word | observed coordinates) = sum_p coeff[p] log_q(p).

    Each observed fibre is uniform, with probability size/total. The projected
    observed variable itself is NOT assumed uniform on its image.
    """
    total = sum(counts)
    signature = Counter()
    for size in counts:
        for p, exponent in prime_factors(size).items():
            signature[p] += Fraction(size * exponent, total)
    return dict(signature)


def fidelity():
    # Source equation (50), Construction C, Proposition 3 and Corollary 5.
    table = ((0, 1, 2, 3, 4), (1, 0, 3, 4, 2), (2, 3, 4, 0, 1),
             (3, 4, 1, 2, 0), (4, 2, 0, 1, 3))
    assert all(set(row) == set(range(5)) for row in table)
    assert all({table[i][j] for i in range(5)} == set(range(5)) for j in range(5))

    def fold(xs):
        value = xs[0]
        for x in xs[1:]:
            value = table[value][x]
        return value

    m, t = 3, 2
    initial = [x + (fold(x),) for x in product(range(3), repeat=m)]
    final = [x + (fold((fold(x[:m]), fold(x[m:]))),)
             for x in product(range(3), repeat=t*m)]
    assert len(set(initial)) == 27 and len(set(final)) == 729
    for code in (initial, final):
        n = len(code[0])
        attained = set()
        for J in subsets(n):
            comp = tuple(i for i in range(n) if i not in J)
            counts = Counter(tuple(x[i] for i in comp) for x in code)
            bound = 3 ** max(0, len(J)-1)
            assert max(counts.values()) <= bound
            if all(v == bound for v in counts.values()):
                assert entropy_signature(list(counts.values())) == (
                    {3: Fraction(len(J)-1)} if len(J) > 1 else {})
                attained.add(len(J))
        assert attained == set(range(n+1))
        d = min(sum(a != b for a, b in zip(x, y))
                for i, x in enumerate(code) for y in code[i+1:])
        assert d == 2
    # h_3 = 2 log_5 3, delta_(3 log_5 3) = 4, U=3, R=1.
    assert [max(0, j-1) for j in range(8)] == [0, 0, 1, 2, 3, 4, 5, 6]
    assert min(j for j in range(8) if max(0, j-1) >= 3) == 4
    assert 5**t < 3**(t*m) and 1 < 3 < 5
    return {'source': 'Construction C, Proposition 3, Corollary 5',
            'initial': [4, 27, 2], 'final': [7, 729, 2],
            'coordinate_subsets_checked': 16+128, 'h_m_log_5_3_coefficient': 2,
            'delta': 4, 'unchanged_per_input': 3, 'read_per_input': 1}


def determinant_mod(matrix, q):
    a = [list(row) for row in matrix]
    det = 1
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j] % q), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            det = -det
        v = a[j][j] % q
        det = det * v % q
        inv = pow(v, -1, q)
        for i in range(j+1, len(a)):
            factor = a[i][j] * inv % q
            for k in range(j, len(a)):
                a[i][k] = (a[i][k] - factor*a[j][k]) % q
    return det % q


def cauchy_case(r, t, m, q):
    assert r >= 2 and t >= 2 and m >= r+1 and q > t*m+r
    assert all(q % d for d in range(2, q) if d*d <= q)
    b, K = q-1, t*m
    A = np.array([[pow(r+j-l, -1, q) for j in range(K)]
                  for l in range(r)], dtype=np.int64)
    minors = 0
    for size in range(1, r+1):
        for rows in combinations(range(r), size):
            for cols in combinations(range(K), size):
                assert determinant_mod([[int(A[i,j]) for j in cols] for i in rows], q)
                minors += 1
    fibre_cases = 0
    assignments = 0
    lengths = []
    # Check each initial block and the final code, every selected coordinate set.
    for cols in [tuple(range(i*m, (i+1)*m)) for i in range(t)] + [tuple(range(K))]:
        B = A[:, cols]
        k = len(cols)
        attained = set()
        for unknown in subsets(k):
            X = words(b, len(unknown))
            parities = X @ B[:, unknown].T % q
            assert len(X) == b**len(unknown)
            for observed_parities in subsets(r):
                # Outside systematic values only translate these parity values;
                # translation permutes fibres and leaves all sizes unchanged.
                keys = np.zeros(len(X), dtype=np.int64)
                for l in observed_parities:
                    keys = keys*q + parities[:, l]
                counts = np.bincount(keys, minlength=q**len(observed_parities))
                counts = counts[counts != 0]
                assert int(counts.sum()) == b**len(unknown)
                j = len(unknown) + r-len(observed_parities)
                exponent = max(0, j-r)
                bound = b**exponent
                assert int(counts.max()) <= bound, (q, cols, unknown, observed_parities)
                # A singleton fibre at j=r checks the distance lower bound for
                # all codeword pairs, not just distances from zero.
                if j <= r:
                    assert (counts == 1).all()
                if not observed_parities:
                    assert len(counts) == 1 and int(counts[0]) == bound
                    attained.add(j)
                if j == 0:
                    attained.add(0)
                # Exact entropy is a weighted sum of log_q(fibre size).
                # The size bound proves its upper bound; these constant-size
                # witnesses attain it with coefficient exponent*log_q(b).
                fibre_cases += 1
                assignments += len(X)
        assert attained >= set(range(r, k+r+1)) | {0}
        assert all(A[l, cols[0]] != 0 for l in range(r))
        witness = [1] + [0]*(k-1) + [int(B[l,0]) for l in range(r)]
        assert sum(v != 0 for v in witness) == r+1
        lengths.append([k+r, b**k, r+1])
    # Conversion checked on every possible final systematic word (= input tuple).
    X = words(b, K)
    total = np.zeros((len(X), r), dtype=np.int64)
    for i in range(t):
        total += X[:, i*m:(i+1)*m] @ A[:, i*m:(i+1)*m].T % q
    assert np.array_equal(total % q, X @ A.T % q)
    assert len(X) == b**K  # systematic coordinates make the conversion bijective
    # All logarithmic comparisons are reduced to exact powers, q>1.
    assert q**(m-1) < b**m < q**m  # ma>m-1>=r
    assert q**(r-1) < b**r < q**r  # ceil(r log_q b)=r
    assert q**(t*r) < b**K         # tr<k_F
    profile_coefficients = [max(0, j-r) for j in range(K+r+1)]
    delta = min(j for j, c in enumerate(profile_coefficients) if c >= (t-1)*m)
    assert delta == (t-1)*m+r
    assert min(m+r, K+r-delta) == m
    assert m-profile_coefficients[m] == r
    assert profile_coefficients[m] == m-r > 0
    # log_q(q-1) is irrational: q and q-1 are coprime and q is prime.
    # This exact check excludes q-power cardinality; irrationality in the
    # unbounded family is part of the paper proof, not a numerical estimate.
    for k in [m, K, m-r]:
        assert pow(b, k, q) != 0
    return {'parameters': {'r': r, 't': t, 'm': m, 'q': q},
            'matrix': A.tolist(), 'codes': lengths, 'nonempty_minors': minors,
            'coordinate_subsets_checked': fibre_cases,
            'systematic_assignments_enumerated_with_multiplicity': assignments,
            'conversion_input_tuples_checked': len(X),
            'final_elp_log_q_qminus1_coefficients': profile_coefficients,
            'delta': delta, 'theta_log_q_qminus1_coefficient': m-r,
            'unchanged_per_input': m, 'read_per_input': r,
            'nontriviality_integer_comparison': [q**(t*r), b**K]}


def main():
    print(json.dumps({'fidelity': fidelity()}, sort_keys=True))
    for case in [(2, 2, 3, 11), (2, 2, 3, 13)]:
        print(json.dumps({'expected_settlement_check': cauchy_case(*case)}, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
