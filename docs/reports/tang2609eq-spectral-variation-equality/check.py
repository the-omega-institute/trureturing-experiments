#!/usr/bin/env python3
"""Exact tests of the proposed answer to Tang, arXiv:2609.09177v2, Question 5.

The computations are finite tests, not a proof of the all-dimensional assertion.
SymPy uses rational characteristic polynomials and certified rational real-root
isolating intervals. No floating-point value is used as mathematical evidence.
"""
from collections import Counter
from itertools import permutations, product
from random import Random
import sys
import sympy as s

SEED = 260909177
R = s.Rational
I = s.I


def abs2(z):
    z = s.expand(z * s.conjugate(z))
    assert z.is_Rational, z
    return z


def d2(alpha, beta):
    return min(max(abs2(alpha[i] - beta[j]) for i, j in enumerate(p))
               for p in permutations(range(len(alpha))))


def separated(a, b, common):
    """The graph has separate labels for a,b and every common eigenvalue."""
    values = [a, b, *common]
    radius2 = abs2(a - b)
    reached, todo = {0}, [0]
    while todo:
        u = todo.pop()
        for v in range(len(values)):
            if v != u and v not in reached and abs2(values[u] - values[v]) < radius2:
                reached.add(v)
                todo.append(v)
    return 1 not in reached


def diagonal_rhs(a, b):
    n = len(a)
    if n == 1:
        return True
    if n == 2:
        return max(abs2(a[0]-b[0]), abs2(a[1]-b[1])) <= max(
            abs2(a[0]-b[1]), abs2(a[1]-b[0]))
    changed = [j for j in range(n) if a[j] != b[j]]
    if not changed:
        return True
    if len(changed) != 1:
        return False
    j = changed[0]
    return separated(a[j], b[j], [a[k] for k in range(n) if k != j])


def exhaustive_diagonal():
    # Enumerate every pair of diagonal matrices over {0,1,2} through n=4.
    for n in range(1, 5):
        counts = Counter()
        perms = list(permutations(range(n)))
        for a in product(range(3), repeat=n):
            for b in product(range(3), repeat=n):
                distances = [[(x-y)**2 for y in b] for x in a]
                match2 = min(max(distances[i][p[i]] for i in range(n)) for p in perms)
                total = sum(sorted(((a[i]-b[i])**2 for i in range(n)), reverse=True)[:(n+1)//2])
                expected = diagonal_rhs(a, b)
                assert match2 <= total, (a, b, match2, total)
                assert (match2 == total) == expected, (a, b, match2, total)
                assert (match2 < total) == (not expected)
                counts['equality' if expected else 'strict'] += 1
        print(f'EXHAUSTIVE_DIAGONAL n={n} cases={sum(counts.values())} '
              f'equality={counts["equality"]} strict={counts["strict"]}')


def zero(M):
    return all(s.cancel(x) == 0 for x in M)


def eigenlist(A):
    values = []
    for z, count in A.eigenvals().items():
        assert s.expand(z).is_algebraic
        values.extend([s.expand(z)] * count)
    assert len(values) == A.rows
    return values


def joint_pairs(A, B):
    """Recover simultaneous diagonal entries from eigenspace restrictions."""
    assert zero(A*B-B*A)
    pairs = []
    for a, multiplicity, basis in A.eigenvects():
        V = s.Matrix.hstack(*basis)
        assert V.cols == multiplicity
        W, parameters = V.gauss_jordan_solve(B*V)
        assert parameters.rows == 0
        assert zero(V*W-B*V)
        for b, count in W.eigenvals().items():
            pairs.extend([(a, b)] * count)
    assert len(pairs) == A.rows
    return pairs


def matrix_rhs(A, B):
    n, M = A.rows, A-B
    if n == 1:
        return True
    if n == 2:
        if not zero(A*B-B*A):
            return False
        pairs = joint_pairs(A, B)
        return diagonal_rhs([x[0] for x in pairs], [x[1] for x in pairs])
    if zero(M):
        return True
    if M.rank() != 1:
        return False
    z = M.columnspace()[0]
    j = next(j for j in range(n) if z[j] != 0)
    a, b = s.cancel((A*z)[j]/z[j]), s.cancel((B*z)[j]/z[j])
    if not zero(A*z-a*z) or not zero(B*z-b*z):
        return False
    norm2 = (z.H*z)[0]
    assert norm2 > 0
    if not zero(M-(a-b)*(z*z.H)/norm2):
        return False
    common_a, common_b = eigenlist(A), eigenlist(B)
    common_a.remove(a)
    common_b.remove(b)
    assert Counter(common_a) == Counter(common_b)
    return separated(a, b, common_a)


def spectral_comparison(M, match2):
    """Compare the literal top-m Gram eigenvalue sum with a rational number.

    Rational interval arithmetic certifies strictness; singleton intervals
    certify equality. An unresolved comparison raises rather than using floats.
    """
    n, m = M.rows, (M.rows+1)//2
    x = s.Symbol('x')
    gram = M.H*M
    assert zero(gram-gram.H)
    p = s.Poly(gram.charpoly(x).as_expr(), x, domain=s.QQ)
    for digits in (8, 16, 32, 64):
        intervals = []
        for (lo, hi), multiplicity in p.intervals(eps=R(1, 10**digits)):
            assert lo >= 0, (p, lo, hi)
            intervals.extend([(lo, hi)] * multiplicity)
        assert len(intervals) == n
        top = intervals[-m:]
        lo, hi = sum(t[0] for t in top), sum(t[1] for t in top)
        if lo > match2:
            return False, (lo, hi), p.as_expr()
        if lo == hi == match2:
            return True, (lo, hi), p.as_expr()
        if hi < match2:
            raise AssertionError(('source inequality failed', match2, lo, hi, p))
    raise AssertionError(('exact comparison unresolved', p, match2, lo, hi))


matrix_counts = Counter()


def check_matrix(label, A, B, verbose=False, expected=None):
    assert A.rows == A.cols == B.rows == B.cols
    assert zero(A*A.H-A.H*A) and zero(B*B.H-B.H*B), label
    alpha, beta = eigenlist(A), eigenlist(B)
    # These fixtures have Gaussian-rational spectra; matching is exact rational arithmetic.
    match2 = d2(alpha, beta)
    is_equal, interval, polynomial = spectral_comparison(A-B, match2)
    prediction = matrix_rhs(A, B)
    assert is_equal == prediction, (label, A, B, match2, interval, polynomial)
    if expected is not None:
        assert is_equal == expected, label
    matrix_counts[(A.rows, 'equality' if is_equal else 'strict')] += 1
    if verbose:
        print(f'MATRIX {label}: n={A.rows} rank_difference={(A-B).rank()} '
              f'commute={zero(A*B-B*A)} d_squared={match2} '
              f'S_m_interval=[{interval[0]},{interval[1]}] equality={is_equal}')
    return is_equal


def rotation(n, i, j, t):
    U = s.eye(n)
    c, q = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    U[i, i] = U[j, j] = c
    U[i, j], U[j, i] = -q, q
    assert zero(U.H*U-s.eye(n))
    return U


def unitary(n, rng):
    U = s.eye(n)
    for _ in range(3):
        i, j = rng.sample(range(n), 2)
        U = U * rotation(n, i, j, rng.choice([R(1, 2), R(1, 3), R(2, 3), R(1)]))
    U = U * s.diag(*(I**rng.randrange(4) for _ in range(n)))
    assert zero(U.H*U-s.eye(n))
    return U


def structured():
    # Source Section 1: A=diag(t,0,...,0), B=0, t>0 gives d^2=S_m=t^2.
    for n in range(1, 5):
        for t in (R(1, 3), R(2)):
            A, B = s.diag(t, *([0]*(n-1))), s.zeros(n)
            assert d2(eigenlist(A), eigenlist(B)) == t*t
            assert check_matrix(f'source_rank_one_t={t}', A, B, True, True)
    # Dimension two permits rank-two equality and has a separate rematching branch.
    check_matrix('n2_rank_two_equality', s.diag(0, 2), s.diag(1, 3), True, True)
    check_matrix('n2_crossed_improves', s.diag(0, 1), s.diag(2, 1), True, False)
    eps = R(1, 100)
    for n in (2, 3, 4):
        tail = list(range(30, 30+n-2))
        A, B = s.diag(0, 10, *tail), s.diag(1, 10, *tail)
        check_matrix('separated_equality', A, B, True, True)
        U = rotation(n, 0, 1, eps)
        check_matrix('near_equality_noncommuting', A, U*B*U.H, True, False)
        # Rational-entry, noncommuting, rank-one difference; spectra {0,3} and {1,6}.
        C = s.diag(0, 3, *tail)
        D = s.diag(2, 5, *tail)
        D[0, 1] = D[1, 0] = 2
        check_matrix('noncommuting_rank_one', C, D, True, False)
        if n >= 3:
            check_matrix('rank_two_near_equality', A, s.diag(1, 10+eps, *tail), True, False)
            a, b = [0, 1-eps, 3] + tail[:n-3], [1, 1-eps, 3] + tail[:n-3]
            check_matrix('connected_rank_one_near_equality', s.diag(*a), s.diag(*b), True, False)
            check_matrix('identical_normal_pair', A, A, True, True)
    # Repeated common values retain distinct labels in the threshold graph.
    check_matrix('repeated_common_labels', s.diag(0, 3, 3, 3), s.diag(1, 3, 3, 3), True, True)
    check_matrix('rank_one_path_0_1_2', s.diag(0, 1, 2), s.diag(2, 1, 2), True, False)


def random_normal():
    rng = Random(SEED)
    for n in (2, 3, 4):
        for k in range(12):
            alpha = [s.Integer(rng.randrange(-3, 4)) + I*rng.randrange(-2, 3) for _ in range(n)]
            beta = alpha[:] if k == 0 else [
                s.Integer(rng.randrange(-3, 4)) + I*rng.randrange(-2, 3) for _ in range(n)]
            U, V = unitary(n, rng), unitary(n, rng)
            A, B = U*s.diag(*alpha)*U.H, V*s.diag(*beta)*V.H
            check_matrix(f'random_independent_n{n}_{k}', A, B)
        for k in range(4):
            U = unitary(n, rng)
            alpha = [s.Integer(3*j) + I*j for j in range(n)]
            beta = alpha[:]
            beta[0] += 1
            check_matrix(f'random_common_basis_n{n}_{k}',
                         U*s.diag(*alpha)*U.H, U*s.diag(*beta)*U.H, expected=True)
        print(f'RANDOM_NORMAL n={n} independent_bases=12 common_basis=4 seed={SEED}')


def main():
    print(f'Python {sys.version_info.major}.{sys.version_info.minor}; SymPy {s.__version__}; exact arithmetic only')
    exhaustive_diagonal()
    structured()
    random_normal()
    for n in range(1, 5):
        print(f'MATRIX_TOTAL n={n} equality={matrix_counts[n,"equality"]} strict={matrix_counts[n,"strict"]}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
