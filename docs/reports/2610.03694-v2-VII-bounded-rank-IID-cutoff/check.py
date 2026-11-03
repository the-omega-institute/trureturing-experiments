#!/usr/bin/env python3
"""Exact preregistration checks for arXiv:2610.03694v2, Sections II–III and VII.

The universal cutoff is NOT proved by these finite checks. The delay families
are linear realizations with measurement completeness, not certified GPTs.
Run with Python 3 and SymPy. No floating-point arithmetic is used.
"""
from fractions import Fraction as F
from itertools import product
import sympy as sp


def smat(rows):
    return sp.Matrix(rows)


J = sp.diag(1, -1, 1, 1)
P = sp.diag(0, 0, 0, 1)
KS = [sp.eye(4), sp.diag(1, -1, -1, 1),
      sp.diag(-1, 1, -1, 1), sp.diag(-1, -1, 1, 1)]


def targets(c, s):
    rz = smat([[c, s, 0, 0], [-s, c, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    rx = smat([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])
    return [[k * rot * J for k in KS] for rot in [sp.eye(4), rz, rx]]


def fidelity():
    """Derive equation (1) by contracting two Bell states with Erwin's effect."""
    pauli = [sp.eye(2), smat([[0, 1], [1, 0]]),
             smat([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
    obs = pauli[1:] + [pauli[0]]
    xi = [sp.eye(2), (3 * sp.eye(2) + sp.I * pauli[3]) / sp.sqrt(10),
          (3 * sp.eye(2) + sp.I * pauli[1]) / sp.sqrt(10)]
    cs = targets(sp.Rational(4, 5), sp.Rational(3, 5))
    for r, mu in product(range(3), range(4)):
        measurement_ket = (pauli[mu] * xi[r].conjugate().T).reshape(4, 1) / sp.sqrt(2)
        # Phi(A,E1) tensor Phi(E2,B) contracted with the bra on E1,E2.
        post = measurement_ket.conjugate() / 2
        probability = sp.simplify((post.conjugate().T * post)[0])
        assert probability == sp.Rational(1, 4)
        state = post * post.conjugate().T / probability
        corr = sp.Matrix(4, 4, lambda i, j: sp.simplify(
            sp.trace(state * sp.kronecker_product(obs[i], obs[j]))))
        assert corr == cs[r][mu], (r, mu, corr, cs[r][mu])
    assert sum(KS, sp.zeros(4)) == 4 * P
    assert sp.Matrix.hstack(*(sp.Matrix(k).reshape(16, 1) for k in KS)).rank() == 4
    print('FIDELITY: equation (1), all 12 setting/outcome pairs; Bell probabilities = 1/4')
    print('FIDELITY: four Pauli sign matrices independent; sum K_mu = 4P')


def rational_matrix(m):
    return tuple(tuple(F(x) for x in row) for row in m.tolist())


def identity(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def multiply(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(len(b))
                            if a[i][k] and b[k][j]), F(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def all_words(n):
    for outcomes in product(range(4), repeat=n):
        yield tuple((0, mu) for mu in outcomes), None
        for pos in range(n):
            for r in (1, 2):
                yield tuple((r if k == pos else 0, mu) for k, mu in enumerate(outcomes)), pos


def restricted_products():
    """All allowed words through length 4 for two exact rational rotations."""
    count = 0
    for c, s in [(sp.Rational(4, 5), sp.Rational(3, 5)),
                 (sp.Rational(12, 13), sp.Rational(5, 13))]:
        assert c * c + s * s == 1 and 0 < s <= c
        cs = targets(c, s)
        q = [[rational_matrix(cm.T * J) for cm in row] for row in cs]
        k = [rational_matrix(x) for x in KS]
        jmat = rational_matrix(J)
        for n in range(5):
            for word, pos in all_words(n):
                actual = jmat
                for r, mu in word:
                    actual = multiply(q[r][mu], actual)
                if not n:
                    expected = jmat
                else:
                    pivot = 0 if pos is None else pos
                    left = identity(4)
                    right = identity(4)
                    for _, mu in word[:pivot]:
                        left = multiply(left, k[mu])
                    for _, mu in word[pivot + 1:]:
                        right = multiply(right, k[mu])
                    r, mu = word[pivot]
                    expected = multiply(multiply(left, rational_matrix(cs[r][mu])), right)
                # Readout matrix has input index first, hence transpose.
                assert actual == tuple(zip(*expected)), (c, s, word)
                assert actual[3][3] == 1
                prob = F(1, 4 ** n)
                assert prob > 0 and 4 ** n * prob == 1
                count += 1
        # Negative control: reversed chronological order is distinguishable.
        assert multiply(q[0][1], multiply(q[1][0], jmat)) != multiply(
            q[1][0], multiply(q[0][1], jmat))
    print(f'PRODUCTS: {count} words, lengths 0..4, two rotations; targets and p_w = 4^(-N) agree')
    print('NEGATIVE_CONTROL: reversing map order fails an allowed length-2 word')


def finite_normalization():
    """Each finite induction step uniquely determines the four new probabilities."""
    count = 0
    cs = targets(sp.Rational(4, 5), sp.Rational(3, 5))
    for n, r, mu, side in product(range(4), range(3), range(4), ('left', 'right')):
        cmat = cs[r][mu]
        matrices = [k * cmat if side == 'left' else cmat * k for k in KS]
        a = sp.Matrix.hstack(*(m.reshape(16, 1) for m in matrices))
        rhs = (sp.Rational(1, 4 ** n) * P).reshape(16, 1)
        solution = sp.ones(4, 1) / (4 ** (n + 1))
        assert a.rank() == 4
        assert a * solution == rhs
        count += 1
    print(f'NORMALIZATION: {count} exact, uniquely solvable endpoint-padding steps through length 4')


def basis(m):
    cols = m.columnspace()
    return sp.Matrix.hstack(*cols) if cols else sp.zeros(m.rows, 0)


def exact_span_scan(ts, seeds, readout, max_length):
    """Span all base words and all one-marked words; enumerate linear bases exactly."""
    n = ts[0][0].rows
    ls = [[sp.diag(t, cm.T * J) for t, cm in zip(row, crow)]
          for row, crow in zip(ts, targets(sp.Rational(4, 5), sp.Rational(3, 5)))]
    z = seeds.col_join(J)
    f = readout.row_join(-sp.eye(4))
    base, marked = basis(z), sp.zeros(n + 4, 0)
    first = None
    for length in range(max_length + 1):
        if first is None and (f * base != sp.zeros(4, base.cols)
                              or f * marked != sp.zeros(4, marked.cols)):
            first = length
        if length == max_length:
            break
        old_base = base
        base = basis(sp.Matrix.hstack(base, *(a * base for a in ls[0])))
        marked = basis(sp.Matrix.hstack(marked, *(a * marked for a in ls[0]),
                                       *(a * old_base for row in ls[1:] for a in row)))
    return first, base.cols, marked.cols


def realization_family(d, epsilon, delayed):
    """Rank-d rhohat = identity and unit e_3, with sum_mu T_(r,mu) = 4 e_3 e_3^T.

    delayed=False: arbitrary hidden dynamics and one-way hidden injection.
    delayed=True: a hidden delay line feeding back at the marked letter.
    Positivity on GPT cones is not asserted for these algebraic models.
    """
    cs = targets(sp.Rational(4, 5), sp.Rational(3, 5))
    seeds = J.col_join(sp.zeros(d - 4, 4))
    readout = sp.eye(4).row_join(sp.zeros(4, d - 4))
    ts = []
    for r in range(3):
        row = []
        for mu in range(4):
            t = sp.zeros(d)
            t[:4, :4] = cs[r][mu].T * J
            sign = KS[mu][0, 0]  # A nontrivial character, summing to zero.
            if delayed:
                if d == 4:
                    if r == 1:
                        t[0, 0] += epsilon * sign
                else:
                    if r == 0:
                        t[4, 0] = sign
                        for k in range(4, d - 1):
                            t[k + 1, k] = sign
                    if r == 1:
                        t[0, d - 1] = epsilon * sign
            else:
                for i in range(4, d):
                    t[i, 0] = epsilon * sign * (i + r)
                    t[i, i] = epsilon * sign * (r + 1)
            row.append(t)
        complete = sp.zeros(d)
        complete[3, 3] = 4
        assert sum(row, sp.zeros(d)) == complete
        ts.append(row)
    assert sp.eye(d).rank() == d
    assert readout * seeds == J
    return ts, seeds, readout


def cutoff_cases():
    cases = 0
    sharpness = []
    for d in range(4, 10):
        for delayed, epsilon in product((False, True), (-1, 0, 1)):
            ts, seeds, readout = realization_family(d, epsilon, delayed)
            first, u, v = exact_span_scan(ts, seeds, readout, 2 * d)
            assert first is None or first <= d - 3
            if delayed and epsilon:
                assert first == d - 3
                sharpness.append((d, first))
            else:
                assert first is None
            cases += 1
    print(f'CUTOFF_CASES: {cases} complete exact linear datasets, ranks 4..9, scanned through 2D')
    print('CUTOFF_CASES: every dataset passing through D-3 also passed the full scan')
    print('TRIAGE_ALGEBRAIC_SHARPNESS: first failure D-3 for D=4..9 and epsilon=-1,+1')
    print('TRIAGE_BOUNDARY: these delay lines are not certified GPT realizations; physical sharpness is open')


def main():
    fidelity()
    restricted_products()
    finite_normalization()
    cutoff_cases()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
