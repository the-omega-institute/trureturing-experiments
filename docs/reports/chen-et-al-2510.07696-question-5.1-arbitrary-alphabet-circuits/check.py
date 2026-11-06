#!/usr/bin/env python3
"""Exact finite checks for arXiv:2510.07696v1, Question 5.1.

This is experimental evidence for the nominated argument, not a proof of the
asymptotic circuit lower bound. Python's integers and Fraction are the only
numeric evidence. No third-party packages or random samples are used.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product


def exactly_accepts(transitions, initial, final, word, alphabet):
    """Track whether an accepting word of the given length diverges from word.

    Multiple paths reading the same word remain allowed. State and divergence
    are propagated together; thus no list of all same-length words is needed.
    """
    current = {(initial, False)}
    for target in word:
        current = {
            (dest, divergent or letter != target)
            for state, divergent in current
            for letter in alphabet
            for dest in transitions.get((state, letter), ())
        }
    return (final, False) in current and (final, True) not in current


def fidelity():
    # Source Proposition 2.1: for nonempty x, m > |x|/q implies x^m in L_q.
    # Its displayed cycle is checked literally; this verifies the stated upper
    # bound, rather than assuming that this particular NFA is state-minimal.
    count = 0
    for r in (2, 3, 4):
        alphabet = range(r)
        for length in range(1, 4):
            for x in product(alphabet, repeat=length):
                transitions = {(i, a): {(i + 1) % length} for i, a in enumerate(x)}
                for q in (Fraction(1, 4), Fraction(1, 3), Fraction(2, 5)):
                    m = length // q + 1
                    assert m > Fraction(length, 1) / q
                    w = x * m
                    assert exactly_accepts(transitions, 0, 0, w, alphabet)
                    assert length < q * len(w)
                    count += 1
    # Nonemptiness supplier, including strict-threshold equality controls.
    for r in (2, 3, 4):
        for a in range(r):
            for q in (Fraction(1, 4), Fraction(1, 3), Fraction(2, 5)):
                for n in range(13):
                    assert exactly_accepts({(0, a): {0}}, 0, 0, (a,) * n, range(r))
                    assert (Fraction(1) < q * n) == (n > Fraction(1) / q)
    print(f"FIDELITY: Proposition 2.1 cycle witnesses, {count} exact cases; nonempty x of lengths 1..3, alphabet sizes 2..4, q=1/4,1/3,2/5")
    print("NONEMPTINESS: one-state constant-word witnesses and strict-threshold controls passed")


def exponent_vectors(n, degree):
    return [e for e in product(range(degree + 1), repeat=n) if sum(e) <= degree]


def grid_support():
    # Exhaustive coefficient vectors, not random testing. Unreduced polynomials
    # are included: e.g. X^2-X can vanish on a grid without being syntactically 0.
    cases = [(2, 0, 3), (2, 1, 4), (2, 2, 3), (2, 3, 2),
             (3, 0, 2), (3, 1, 4), (3, 2, 2),
             (5, 0, 2), (5, 1, 3), (5, 2, 1)]
    total_polynomials = total_grids = nonzero_grids = vanishing_grids = 0
    for p, n, cap in cases:
        monomials = exponent_vectors(n, cap)
        grids = []
        for r in range(2, p + 1):
            for s in combinations(range(p), r):
                points = list(product(s, repeat=n))
                matrix = [
                    [prod_mod((pow(x, a, p) for x, a in zip(point, exp)), p)
                     for exp in monomials]
                    for point in points
                ]
                grids.append((r, matrix))
        count = supported = 0
        for coefficients in product(range(p), repeat=len(monomials)):
            degree = max((sum(e) for a, e in zip(coefficients, monomials) if a), default=0)
            count += 1
            for r, matrix in grids:
                values = [sum(a * m for a, m in zip(coefficients, row)) % p for row in matrix]
                support = sum(v != 0 for v in values)
                total_grids += 1
                if support:
                    # |support| / r^n >= r^(-degree), with integer cross-multiplication.
                    assert support * r ** degree >= r ** n
                    nonzero_grids += 1
                    supported += 1
                else:
                    vanishing_grids += 1
                # Both directions of the complement zero-test, for every point.
                transformed = [(1 - pow(v, p - 1, p)) % p for v in values]
                assert all((u != 0) == (v == 0) for u, v in zip(transformed, values))
        total_polynomials += count
        print(f"GRID: F_{p}, n={n}, totalDegree<={cap}, polynomials={count}, grids-per-polynomial={len(grids)}, nonzero-grid-pairs={supported}")
    assert vanishing_grids > 0 and nonzero_grids > 0
    # Outside-grid negative control: X(X-1) vanishes on {0,1}, not on all F5.
    assert [x * (x - 1) % 5 for x in (0, 1)] == [0, 0]
    assert 2 * (2 - 1) % 5 != 0
    print(f"GRID_TOTAL: polynomials={total_polynomials}, polynomial-grid-pairs={total_grids}, nonzero={nonzero_grids}, vanishing={vanishing_grids}")
    print("COMPLEMENT: 1-P^(p-1) nonzero iff P=0 at every enumerated grid point")


def prod_mod(values, p):
    answer = 1
    for v in values:
        answer = answer * v % p
    return answer


def restricted_growth_sequences(length):
    def extend(prefix):
        if len(prefix) == length:
            yield tuple(prefix)
        else:
            for state in range(max(prefix) + 2):
                yield from extend(prefix + [state])
    yield from extend([0])


@lru_cache(None)
def exact_complexity(word, r):
    # Any exactly accepting NFA has an accepting path. Delete other edges and
    # final states: exact acceptance persists. Canonical renamings of those
    # paths are precisely the restricted growth sequences enumerated here.
    answer = len(word) + 1
    for path in restricted_growth_sequences(len(word) + 1):
        states = max(path) + 1
        if states >= answer:
            continue
        transitions = {}
        for s, a, t in zip(path, word, path[1:]):
            transitions.setdefault((s, a), set()).add(t)
        if exactly_accepts(transitions, 0, path[-1], word, range(r)):
            answer = states
    return answer


def nfa_and_boolean_restriction():
    words = pairs = literal_checks = 0
    for n in range(5):
        for word in product(range(2), repeat=n):
            base = exact_complexity(word, 2)
            for r in (2, 3, 4):
                for a, b in permutations(range(r), 2):
                    renamed = tuple((a, b)[i] for i in word)
                    assert exact_complexity(renamed, r) == base
                    words += 1
    # All ordered codeword pairs in blocks of 1..3 bits: every encoded bit is
    # a constant, z or 1-z. There is no condition on invalid block encodings.
    for t in (1, 2, 3):
        codes = list(product(range(2), repeat=t))
        for a, b in permutations(codes, 2):
            pairs += 1
            for i in range(t):
                alpha, beta = a[i], b[i] - a[i]
                assert (alpha, beta) in {(0, 0), (1, 0), (0, 1), (1, -1)}
                for z in (0, 1):
                    assert alpha + beta * z == (a, b)[z][i]
                    literal_checks += 1
    # One small ambiguous-path NFA: accepts 00 along two paths. This must count
    # as exact acceptance, rather than being rejected by a unique-path convention.
    ambiguous = {(0, 0): {0, 1}, (1, 0): {1}}
    assert exactly_accepts(ambiguous, 0, 1, (0, 0), range(2))
    print(f"ALPHABET_INVARIANCE: {words} exhaustive binary-word/ordered-subalphabet comparisons, lengths 0..4, ambient sizes 2..4")
    print(f"BOOLEAN_BLOCKS: {pairs} distinct ordered codeword pairs, block widths 1..3, {literal_checks} affine-literal evaluations")
    print("NFA_SEMANTICS: epsilon-free, one initial state; multiple accepting paths allowed")


def main():
    fidelity()
    grid_support()
    nfa_and_boolean_restriction()
    print("LIMIT: finite checks validate the listed suppliers, not the density limit or the unbounded circuit lower bound")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
