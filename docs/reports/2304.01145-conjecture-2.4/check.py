#!/usr/bin/env python3
"""Exact checks for Athreya--Mukherjee--Mukherjee Conjecture 2.4.
All numerical evidence is integer/Fraction/SymPy exact arithmetic; no floats.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
import sympy as sp


def passage(n: int, r: int, s: int) -> list[Fraction]:
    """e_j = expected time to go from overlap j to j+1."""
    e = [Fraction(n - r, s)]
    for j in range(1, s):
        e.append(Fraction(r * (n - r), (r - j) * (s - j))
                 + Fraction(j * (n - r - s + j), (r - j) * (s - j)) * e[-1])
    return e


def direct_overlap(n: int, r: int, s: int) -> list[Fraction]:
    """Solve first-step equations for E_j, the time to hit overlap s."""
    den = r * (n - r)
    rows = []
    rhs = []
    for j in range(s):
        up = Fraction((r - j) * (s - j), den)
        down = Fraction(j * (n - r - s + j), den)
        stay = 1 - up - down
        row = [Fraction(0) for _ in range(s)]
        row[j] = 1 - stay
        if j + 1 < s:
            row[j + 1] -= up
        if j > 0:
            row[j - 1] -= down
        rows.append(row)
        rhs.append(Fraction(1))
    M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in rows])
    b = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in rhs])
    sol = M.inv() * b
    return [Fraction(int(x.p), int(x.q)) for x in sol]


def log_upper(q: Fraction, terms: int = 80) -> Fraction:
    """Rational upper bound: log(q) = 2 atanh(x), geometric tail bound."""
    assert q >= 1
    x = (q - 1) / (q + 1)
    partial = 2 * sum((x ** (2*k+1)) / (2*k+1) for k in range(terms))
    remainder = 2 * x ** (2*terms+1) / ((2*terms+1) * (1-x*x))
    return partial + remainder


def verify_actual_subset_chain(n: int, r: int, s: int, e: list[Fraction]) -> int:
    """Enumerate actual replacements, verify lumping and exact hitting residuals."""
    target = set(range(s))
    hit = [sum(e[j:], Fraction(0)) for j in range(s)] + [Fraction(0)]
    count = 0
    for H0 in combinations(range(n), r):
        H = set(H0)
        j = len(H & target)
        overlaps = [0] * (s+1)
        for a in H:
            for b in set(range(n)) - H:
                k = len(((H - {a}) | {b}) & target)
                overlaps[k] += 1
        den = r * (n-r)
        assert sum(overlaps) == den
        if j < s:
            assert overlaps[j+1] == (r-j)*(s-j)
            if j: assert overlaps[j-1] == j*(n-r-s+j)
            assert hit[j] == 1 + sum((Fraction(c,den)*hit[k] for k,c in enumerate(overlaps)), Fraction(0))
        count += 1
    return count


def check_tied_records() -> int:
    """Enumerate small groupings/permutations and certify the record inequality."""
    cases = 0
    for m in range(2, 7):
        # All ordered compositions via cuts, first group is discovered at time zero.
        for cuts in product((False, True), repeat=m-1):
            groups = [0]
            for cut in cuts: groups.append(groups[-1] + int(cut))
            sizes = [groups.count(i) for i in range(groups[-1]+1)]
            C = max(sizes)
            expected = sum((Fraction(sizes[i], sum(sizes[i:])) for i in range(1,len(sizes))), Fraction(0))
            records = 0
            for perm in permutations(range(m)):
                highest = 0
                for a in perm:
                    g = groups[a]
                    if g > highest: records += 1; highest = g
            from math import factorial
            assert Fraction(records, factorial(m)) == expected
            assert expected <= harmonic(m)
            if m > C:
                assert log_upper(Fraction(m,C)) <= expected
            cases += 1
    return cases


def harmonic(m: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, m + 1)), Fraction(0))


def label_masks(n: int, r: int, s: int):
    sets = list(combinations(range(n), r))
    labels = list(combinations(range(n), s))
    index = {a: i for i, a in enumerate(labels)}
    hm = []
    for H in sets:
        hs = set(H)
        mask = 0
        for i, A in enumerate(labels):
            if set(A).issubset(hs):
                mask |= 1 << i
        hm.append(mask)
    return sets, labels, hm, index


def exact_cover_expectation(n: int, r: int, s: int) -> Fraction:
    """Exact E[steps after time 0 to collect all labels], by mask layers."""
    sets, labels, hmask, _ = label_masks(n, r, s)
    full = (1 << len(labels)) - 1
    # All masks containing at least one state label; solve each mask layer.
    masks = sorted({m for m in range(full + 1) if any((hm & ~m) == 0 for hm in hmask)},
                   key=lambda m: (m.bit_count(), m), reverse=True)
    values: dict[tuple[int, int], Fraction] = {}
    for mask in masks:
        if mask == full:
            continue
        hs = [i for i, hm in enumerate(hmask) if hm & ~mask == 0]
        pos = {i: k for k, i in enumerate(hs)}
        if not hs:
            continue
        A = [[Fraction(0) for _ in hs] for _ in hs]
        b = [Fraction(1) for _ in hs]
        for i in hs:
            row = pos[i]
            A[row][row] = 1
            H = set(sets[i])
            for x in H:
                for y in range(n):
                    if y in H:
                        continue
                    H2 = tuple(sorted((H - {x}) | {y}))
                    j = sets.index(H2)
                    m2 = mask | hmask[j]
                    p = Fraction(1, r * (n - r))
                    if m2 == mask:
                        A[row][pos[j]] -= p
                    elif m2 != full:
                        b[row] += p * values[(m2, j)]
        M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in A])
        B = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in b])
        sol = M.inv() * B
        for i, x in zip(hs, sol):
            values[(mask, i)] = Fraction(int(x.p), int(x.q))
    # Start H is uniform; mask is exactly its collected labels.
    ans = sum((Fraction(0) if hmask[i] == full else values[(hmask[i], i)] for i in range(len(sets))), Fraction(0)) / len(sets)
    return ans


def main() -> None:
    checked = 0
    subset_states = 0
    passage_rows = []
    # All overlap levels exist when n >= r+s; smaller phantom levels are excluded.
    for n in range(4, 13):
        for r, s in ((2, 1), (3, 2), (3, 1)):
            if n < r+s or s > r:
                continue
            e = passage(n, r, s)
            direct = direct_overlap(n, r, s)
            assert e == [direct[j] - (direct[j + 1] if j + 1 < s else Fraction(0)) for j in range(s)]
            assert all(x > 0 for x in e)
            subset_states += verify_actual_subset_chain(n,r,s,e)
            passage_rows.append((n,r,s,e[-1],e[-1]/n**s))
            checked += 1
    # Symbolic leading coefficient in n, exactly as in Lemma 2.
    n = sp.symbols("n")
    for r, s in ((2, 1), (3, 2), (3, 1), (4, 2)):
        e = [sp.Rational(1, s) * (n - r)]
        for j in range(1, s):
            e.append(sp.cancel(sp.Rational(r, (r-j)*(s-j)) * (n-r)
                               + sp.Rational(j, (r-j)*(s-j)) * (n-r-s+j) * e[-1]))
        lead = sp.LC(sp.Poly(e[-1], n))
        expected = sp.Rational(1, s)
        for i in range(1, s):
            expected *= sp.Rational(i, (r-i)*(s-i))
        assert lead == expected
    # Independent small-state cover computations; test the record-cover squeeze.
    cover_rows = []
    for n0 in range(3, 6):
        for r, s in ((2, 1), (3, 2)):
            if r >= n0:
                continue
            actual = exact_cover_expectation(n0, r, s)
            m, C = comb(n0, s), comb(r, s)
            e = passage(n0, r, s)
            hminus, hplus = e[-1], sum(e, Fraction(0))
            lower = hminus * log_upper(Fraction(m, C))
            upper = hplus * harmonic(m)
            assert lower <= actual <= upper, (n0, r, s, lower, actual, upper)
            cover_rows.append((n0, r, s, actual, lower, upper))
    tied_cases = check_tied_records()
    # Source fidelity: r=s=2 is the Johnson graph, 2(n-2)-regular on C(n,2) vertices.
    for n0 in range(3, 8):
        r = s = 2
        vertices = list(combinations(range(n0), r))
        degree = sum(1 for V in vertices if len(set(vertices[0]) & set(V)) == r - 1)
        assert degree == r * (n0 - r)
    print(f"recurrence_cases={checked}")
    print(f"actual_subset_states_checked={subset_states}")
    print(f"tied_group_compositions_checked={tied_cases}")
    print("symbolic_leading_coefficients={(2,1):1, (3,2):1/4, (3,1):1, (4,2):1/6}")
    for n0,r,s,e,ratio in passage_rows:
        print(f"passage n={n0} r={r} s={s}: e_last={e}; e_last/n^s={ratio}; exact_state_residuals=0")
    for row in cover_rows:
        n0, r, s, actual, lower, upper = row
        print(f"cover n={n0} r={r} s={s}: E_steps={actual}; exact record-squeeze inequality verified (certified_log_upper_lower_bound <= E_steps <= {upper})")
    print("fidelity_johnson_graph_degree=verified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
