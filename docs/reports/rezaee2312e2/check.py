#!/usr/bin/env python3
"""Exact checks for Rezaee, arXiv:2312.04520v2, Conjecture E(2).

Requires Python >= 3.12 and Singular. No third-party program source is copied.
Tangent dimensions are computed twice: Taylor syzygy coefficient equations
over Q and Singular's independently computed syzygies over Q[x,y,z]/I.
These rational matrices have the same rank after scalar extension to C.
Finite checks do not prove the conjecture or its global maximality.
"""

from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import re
import subprocess
import sys
import tempfile


assert sys.version_info >= (3, 12)


def degree(a):
    return sum(a)


def monomials(d):
    return [(a, b, d - a - b) for a in range(d + 1) for b in range(d - a + 1)]


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def contains(gens, a):
    return any(divides(g, a) for g in gens)


def minimum(gens):
    return sorted(g for g in set(gens) if not any(h != g and divides(h, g) for h in gens))


def basis(gens):
    bounds = [min(g[i] for g in gens if all(g[j] == 0 for j in range(3) if j != i))
              for i in range(3)]
    return sorted(a for a in product(*(range(b) for b in bounds)) if not contains(gens, a))


def source_ideal(k):
    return minimum([a for a in monomials(k) if a[:2] != (0, 0)] + [(0, 0, k + 2)])


def canonical_ideal(k):
    b = set(a for d in range(k) for a in monomials(d)) | {(0, 0, k), (0, 1, k - 1)}
    border = {tuple(a[i] + (i == j) for i in range(3)) for a in b for j in range(3)} - b
    gens = minimum(border)
    assert set(basis(gens)) == b
    return gens


def tangent(gens):
    """Dimension of Hom(I,R/I), including all homogeneous degrees.

    Pair syzygies (lcm/gi) ei - (lcm/gj) ej generate the relation module
    for a monomial ideal. Each coefficient equation is u=v or u=0.
    A disjoint-set computation solves this exact integer equation system;
    it does not use a numerical rank threshold or modular rank estimate.
    """
    b = basis(gens)
    indices = {(i, a): n for n, (i, a) in enumerate(product(range(len(gens)), b), 1)}
    parent = list(range(len(indices) + 1))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def join(i, j):
        i, j = root(i), root(j)
        if i != j:
            parent[max(i, j)] = min(i, j)

    for i, j in combinations(range(len(gens)), 2):
        lcm = tuple(max(x, y) for x, y in zip(gens[i], gens[j]))
        terms = []
        for n in [i, j]:
            shift = tuple(x - y for x, y in zip(lcm, gens[n]))
            terms.append({tuple(x + y for x, y in zip(a, shift)): indices[n, a]
                          for a in b if tuple(x + y for x, y in zip(a, shift)) in b})
        for a in terms[0].keys() | terms[1].keys():
            join(terms[0].get(a, 0), terms[1].get(a, 0))
    classes = {}
    for (i, a), n in indices.items():
        r = root(n)
        if r != root(0):
            d = degree(a) - degree(gens[i])
            assert r not in classes or classes[r] == d
            classes[r] = d
    return len(classes), dict(sorted(Counter(classes.values()).items()))


def singular_monomial(a):
    return "*".join(f"{v}^{n}" for v, n in zip("xyz", a) if n) or "1"


def singular_check(gens):
    """Compute the kernel of the transposed presentation over the quotient.

    Singular's syz is over the quotient ring A, and vdim(A^g / kernel)
    gives g*dim(A) - dim(kernel). This is a second computation, not a
    reimplementation of the Taylor-pair equations used above.
    """
    program = f"""ring r=0,(x,y,z),dp;
ideal I={','.join(map(singular_monomial, gens))};
matrix T=transpose(syz(I));
int n=vdim(std(I));
int g=size(I);
qring A=std(I);
matrix M=imap(r,T);
module K=syz(M);
print("RESULT="+string(n)+","+string(g*n-vdim(std(K))));
quit;
"""
    with tempfile.TemporaryDirectory(prefix="singular-", dir=Path(__file__).parent) as td:
        p = Path(td) / "check.sing"
        p.write_text(program)
        out = subprocess.run(["Singular", "-q", str(p)], text=True, capture_output=True,
                             check=True, timeout=120)
    assert not out.stderr, out.stderr
    assert "?" not in out.stdout, out.stdout
    m = re.fullmatch(r"RESULT=(\d+),(\d+)\s*", out.stdout)
    assert m, out.stdout
    return tuple(map(int, m.groups()))


def staircase_candidates(k):
    """All Borel staircases with colength L_k+2 and smallest x exponent k.

    Borel-fixedness forces every monomial of degree <k to be standard.
    Exactly two additional monomials are standard. Order-ideal closure
    bounds their degree by k+1. We enumerate that finite region exactly.
    A standard staircase is closed under replacing a smaller-index
    variable by a larger-index variable (the complement of strong stability).
    """
    low = set(a for d in range(k) for a in monomials(d))
    region = monomials(k) + monomials(k + 1)
    answer = []
    for extra in combinations(region, 2):
        b = low | set(extra)
        if (k, 0, 0) in b:
            continue
        good = True
        for a in extra:
            for i in range(3):
                if not a[i]:
                    continue
                prev = list(a)
                prev[i] -= 1
                if tuple(prev) not in b:
                    good = False
                for j in range(i + 1, 3):
                    moved = prev.copy()
                    moved[j] += 1
                    if tuple(moved) not in b:
                        good = False
        if good:
            answer.append(tuple(sorted(extra)))
    expected = {tuple(sorted([(0, 0, k), (0, 0, k + 1)])),
                tuple(sorted([(0, 0, k), (0, 1, k - 1)]))}
    assert set(answer) == expected, (k, answer)
    return answer


def main():
    print("Exact arithmetic: integer coefficient equations and Singular characteristic 0")
    print("Source fidelity: Corollary 1.9 and the displayed k=3 E(2) example")
    for k in range(1, 7):
        fat = monomials(k)
        length = comb(k + 2, 3)
        formula = comb(k + 1, 2) * comb(k + 2, 2)
        assert len(basis(fat)) == length
        assert tangent(fat)[0] == formula
        assert singular_check(fat) == (length, formula)
        print(f"FIDELITY m^{k}: colength={length} tangent={formula} dual_check=OK")
    example = [(3, 0, 0), (0, 3, 0), (0, 0, 5), (0, 1, 2), (1, 0, 2),
               (0, 2, 1), (1, 1, 1), (2, 0, 1), (1, 2, 0), (2, 1, 0)]
    assert minimum(example) == source_ideal(3)
    assert singular_check(example) == (12, 66)
    print("FIDELITY explicit source example: colength=12 tangent=66")
    for k in range(1, 7):
        length = comb(k + 2, 3) + 2
        formula = comb(k + 1, 2) * comb(k + 2, 2) + 6
        r, s = comb(k + 2, 2) - 2, comb(k + 1, 2)
        predicted = {-2: 2 * s - 2, -1: r * (s - 2) + 8, 0: 2 * r}
        predicted = {d: v for d, v in predicted.items() if v}
        for name, gens in [("I", source_ideal(k)), ("J", canonical_ideal(k))]:
            dim, grades = tangent(gens)
            assert len(basis(gens)) == length
            assert dim == formula, (k, name, dim, formula)
            assert singular_check(gens) == (length, dim)
            if name == "J":
                assert grades == predicted, (k, grades, predicted)
            print(f"E2 k={k} {name}: colength={length} tangent={dim} grades={grades} dual_check=OK")
        candidates = staircase_candidates(k)
        print(f"BOREL k={k}: count={len(candidates)} extras={candidates}")
    # A fidelity diagnostic, explicitly not a refutation of the intended problem.
    # The TeX notation takes exponents of all monomials of I, but the later
    # upper-boundary clause and pictures use the compact minimal-generator hull.
    # z^(k+3) belongs to I and to the first hull, but is outside the second:
    # every minimal generator has z-coordinate <= k+2, preserved by convex sums.
    for k in range(1, 7):
        gens = source_ideal(k)
        a = (0, 0, k + 3)
        assert contains(gens, a)
        assert a[2] > max(g[2] for g in gens)
        print(f"CONVENTION_MISMATCH k={k}: z^{k+3} is an ideal monomial outside the minimal-generator hull")
    print("GLOBAL_MAXIMALITY=NOT_CHECKED")
    print("CONVEX_HULL_CLAUSES=REQUIRE_SOURCE_CONVENTION_CLARIFICATION")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
