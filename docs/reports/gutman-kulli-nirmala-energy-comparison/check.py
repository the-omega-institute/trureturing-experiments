#!/usr/bin/env python3
"""Exact checks of Gutman--Kulli Nirmala energy on graphs with at most four vertices.

Only integer, rational, and symbolic algebraic arithmetic is used. Real roots of
rational characteristic polynomials are enclosed by SymPy's exact Sturm isolation.
No numerical eigenvalue approximation is an acceptance criterion.
"""
from itertools import combinations, permutations
import sympy as sp

X = sp.Symbol("X")
EPS = sp.Rational(1, 10**8)


def require(test, message):
    if not test:
        raise AssertionError(message)


def matrices(n, edges):
    degrees = [sum(v in edge for edge in edges) for v in range(n)]
    adjacency = sp.zeros(n)
    nirmala = sp.zeros(n)
    for i, j in edges:
        adjacency[i, j] = adjacency[j, i] = 1
        nirmala[i, j] = nirmala[j, i] = sp.sqrt(degrees[i] + degrees[j])
    return degrees, adjacency, nirmala


def canonical_edges(n, edges):
    return min(tuple(sorted(tuple(sorted((p[i], p[j]))) for i, j in edges))
               for p in permutations(range(n)))


def energy_interval(poly, n):
    """Enclose sum |lambda|, with exact algebraic multiplicity."""
    poly = sp.Poly(poly, X, domain=sp.QQ)
    intervals = poly.intervals(eps=EPS)
    require(sum(m for _, m in intervals) == n, "not all eigenvalues are real")
    lower = sp.S.Zero
    upper = sp.S.Zero
    for (left, right), multiplicity in intervals:
        if right <= 0:
            lo, hi = -right, -left
        elif left >= 0:
            lo, hi = left, right
        else:
            lo, hi = sp.S.Zero, max(-left, right)
        lower += multiplicity * lo
        upper += multiplicity * hi
    return lower, upper


def sqrt_interval(c):
    """Exact rational enclosure of the positive root of X^2-c."""
    candidates = sp.Poly(X**2 - c, X, domain=sp.QQ).intervals(eps=EPS)
    positive = [(a, b) for (a, b), _ in candidates if a > 0]
    require(len(positive) == 1, "missing positive square root")
    return positive[0]


def scaled_nirmala_interval(nirmala, n):
    # On these graphs, one of the following exact rescalings gives a rational
    # characteristic polynomial. Its roots are precisely the scaled eigenvalues.
    for c in range(1, 9):
        expression = sp.expand((nirmala / sp.sqrt(c)).charpoly(X).as_expr())
        if all(a.is_Rational for a in sp.Poly(expression, X).all_coeffs()):
            lo, hi = energy_interval(expression, n)
            sl, su = sqrt_interval(c)
            return c, expression, (lo * sl, hi * su)
    raise AssertionError("no rational characteristic-polynomial rescaling")


def main():
    labeled = 0
    classes = 0
    nonempty_classes = 0
    sharp_classes = 0
    tree_classes = 0
    for n in range(5):
        pairs = list(combinations(range(n), 2))
        seen = set()
        for mask in range(1 << len(pairs)):
            labeled += 1
            edges = {pairs[k] for k in range(len(pairs)) if (mask >> k) & 1}
            canonical = canonical_edges(n, edges)
            if canonical in seen:
                continue
            seen.add(canonical)
            classes += 1
            degrees, adjacency, nirmala = matrices(n, edges)
            require(adjacency == adjacency.T and nirmala == nirmala.T,
                    "source matrices are not symmetric")
            require(all(nirmala[i, i] == 0 for i in range(n)), "nonzero diagonal")
            # Source Lemma 1: sum nu_i=0 and sum nu_i^2=2 Zg(G),
            # evaluated via trace and trace of the square.
            zagreb = sum(d**2 for d in degrees)
            require(nirmala.trace() == 0, "source Lemma 1, first identity")
            require(sp.simplify((nirmala * nirmala).trace() - 2 * zagreb) == 0,
                    "source Lemma 1, second identity")
            ea_lo, ea_hi = energy_interval(adjacency.charpoly(X).as_expr(), n)
            c, poly, (en_lo, en_hi) = scaled_nirmala_interval(nirmala, n)
            if not edges:
                require(ea_lo == ea_hi == en_lo == en_hi == 0,
                        "edgeless boundary")
                continue
            nonempty_classes += 1
            delta = min(d for d in degrees if d > 0)
            require(ea_lo > 0, "energy positivity for a graph with an edge")
            # Independent energy intervals certify the strict source comparison.
            require(en_lo > ea_hi, "expected strict comparison")
            positive_degrees = {d for d in degrees if d > 0}
            if len(positive_degrees) == 1:
                # Exact sharpness: isolates append zero rows and columns and
                # every remaining component has the same regular degree.
                require(nirmala == sp.sqrt(2 * delta) * adjacency,
                        "source regular scaling identity")
                sharp_classes += 1
                bound = "sharp exact scaling"
            else:
                require(en_lo**2 > 2 * delta * ea_hi**2,
                        "proposed sqrt(2 delta_+) lower bound")
                bound = "strict exact rational isolation"
            # Connectedness by a finite reachability expansion; n-1 edges then
            # identifies a tree, hence checks the source's Theorem 4 cases too.
            reached = {0}
            for _ in range(n):
                reached |= {v for u, v in edges if u in reached}
                reached |= {u for u, v in edges if v in reached}
            if len(reached) == n and len(edges) == n - 1:
                tree_classes += 1
            print(f"n={n} edges={canonical} delta_positive={delta} "
                  f"rescale_square={c} charpoly={sp.factor(poly)}; {bound}")
        require(len(seen) == [1, 1, 2, 4, 11][n], "unlabelled graph count")
    require((labeled, classes, nonempty_classes) == (76, 19, 14), "coverage counts")
    print(f"coverage: labeled_graphs={labeled}; isomorphism_classes={classes}; "
          f"classes_with_edges={nonempty_classes}")
    print(f"source fidelity: Lemma 1 on {classes} classes; "
          f"Theorem 4 on {tree_classes} tree classes; "
          f"regular positive-degree scaling on {sharp_classes} classes")
    print("expected settlement: strict comparison and sqrt(2 delta_+) lower bound "
          "certified for every graph on at most four vertices")
    print("boundary: all edgeless graphs in this range have both energies zero")
    print("scope: finite checks only; the universal claim needs a general proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
