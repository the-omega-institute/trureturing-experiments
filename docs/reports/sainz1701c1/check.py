"""Exact checks for Sainz–Wolfe Conjecture 1, arXiv:1701.05171v2."""

from itertools import combinations, product

from sympy import Matrix, Rational, diag, eye, kronecker_product, simplify, sqrt, symbols, zeros

BITS = (0, 1)


def require(condition, explanation):
    if not condition:
        raise AssertionError(explanation)


def incidence(edge, size):
    return Matrix([[int(vertex in edge) for vertex in range(size)] + [1]])


def fr_contexts():
    alice_first = {
        frozenset((a, b, x, y_of_a[a]) for a, b in product(BITS, repeat=2))
        for x, y_of_a in product(BITS, product(BITS, repeat=2))
    }
    bob_first = {
        frozenset((a, b, x_of_b[b], y) for a, b in product(BITS, repeat=2))
        for y, x_of_b in product(BITS, product(BITS, repeat=2))
    }
    return alice_first | bob_first


def check_source():
    # Appendix E incidence tables, in printed column order and zero-based indexing.
    original = (
        (0, 1, 2, 3), (0, 1, 6, 7), (0, 2, 9, 11), (1, 3, 8, 10),
        (2, 3, 4, 5), (4, 6, 13, 15), (5, 7, 12, 14), (8, 9, 14, 15),
    )
    extra = ((12, 13, 14, 15), (10, 11, 12, 13), (8, 9, 10, 11), (4, 5, 6, 7))
    columns = [(a, b, x, y) for x, y, a, b in product(BITS, repeat=4)]
    require({frozenset(columns[i] for i in e) for e in original + extra} == fr_contexts(),
            "Appendix E's H₂ is the full bipartite Foulis–Randall product")
    rows = [incidence(e, 16) for e in original]
    new_rows = [incidence(e, 16) for e in extra]
    old_rref, old_pivots = Matrix.vstack(*rows).rref()
    new_rref, new_pivots = Matrix.vstack(*rows, *new_rows).rref()
    rank = len(old_pivots)
    require(old_pivots == new_pivots and old_rref[:rank, :] == new_rref[:rank, :],
            "Appendix E's augmented row spaces coincide")
    for new_row in new_rows:
        require(any(sum((sign * row for sign, row in zip(signs, rows)), zeros(1, 17))
                    == new_row for signs in product((-1, 0, 1), repeat=len(rows))
                    if sum(signs) == 1), "Appendix E signed row combination")
    print(f"FIDELITY Appendix E: 8 and 12 contexts have the same augmented RREF (rank {rank})")


def graph_edges(contexts):
    return {frozenset(pair) for edge in contexts for pair in combinations(edge, 2)}


def matrix_sum(matrices, size):
    return sum(matrices, zeros(size))


def refine(rho, projections, c, d, alpha, beta):
    size = rho.rows
    p_c = matrix_sum((projections[v] for v in c), size)
    p_d = matrix_sum((projections[v] for v in d), size)
    na, nb = len(alpha), len(beta)
    extended = {v: kronecker_product(value, eye(na), eye(nb))
                for v, value in projections.items()}
    for label, count, prefix in (("A", na, p_d), ("B", nb, p_c)):
        for i in range(count):
            coordinate = eye(count)[:, i]
            basis_projection = coordinate * coordinate.T
            a_factor = basis_projection if label == "A" else eye(na)
            b_factor = basis_projection if label == "B" else eye(nb)
            extended[(label, i)] = kronecker_product(prefix, a_factor, b_factor)
    return kronecker_product(rho, diag(*alpha), diag(*beta)), extended


def check_small_family():
    cases = 0
    t = symbols("t", real=True)
    for sizes in product(range(1, 4), repeat=4):
        na, nb, nc, nd = sizes
        a = set(range(na))
        b = set(range(na, na + nb))
        c = set(range(na + nb, na + nb + nc))
        d = set(range(na + nb + nc, sum(sizes)))
        original = [a | c, b | d, c | d]
        new = a | b
        require(incidence(new, sum(sizes)) == incidence(original[0], sum(sizes))
                + incidence(original[1], sum(sizes)) - incidence(original[2], sum(sizes)),
                "augmented incidence identity")
        require(graph_edges([*original, new]) - graph_edges(original)
                == {frozenset(pair) for pair in product(a, b)}, "exactly the cross edges")
        # A spare direction gives a nonzero context slack. The state has no
        # support there, as guaranteed by normalization. Include both zero masses.
        core = {v: eye(3)[:, int(v in d)] * eye(3)[:, int(v in d)].T
                for v in (min(c), min(d))}
        core.update({v: zeros(3) for v in c | d if v not in core})
        alpha = [Rational(2 * (i + 1), na * (na + 1)) for i in range(na)]
        beta = [Rational(2 * (i + 1), nb * (nb + 1)) for i in range(nb)]
        rho, projections = refine(diag(1 - t, t, 0), core, c, d, alpha, beta)
        fresh_a = {("A", i) for i in range(na)}
        fresh_b = {("B", i) for i in range(nb)}
        size = rho.rows
        target = kronecker_product(diag(1, 1, 0), eye(na), eye(nb))
        for e in (c | d, c | fresh_a, d | fresh_b, fresh_a | fresh_b):
            require(matrix_sum((projections[v] for v in e), size) == target,
                    "all context sums have the prescribed positive slack")
        require(simplify(rho.trace()) == 1, "trace one")
        for (label, weights, mass) in (("A", alpha, t), ("B", beta, 1 - t)):
            for i, weight in enumerate(weights):
                p = projections[(label, i)]
                require(p.H == p and p * p == p, "fresh orthogonal projection")
                born = (rho * p).trace()
                require(simplify(born - mass * weight) == 0, "symbolic Born mass")
                for endpoint in BITS:
                    require(born.subs(t, endpoint) == (mass * weight).subs(t, endpoint),
                            "endpoint Born mass")
        cases += 1
    require(cases == 3**4, "all requested size cases executed")
    print(f"FAMILY: {cases} block-size tuples with symbolic masses and both endpoints")


def check_quantum_witness():
    z = diag(1, -1)
    x = Matrix([[0, 1], [1, 0]])
    alice = (z, x)
    bob = ((z + x) / sqrt(2), (z - x) / sqrt(2))
    phi = Matrix([1, 0, 0, 1]) / sqrt(2)
    rho = Rational(3, 4) * phi * phi.H + eye(4) / 16
    require(rho.eigenvals() == {Rational(13, 16): 1, Rational(1, 16): 3},
            "positive full-rank core density")
    core = {(a, b, i, j): kronecker_product((eye(2) + (-1)**a * alice[i]) / 2,
                                            (eye(2) + (-1)**b * bob[j]) / 2)
            for a, b, i, j in product(BITS, repeat=4)}
    c = {(0, b, 0, 0) for b in BITS}
    d = {(1, b, 0, 0) for b in BITS}
    a = {("A", i) for i in BITS}
    b = {("B", i) for i in BITS}
    base = fr_contexts()
    old = base | {frozenset(c | a), frozenset(d | b)}
    new = old | {frozenset(a | b)}
    require(graph_edges(new) - graph_edges(old) == {frozenset(pair) for pair in product(a, b)},
            "four new orthogonality edges")
    reachable = {next(iter(core))}
    while True:
        joined = reachable | set().union(*(e for e in old if e & reachable))
        if joined == reachable:
            break
        reachable = joined
    require(reachable == set(core) | a | b, "connected orthogonality graph")
    uniform = (Rational(1, 2),) * 2
    new_rho, projections = refine(rho, core, c, d, uniform, uniform)
    require(new_rho.trace() == 1 and new_rho.rows == 16, "extended density and dimension")
    require(all(p != zeros(16) and p.H == p and (p * p).applyfunc(simplify) == p
                for p in projections.values()), "nonzero orthogonal projectors")
    require(len({tuple(p) for p in projections.values()}) == len(projections) == 20,
            "20 distinct projectors")
    require(all(matrix_sum((projections[v] for v in e), 16).applyfunc(simplify) == eye(16)
                for e in new), "all context sums equal identity")
    p = {v: simplify((new_rho * projector).trace()) for v, projector in projections.items()}
    require(all(value.is_positive for value in p.values()), "all 20 probabilities positive")
    require(all(p[v] == Rational(1, 4) for v in a | b), "four fresh probabilities are 1/4")
    for a_out, b_out, i, j in core:
        expected = (1 + (-1)**(a_out + b_out + i * j) * 3 / (4 * sqrt(2))) / 4
        require(simplify(p[(a_out, b_out, i, j)] - expected) == 0, "source-model Born formula")

    def chsh(behavior):
        return simplify(sum((-1)**(a_out + b_out + i * j) * behavior.get((a_out, b_out, i, j), 0)
                            for a_out, b_out, i, j in core))

    deterministic = []
    for outputs in product(tuple(product(BITS, repeat=2)), repeat=4):
        selected = {(*out, i, j) for out, (i, j) in zip(outputs, product(BITS, repeat=2))}
        if all(len(selected & e) == 1 for e in base):
            deterministic.append({v: 1 for v in selected})
    require(len(deterministic) == 16 and {chsh(model) for model in deterministic} == {-2, 2},
            "exact CHSH bound for all deterministic FR models")
    score = chsh(p)
    require(score == 3 * sqrt(2) / 2 and (score - 2).is_positive, "strict CHSH violation")
    print(f"WITNESS: 20 positive probabilities, 4 added edges, Hilbert dimension 16, CHSH {score} > 2")


if __name__ == "__main__":
    check_source()
    check_small_family()
    check_quantum_witness()
    print("ALL_CHECKS_OK")
