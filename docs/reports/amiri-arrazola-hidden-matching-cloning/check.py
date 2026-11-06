#!/usr/bin/env python3
"""Exact checks for Amiri--Arrazola arXiv:1610.06345v1, Eqs. (9), (20)--(21).

Only integers, fractions and exact SymPy polynomials are evidence. Finite checks
support, but do not prove, the proposed all-dimensions spectral argument.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product, combinations
import sympy as sp

T = sp.Symbol('T')


def signs(n):
    return product((-1, 1), repeat=n)


def moment_table(n):
    """Compute the fourth sign moment by summing the entire source ensemble."""
    table = {}
    all_signs = list(signs(n))
    for q in product(range(n), repeat=4):
        value = F(sum(s[q[0]] * s[q[1]] * s[q[2]] * s[q[3]] for s in all_signs), 2**n)
        a, b, c, d = q
        pairing = int(a == c and b == d) + int(a == d and b == c) + int(a == b and c == d) - 2 * int(a == b == c == d)
        assert value == pairing, (n, q, value, pairing)
        table[q] = int(value)
    return table


def parity(q):
    return tuple(sorted(a for a, count in Counter(q).items() if count % 2))


def operator(n, table=None):
    """Integer L = 2 n^2 Q. For n<=8, entries use enumerated sign moments."""
    def m(a, b, c, d):
        if table is not None:
            return table[a, b, c, d]
        return int(all(v % 2 == 0 for v in Counter((a, b, c, d)).values()))
    basis = list(product(range(n), repeat=3))
    rows = {q: {} for q in basis}
    # Either third output coordinate or second output coordinate agrees.
    for a, b, c in basis:
        row = rows[a, b, c]
        for d, e in product(range(n), repeat=2):
            val = m(a, b, d, e)
            if val:
                row[d, e, c] = row.get((d, e, c), 0) + val
        for d, f in product(range(n), repeat=2):
            val = m(a, c, d, f)
            if val:
                row[d, b, f] = row.get((d, b, f), 0) + val
    return rows


def matrix(rows, basis):
    return sp.Matrix([[rows[a].get(b, 0) for b in basis] for a in basis])


def matching_factorization(n):
    ring = list(range(n))
    out = []
    for _ in range(n - 1):
        out.append([(ring[i], ring[n - 1 - i]) for i in range(n // 2)])
        ring = [ring[0], ring[-1], *ring[1:-1]]
    assert Counter(tuple(sorted(e)) for match in out for e in match) == Counter(combinations(range(n), 2))
    return out


def source_error_identity(n):
    """Reproduce Eq. (9) with its n-1 disjoint matchings and all input strings."""
    matches = matching_factorization(n)
    for s in signs(n):
        # A nontrivial exact density: (2/3) rho_x + (1/3) |0><0|.
        eta = [[F(2 * s[i] * s[j], 3 * n) + F(int(i == j == 0), 3) for j in range(n)] for i in range(n)]
        fidelity = sum(F(s[i] * s[j], n) * eta[i][j] for i, j in product(range(n), repeat=2))
        errors = [sum((eta[i][i] + eta[j][j] - s[i]*s[j]*(eta[i][j]+eta[j][i])) / 2 for i, j in match) for match in matches]
        average_error = sum(errors, F(0)) / (n - 1)
        assert average_error == F(n, 2 * (n - 1)) * (1 - fidelity)
    print(f'SOURCE_EQ9 n={n} strings={2**n} fidelity={fidelity} error={average_error} OK')


def check_spectrum(n, rows):
    sectors = defaultdict(list)
    for q in rows:
        sectors[parity(q)].append(q)
    assert len(sectors) == n + sp.binomial(n, 3)
    assert sum(map(len, sectors.values())) == n**3
    for a, row in rows.items():
        for b, val in row.items():
            assert rows[b].get(a, 0) == val
            assert parity(a) == parity(b)
    singleton_poly = T * (T-2) * (T-n) * (T-n-2) * (T*(T-1)*(T-3))**(n-2)
    # Every singleton block is checked, not only a representative.
    for k in range(n):
        block = sectors[(k,)]
        assert len(block) == 3*n-2
        assert sp.expand(matrix(rows, block).charpoly(T).as_expr() - singleton_poly) == 0
    distinct_poly = T * (T-1)**2 * (T-3)**2 * (T-4)
    for key, block in sectors.items():
        if len(key) == 3:
            assert len(block) == 6
            assert sp.expand(matrix(rows, block).charpoly(T).as_expr() - distinct_poly) == 0
    # All roots listed above are real and between 0 and n+2 (n>=2).
    # Hence the exact symmetric operator Q has norm (n+2)/(2 n^2).
    norm_q = F(n+2, 2*n*n)
    assert n * norm_q == F(1, 2) + F(1, n)
    return norm_q


def top_vectors(n):
    ws = []
    for k in range(n):
        w = {}
        for i in range(n):
            if i != k:
                w[i, i, k] = F(1)
                w[i, k, i] = F(1)
                w[k, i, i] = F(2, n)
        w[k, k, k] = F(2*(n-1), n)
        ws.append(w)
    return ws


def check_attainment(n, rows):
    ws = top_vectors(n)
    D = F(2*(n-1)*(n+2), n)
    for k, w in enumerate(ws):
        assert sum(v*v for v in w.values()) == D
        for a, row in rows.items():
            assert sum(v*w.get(b, F(0)) for b, v in row.items()) == (n+2)*w.get(a, F(0))
        for l in range(k):
            assert sum(v*ws[l].get(a, F(0)) for a, v in w.items()) == 0
        tr = [[sum(w.get((a,y,z), F(0))*w.get((b,y,z), F(0)) for y,z in product(range(n),repeat=2))/D for b in range(n)] for a in range(n)]
        for a,b in product(range(n),repeat=2):
            expected = (F(2,n+2) if a==k else F(n,(n-1)*(n+2))) if a==b else 0
            assert tr[a][b] == expected
    # J=sum |w_k><w_k|/D is positive by its displayed Gram decomposition.
    # The partial trace is exactly I, so it specifies a CPTP map.
    for a,b in product(range(n),repeat=2):
        tr_j = sum(sum(w.get((a,y,z),F(0))*w.get((b,y,z),F(0)) for y,z in product(range(n),repeat=2))/D for w in ws)
        assert tr_j == int(a==b)
    if n <= 6:
        for s in signs(n):
            # For each Kraus index, reshape w as a map from input a to Y,Z.
            out_vectors = [[[sum(w.get((a,y,z),F(0))*s[a] for a in range(n)) for z in range(n)] for y in range(n)] for w in ws]
            fy = sum(sum(sum(F(s[y])*v[y][z] for y in range(n))**2 for z in range(n)) for v in out_vectors)/(D*n*n)
            fz = sum(sum(sum(F(s[z])*v[y][z] for z in range(n))**2 for y in range(n)) for v in out_vectors)/(D*n*n)
            assert fy == fz == F(1,2)+F(1,n)
    return D


def symbolic_blocks():
    n = sp.Symbol('n', integer=True, positive=True)
    zero_sum = sp.Matrix([[1,0,1],[0,1,1],[1,1,2]])
    uniform = sp.Matrix([[n,0,1,1],[0,n,1,1],[1,1,2,0],[n-1,n-1,0,2]])
    assert sp.factor(zero_sum.charpoly(T).as_expr()) == T*(T-1)*(T-3)
    assert sp.expand(uniform.charpoly(T).as_expr()-T*(T-2)*(T-n)*(T-n-2)) == 0
    v = sp.Matrix([1,1,sp.Rational(2)/n,2*(n-1)/n])
    assert all(sp.simplify(x)==0 for x in uniform*v-(n+2)*v)
    print('SYMBOLIC_BLOCKS zero_sum={0,1,3} uniform={0,2,n,n+2} top_vector OK')


def main():
    print('EXACT_ARITHMETIC Python fractions; SymPy integer/rational matrices and polynomials')
    print('SOURCE arXiv:1610.06345v1 Eqs. (9), (20), (21)')
    symbolic_blocks()
    for n in (4, 6):
        source_error_identity(n)
    for n in range(2, 15):
        table = moment_table(n) if n <= 8 else None
        rows = operator(n, table)
        norm_q = check_spectrum(n, rows)
        norm_sq = check_attainment(n, rows)
        print(f'n={n} dim={n**3} fourth_moment={"enumerated" if table else "parity"} norm_Q={norm_q} optimum={n*norm_q} top_norm_squared={norm_sq} CPTP_partial_trace=I OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
